// signature.cpp -- Rcpp kernels for truncated path signatures.
//
// Three exported functions share a common workspace structure:
//   sig_terminal_cpp(path, N)  -- terminal signature of one path
//   sig_batch_cpp(paths, N)    -- terminal signatures of many paths
//   sig_path_cpp(path, N)      -- signature at every time step
//
// All three use the exact piecewise-linear formulation
//   S(X) = exp(delta_1) (x) exp(delta_2) (x) ... (x) exp(delta_(T-1))
// with segment-wise tensor exponentials and Chen products. The (d, N)
// dependent precompute (parent positions, split table, letter buffer,
// reciprocals) lives in SigWorkspace so batch and path variants pay
// it once.

#include <Rcpp.h>
#include <climits>
#include <vector>
using namespace Rcpp;

// =====================================================================
// Workspace -- depends only on (d, N), not on path data.
// =====================================================================

struct SigWorkspace {
  int d;
  int N;
  int nWords;
  std::vector<int>    parentPos;
  std::vector<int>    lastLetter;
  std::vector<int>    wordLen;
  std::vector<int>    splitOff;
  std::vector<int>    splitPre;
  std::vector<int>    splitSuf;
  std::vector<double> invK;
};

static SigWorkspace make_workspace(int d, int N) {
  SigWorkspace ws;
  ws.d = d; ws.N = N;

  long long nWordsLL = 0;
  long long dk = 1;
  for (int k = 0; k <= N; ++k) {
    nWordsLL += dk;
    if (nWordsLL > (long long) INT_MAX) {
      Rcpp::stop("signature size overflows int (d=%d, N=%d)", d, N);
    }
    dk *= d;
  }
  ws.nWords = (int) nWordsLL;
  if (N == 0) return ws;

  ws.parentPos.resize(ws.nWords);
  ws.lastLetter.resize(ws.nWords);
  ws.wordLen.resize(ws.nWords);
  ws.parentPos[0] = -1; ws.lastLetter[0] = -1; ws.wordLen[0] = 0;
  for (int w = 1; w < ws.nWords; ++w) {
    ws.parentPos[w]  = (w - 1) / d;
    ws.lastLetter[w] = (w - 1) % d;
    ws.wordLen[w]    = ws.wordLen[ws.parentPos[w]] + 1;
  }

  std::vector<int> letterOff(ws.nWords + 1);
  letterOff[0] = 0;
  for (int w = 0; w < ws.nWords; ++w)
    letterOff[w + 1] = letterOff[w] + ws.wordLen[w];
  std::vector<int> letters(letterOff[ws.nWords]);
  for (int w = 1; w < ws.nWords; ++w) {
    const int p = ws.parentPos[w];
    const int k = ws.wordLen[w];
    for (int i = 0; i < k - 1; ++i)
      letters[letterOff[w] + i] = letters[letterOff[p] + i];
    letters[letterOff[w] + k - 1] = ws.lastLetter[w];
  }

  auto wordPos = [d](const int* L, int k) -> int {
    if (k == 0) return 0;
    long long cumLevel;
    if (d == 1) {
      cumLevel = k;
    } else {
      long long pk = 1;
      for (int i = 0; i < k; ++i) pk *= d;
      cumLevel = (pk - 1) / (d - 1);
    }
    long long offset = 0, power = 1;
    for (int m = k - 1; m >= 0; --m) {
      offset += L[m] * power;
      power  *= d;
    }
    return (int)(cumLevel + offset);
  };

  ws.splitOff.resize(ws.nWords + 1);
  ws.splitOff[0] = 0;
  for (int w = 0; w < ws.nWords; ++w)
    ws.splitOff[w + 1] = ws.splitOff[w] + ws.wordLen[w] + 1;
  const int totalSplits = ws.splitOff[ws.nWords];
  ws.splitPre.resize(totalSplits);
  ws.splitSuf.resize(totalSplits);
  for (int w = 0; w < ws.nWords; ++w) {
    const int k    = ws.wordLen[w];
    const int* L   = &letters[letterOff[w]];
    const int base = ws.splitOff[w];
    for (int j = 0; j <= k; ++j) {
      ws.splitPre[base + j] = wordPos(L,     j);
      ws.splitSuf[base + j] = wordPos(L + j, k - j);
    }
  }

  ws.invK.assign(N + 1, 0.0);
  for (int k = 1; k <= N; ++k) ws.invK[k] = 1.0 / k;
  return ws;
}

// =====================================================================
// Inner kernels.
// =====================================================================

static inline void compute_segment_exp(
    const NumericMatrix& path, int t,
    const SigWorkspace& ws,
    std::vector<double>& E)
{
  E[0] = 1.0;
  for (int w = 1; w < ws.nWords; ++w) {
    const int l     = ws.lastLetter[w];
    const double dl = path(t + 1, l) - path(t, l);
    E[w] = E[ws.parentPos[w]] * dl * ws.invK[ws.wordLen[w]];
  }
}

static inline void chen_product(
    const std::vector<double>& S,
    const std::vector<double>& E,
    const SigWorkspace& ws,
    std::vector<double>& newS)
{
  for (int w = 0; w < ws.nWords; ++w) {
    const int s = ws.splitOff[w], e = ws.splitOff[w + 1];
    double acc = 0.0;
    for (int i = s; i < e; ++i) {
      acc += S[ws.splitPre[i]] * E[ws.splitSuf[i]];
    }
    newS[w] = acc;
  }
}

static void run_path_terminal(
    const NumericMatrix& path,
    const SigWorkspace& ws,
    std::vector<double>& S,
    std::vector<double>& E,
    std::vector<double>& newS)
{
  std::fill(S.begin(), S.end(), 0.0);
  S[0] = 1.0;
  const int T = path.nrow();
  for (int t = 0; t < T - 1; ++t) {
    compute_segment_exp(path, t, ws, E);
    chen_product(S, E, ws, newS);
    S.swap(newS);
  }
}

// =====================================================================
// Exported: sig_terminal_cpp -- one path, terminal signature.
// =====================================================================

// [[Rcpp::export]]
NumericVector sig_terminal_cpp(const NumericMatrix& path, int N) {
  const int T = path.nrow();
  const int d = path.ncol();
  if (T < 1) Rcpp::stop("path must have at least one row");
  if (d < 1) Rcpp::stop("path must have at least one column");
  if (N < 0) Rcpp::stop("depth must be non-negative");

  auto ws = make_workspace(d, N);
  NumericVector out(ws.nWords);
  out[0] = 1.0;
  if (N == 0 || T < 2) return out;

  std::vector<double> S(ws.nWords), E(ws.nWords), newS(ws.nWords);
  run_path_terminal(path, ws, S, E, newS);
  for (int w = 0; w < ws.nWords; ++w) out[w] = S[w];
  return out;
}

// =====================================================================
// Exported: sig_batch_cpp -- many paths sharing (d, N).
//
// Paths come in as a List of NumericMatrix; all must have the same d
// (R wrapper validates). T may vary. Output is nPaths x nWords; row i
// is the terminal signature of paths[[i]].
// =====================================================================

// [[Rcpp::export]]
NumericMatrix sig_batch_cpp(const List& paths, int N) {
  const int nPaths = paths.size();
  if (nPaths < 1) Rcpp::stop("paths must contain at least one path");

  NumericMatrix first = paths[0];
  const int d = first.ncol();
  if (d < 1) Rcpp::stop("paths must have at least one column");
  if (N < 0) Rcpp::stop("depth must be non-negative");

  auto ws = make_workspace(d, N);
  NumericMatrix out(nPaths, ws.nWords);
  for (int i = 0; i < nPaths; ++i) out(i, 0) = 1.0;
  if (N == 0) return out;

  std::vector<double> S(ws.nWords), E(ws.nWords), newS(ws.nWords);
  for (int i = 0; i < nPaths; ++i) {
    NumericMatrix path = paths[i];
    if (path.ncol() != d) {
      Rcpp::stop("paths[[%d]] has ncol = %d, expected %d",
                 i + 1, path.ncol(), d);
    }
    if (path.nrow() < 2) continue;
    run_path_terminal(path, ws, S, E, newS);
    for (int w = 0; w < ws.nWords; ++w) out(i, w) = S[w];
  }
  return out;
}

// =====================================================================
// Exported: sig_path_cpp -- signature at every time step.
//
// Returns a T x nWords matrix. Row t is the signature of the path
// restricted to [1, t+1]. Row 0 is the identity; row T-1 is terminal.
// =====================================================================

// [[Rcpp::export]]
NumericMatrix sig_path_cpp(const NumericMatrix& path, int N) {
  const int T = path.nrow();
  const int d = path.ncol();
  if (T < 1) Rcpp::stop("path must have at least one row");
  if (d < 1) Rcpp::stop("path must have at least one column");
  if (N < 0) Rcpp::stop("depth must be non-negative");

  auto ws = make_workspace(d, N);
  NumericMatrix out(T, ws.nWords);
  for (int t = 0; t < T; ++t) out(t, 0) = 1.0;
  if (N == 0 || T < 2) return out;

  std::vector<double> S(ws.nWords, 0.0), E(ws.nWords, 0.0), newS(ws.nWords, 0.0);
  S[0] = 1.0;
  for (int t = 0; t < T - 1; ++t) {
    compute_segment_exp(path, t, ws, E);
    chen_product(S, E, ws, newS);
    S.swap(newS);
    for (int w = 0; w < ws.nWords; ++w) out(t + 1, w) = S[w];
  }
  return out;
}

// =====================================================================
// Workspace handle for callers that want to amortize make_workspace cost.
// Returned as an Rcpp::XPtr that R can cache; finalizer frees on GC.
// =====================================================================

// [[Rcpp::export]]
SEXP build_sig_workspace(int d, int N) {
  if (d < 1) Rcpp::stop("d must be at least 1");
  if (N < 0) Rcpp::stop("N must be non-negative");
  SigWorkspace* ws = new SigWorkspace(make_workspace(d, N));
  return Rcpp::XPtr<SigWorkspace>(ws, true);
}

// =====================================================================
// Exported: tensor_product_cpp -- truncated tensor product.
//
// Math: (a (x) b)[w] = sum over splits w = u.v of a[u] * b[v].
// This is the same sum chen_product already computes; we expose it as
// a standalone primitive over arbitrary inputs.
// =====================================================================

// [[Rcpp::export]]
NumericVector tensor_product_cpp(const NumericVector& a,
                                 const NumericVector& b,
                                 SEXP ws_xptr) {
  Rcpp::XPtr<SigWorkspace> ws(ws_xptr);
  const int nWords = ws->nWords;
  if ((int) a.size() != nWords) {
    Rcpp::stop("`a` length %d does not match workspace nWords %d",
               (int) a.size(), nWords);
  }
  if ((int) b.size() != nWords) {
    Rcpp::stop("`b` length %d does not match workspace nWords %d",
               (int) b.size(), nWords);
  }
  std::vector<double> aVec(a.begin(), a.end());
  std::vector<double> bVec(b.begin(), b.end());
  std::vector<double> out(nWords);
  chen_product(aVec, bVec, *ws, out);
  return NumericVector(out.begin(), out.end());
}

// =====================================================================
// Shuffle product: separate workspace, its own cache.
// =====================================================================

struct ShuffleWorkspace {
  int d;
  int N;
  int nWords;
  std::vector<int> wordLen;
  std::vector<int> shuffleOff;
  std::vector<int> shuffleUPos;
  std::vector<int> shuffleVPos;
};

static ShuffleWorkspace make_shuffle_workspace(int d, int N) {
  ShuffleWorkspace ws;
  ws.d = d; ws.N = N;

  long long nWordsLL = 0, dk = 1;
  for (int k = 0; k <= N; ++k) {
    nWordsLL += dk;
    if (nWordsLL > (long long) INT_MAX) {
      Rcpp::stop("signature size overflows int (d=%d, N=%d)", d, N);
    }
    dk *= d;
  }
  ws.nWords = (int) nWordsLL;
  ws.wordLen.assign(ws.nWords, 0);

  if (N == 0) {
    ws.shuffleOff  = {0, 1};
    ws.shuffleUPos = {0};
    ws.shuffleVPos = {0};
    return ws;
  }

  std::vector<int> parentPos(ws.nWords), lastLetter(ws.nWords);
  parentPos[0] = -1; lastLetter[0] = -1;
  for (int w = 1; w < ws.nWords; ++w) {
    parentPos[w]  = (w - 1) / d;
    lastLetter[w] = (w - 1) % d;
    ws.wordLen[w] = ws.wordLen[parentPos[w]] + 1;
  }

  std::vector<int> letterOff(ws.nWords + 1);
  letterOff[0] = 0;
  for (int w = 0; w < ws.nWords; ++w)
    letterOff[w + 1] = letterOff[w] + ws.wordLen[w];
  std::vector<int> letters(letterOff[ws.nWords]);
  for (int w = 1; w < ws.nWords; ++w) {
    const int p = parentPos[w];
    const int k = ws.wordLen[w];
    for (int i = 0; i < k - 1; ++i)
      letters[letterOff[w] + i] = letters[letterOff[p] + i];
    letters[letterOff[w] + k - 1] = lastLetter[w];
  }

  auto wordPos = [d](const int* L, int k) -> int {
    if (k == 0) return 0;
    long long cumLevel;
    if (d == 1) {
      cumLevel = k;
    } else {
      long long pk = 1;
      for (int i = 0; i < k; ++i) pk *= d;
      cumLevel = (pk - 1) / (d - 1);
    }
    long long offset = 0, power = 1;
    for (int m = k - 1; m >= 0; --m) {
      offset += L[m] * power;
      power  *= d;
    }
    return (int)(cumLevel + offset);
  };

  ws.shuffleOff.resize(ws.nWords + 1);
  ws.shuffleOff[0] = 0;
  for (int w = 0; w < ws.nWords; ++w)
    ws.shuffleOff[w + 1] = ws.shuffleOff[w] + (1 << ws.wordLen[w]);
  const int total = ws.shuffleOff[ws.nWords];
  ws.shuffleUPos.resize(total);
  ws.shuffleVPos.resize(total);

  std::vector<int> uLet, vLet;
  for (int w = 0; w < ws.nWords; ++w) {
    const int k = ws.wordLen[w];
    const int* L = &letters[letterOff[w]];
    const int base = ws.shuffleOff[w];
    const int numSubsets = 1 << k;
    for (int mask = 0; mask < numSubsets; ++mask) {
      uLet.clear(); vLet.clear();
      for (int i = 0; i < k; ++i) {
        if (mask & (1 << i)) uLet.push_back(L[i]);
        else                 vLet.push_back(L[i]);
      }
      ws.shuffleUPos[base + mask] = wordPos(uLet.data(), (int) uLet.size());
      ws.shuffleVPos[base + mask] = wordPos(vLet.data(), (int) vLet.size());
    }
  }
  return ws;
}

static inline void shuffle_product(
    const std::vector<double>& a,
    const std::vector<double>& b,
    const ShuffleWorkspace& ws,
    std::vector<double>& out)
{
  for (int w = 0; w < ws.nWords; ++w) {
    const int s = ws.shuffleOff[w], e = ws.shuffleOff[w + 1];
    double acc = 0.0;
    for (int i = s; i < e; ++i) {
      acc += a[ws.shuffleUPos[i]] * b[ws.shuffleVPos[i]];
    }
    out[w] = acc;
  }
}

// [[Rcpp::export]]
SEXP build_shuffle_workspace(int d, int N) {
  if (d < 1) Rcpp::stop("d must be at least 1");
  if (N < 0) Rcpp::stop("N must be non-negative");
  ShuffleWorkspace* ws = new ShuffleWorkspace(make_shuffle_workspace(d, N));
  return Rcpp::XPtr<ShuffleWorkspace>(ws, true);
}

// [[Rcpp::export]]
NumericVector shuffle_product_cpp(const NumericVector& a,
                                  const NumericVector& b,
                                  SEXP ws_xptr) {
  Rcpp::XPtr<ShuffleWorkspace> ws(ws_xptr);
  const int nWords = ws->nWords;
  if ((int) a.size() != nWords) {
    Rcpp::stop("`a` length %d does not match workspace nWords %d",
               (int) a.size(), nWords);
  }
  if ((int) b.size() != nWords) {
    Rcpp::stop("`b` length %d does not match workspace nWords %d",
               (int) b.size(), nWords);
  }
  std::vector<double> aVec(a.begin(), a.end());
  std::vector<double> bVec(b.begin(), b.end());
  std::vector<double> out(nWords);
  shuffle_product(aVec, bVec, *ws, out);
  return NumericVector(out.begin(), out.end());
}
