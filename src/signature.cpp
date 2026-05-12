// signature.cpp -- Rcpp kernel for truncated path signatures.
//
// Implements the exact piecewise-linear signature via segment-wise tensor
// exponentials and Chen products:
//
//   S(X) = exp(delta_1) (x) exp(delta_2) (x) ... (x) exp(delta_{T-1})
//
// where each delta_i = X[i+1] - X[i] is a path increment in R^d, exp is
// the truncated tensor exponential, and (x) is the truncated tensor
// product on words up to length N.
//
// Bitwise exact for piecewise-linear paths -- matches R/reference.R to
// machine precision.

#include <Rcpp.h>
#include <climits>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sig_terminal_cpp(const NumericMatrix& path, int N) {
  const int T = path.nrow();
  const int d = path.ncol();

  if (T < 1) Rcpp::stop("path must have at least one row");
  if (d < 1) Rcpp::stop("path must have at least one column");
  if (N < 0) Rcpp::stop("depth must be non-negative");

  // nWords = sum_{k=0..N} d^k, with overflow guard.
  long long nWordsLL = 0;
  long long dk = 1;
  for (int k = 0; k <= N; ++k) {
    nWordsLL += dk;
    if (nWordsLL > (long long) INT_MAX) {
      Rcpp::stop("signature size overflows int (d=%d, N=%d)", d, N);
    }
    dk *= d;
  }
  const int nWords = (int) nWordsLL;

  NumericVector out(nWords);
  out[0] = 1.0;
  if (N == 0 || T < 2) return out;

  // --- word metadata: parent, last letter, length -----------------
  std::vector<int> parentPos(nWords);
  std::vector<int> lastLetter(nWords);
  std::vector<int> wordLen(nWords);
  parentPos[0] = -1; lastLetter[0] = -1; wordLen[0] = 0;
  for (int w = 1; w < nWords; ++w) {
    parentPos[w]  = (w - 1) / d;
    lastLetter[w] = (w - 1) % d;
    wordLen[w]    = wordLen[parentPos[w]] + 1;
  }

  // --- reconstruct letters into a flat buffer ---------------------
  std::vector<int> letterOff(nWords + 1);
  letterOff[0] = 0;
  for (int w = 0; w < nWords; ++w)
    letterOff[w + 1] = letterOff[w] + wordLen[w];
  std::vector<int> letters(letterOff[nWords]);
  for (int w = 1; w < nWords; ++w) {
    const int p = parentPos[w];
    const int k = wordLen[w];
    for (int i = 0; i < k - 1; ++i)
      letters[letterOff[w] + i] = letters[letterOff[p] + i];
    letters[letterOff[w] + k - 1] = lastLetter[w];
  }

  // --- position-of-word helper ------------------------------------
  // pos = (d^k - 1)/(d - 1) + sum_{m=0..k-1} letter[m] * d^{k-1-m}
  // (closed form (d^k - 1)/(d - 1) is replaced by k when d == 1)
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

  // --- split table: for each w, all (prefix, suffix) positions ----
  std::vector<int> splitOff(nWords + 1);
  splitOff[0] = 0;
  for (int w = 0; w < nWords; ++w)
    splitOff[w + 1] = splitOff[w] + wordLen[w] + 1;
  const int totalSplits = splitOff[nWords];
  std::vector<int> splitPre(totalSplits), splitSuf(totalSplits);
  for (int w = 0; w < nWords; ++w) {
    const int k    = wordLen[w];
    const int* L   = &letters[letterOff[w]];
    const int base = splitOff[w];
    for (int j = 0; j <= k; ++j) {
      splitPre[base + j] = wordPos(L,     j);
      splitSuf[base + j] = wordPos(L + j, k - j);
    }
  }

  // --- precomputed 1/k for the segment exponential ----------------
  std::vector<double> invK(N + 1, 0.0);
  for (int k = 1; k <= N; ++k) invK[k] = 1.0 / k;

  // --- main loop: S := S (x) exp(delta_t) -------------------------
  std::vector<double> S(nWords, 0.0), E(nWords, 0.0), newS(nWords, 0.0);
  S[0] = 1.0;
  for (int t = 0; t < T - 1; ++t) {
    // Segment exponential by parent recurrence.
    E[0] = 1.0;
    for (int w = 1; w < nWords; ++w) {
      const int l     = lastLetter[w];
      const double dl = path(t + 1, l) - path(t, l);
      E[w] = E[parentPos[w]] * dl * invK[wordLen[w]];
    }
    // Chen product.
    for (int w = 0; w < nWords; ++w) {
      const int s = splitOff[w], e = splitOff[w + 1];
      double acc = 0.0;
      for (int i = s; i < e; ++i) {
        acc += S[splitPre[i]] * E[splitSuf[i]];
      }
      newS[w] = acc;
    }
    S.swap(newS);
  }

  for (int w = 0; w < nWords; ++w) out[w] = S[w];
  return out;
}
