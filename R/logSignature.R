#' Log-signature of a discrete path
#'
#' Computes the formal logarithm of the truncated signature in the
#' truncated tensor algebra:
#' `logSig(X) = sum_{k = 1}^{N} (-1)^{k+1} (S - e)^{(x) k} / k`,
#' where `S = signature(X, depth)` and `e = (1, 0, 0, ...)`. The series
#' is finite because `(S - e)` has zero level-0 coefficient, so
#' `(S - e)^{(x) k}` has minimum level `k` and contributes nothing
#' above level `N`.
#'
#' @param X a numeric matrix (T x d) or vector (treated as 1D path).
#' @param depth non-negative integer truncation level.
#' @param includeLevelZero if `TRUE` (default), the returned vector
#'   starts with the empty-word entry (always `0`, since `log 1 = 0`).
#'   If `FALSE`, that entry is dropped.
#' @param sep separator for word names (default `","`).
#'
#' @return a named numeric vector of length `sum(d^(0:depth))` (or
#'   `... - 1` when `includeLevelZero = FALSE`).
#'
#' @section Representation:
#' The log-signature lies in the free Lie algebra of rank `d`, which
#' is a strict subspace of the truncated tensor algebra. The returned
#' vector represents the log-signature as a general tensor algebra
#' element in canonical word order; entries are therefore not
#' independent (they satisfy linear relations induced by the Jacobi
#' identity and antisymmetry of the Lie bracket). Projection onto a
#' Hall or Lyndon basis is left for a future release.
#'
#' @section Linear paths:
#' For a linear path from `a` to `b`, the signature is the tensor
#' exponential of the increment `delta = b - a`, so the log-signature
#' is exactly `(0, delta_1, ..., delta_d, 0, 0, ...)` -- level-1
#' coordinates carry the increment, all higher levels are zero.
#'
#' @examples
#' x <- cumsum(rnorm(50))
#' logSignature(x, depth = 4)   # 1D: only level-1 entry is nonzero
#'
#' X <- matrix(cumsum(rnorm(60)), 30, 2)
#' logSignature(X, depth = 3)
#'
#' @seealso [signature()], [tensorProduct()], [tensorInverse()]
#' @export
logSignature <- function(X, depth, includeLevelZero = TRUE, sep = ",") {
  X     <- .validatePath(X)
  depth <- .validateDepth(depth)
  .validateSignatureOptions(sep, includeLevelZero)

  d   <- ncol(X)
  ws  <- .getWorkspace(d, depth)
  out <- log_sig_terminal_cpp(X, ws)
  names(out) <- .wordNames(d, depth, sep = sep)
  if (!includeLevelZero) out <- out[-1L]
  out
}
