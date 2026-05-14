#' Path-valued log-signature
#'
#' Returns the log-signature at every prefix of the path. The result is
#' a `T x p` matrix whose row `t` is `logSignature(X[1:t, , drop = FALSE], depth)`;
#' row 1 is zero (log of the identity), and row `T` matches
#' `logSignature(X, depth)`.
#'
#' @param X a numeric matrix (T x d) or vector (treated as 1D path).
#' @param depth non-negative integer truncation level.
#' @param includeLevelZero if `TRUE` (default), keep the first column
#'   (empty word, always 0). If `FALSE`, drop it.
#' @param sep separator for word names (default `","`).
#'
#' @return a numeric matrix of shape `T x p` (or `T x (p - 1)`). Column
#'   names are word strings; row names are `NULL`.
#'
#' @section Implementation:
#' One-pass C++ kernel: maintains the running signature `S` and, after
#' each Chen update, expands `log(S)` via the truncated Neumann series
#' `sum_{k=1}^N (-1)^{k+1} (S - e)^{(x) k} / k`. Uses the same cached
#' `(d, depth)` workspace as the signature kernel.
#'
#' @examples
#' X  <- matrix(cumsum(rnorm(60)), 30, 2)
#' lp <- logSignaturePath(X, depth = 2)
#' lp[1, ]                                    # zero row (log of identity)
#' all.equal(lp[30, ], logSignature(X, depth = 2))
#'
#' @seealso [logSignature()], [signaturePath()]
#' @export
logSignaturePath <- function(X, depth, includeLevelZero = TRUE, sep = ",") {
  X     <- .validatePath(X)
  depth <- .validateDepth(depth)
  .validateSignatureOptions(sep, includeLevelZero)

  d   <- ncol(X)
  ws  <- .getWorkspace(d, depth)
  out <- log_sig_path_cpp(X, ws)
  colnames(out) <- .wordNames(d, depth, sep = sep)
  if (!includeLevelZero) out <- out[, -1L, drop = FALSE]
  out
}
