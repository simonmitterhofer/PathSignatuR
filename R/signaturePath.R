#' Path-valued truncated signature
#'
#' Returns the signature of every prefix of the path. The result is a
#' `T x p` matrix whose row `t` is the signature of `X[1:t, , drop = FALSE]`;
#' row 1 is the identity `(1, 0, 0, ...)` and row `T` matches
#' `signature(X, depth)`.
#'
#' Useful for rolling / window applications, plotting signature
#' trajectories, and computing partial signatures on any sub-interval as
#' `result[end, ]` "minus" `result[start, ]` (via tensor-algebra inverse,
#' deferred to v0.4's `chenProduct()` / `tensorInverse()`).
#'
#' @param X a numeric matrix (T x d) or vector (treated as 1D).
#' @param depth non-negative integer truncation level.
#' @param includeLevelZero if `TRUE` (default), keep the first column
#'   (empty-word, always 1). If `FALSE`, drop it.
#' @param sep separator for word names (default `","`).
#'
#' @return a numeric matrix of shape `T x p` (or `T x (p - 1)`). Column
#'   names are word strings; row names are `NULL`.
#'
#' @examples
#' X <- matrix(cumsum(rnorm(60)), 30, 2)
#' sp <- signaturePath(X, depth = 2)
#' sp[1, ]              # identity (1, 0, 0, 0, 0, 0, 0)
#' all.equal(sp[30, ], signature(X, depth = 2))
#'
#' @seealso [signature()], [signatureBatch()]
#' @export
signaturePath <- function(X, depth, includeLevelZero = TRUE, sep = ",") {
  X     <- .validatePath(X)
  depth <- .validateDepth(depth)
  .validateSignatureOptions(sep, includeLevelZero)

  d   <- ncol(X)
  out <- sig_path_cpp(X, depth)
  colnames(out) <- .wordNames(d, depth, sep = sep)
  if (!includeLevelZero) out <- out[, -1L, drop = FALSE]
  out
}
