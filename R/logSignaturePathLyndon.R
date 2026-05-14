#' Path-valued log-signature in the Lyndon basis
#'
#' Returns the Lyndon-basis log-signature at every prefix of the path.
#' Row `t` is `logSignatureLyndon(X[1:t, , drop = FALSE], depth)`;
#' row 1 is zero (log of identity), and row `T` matches
#' `logSignatureLyndon(X, depth)`.
#'
#' @param X numeric matrix (T x d) or vector (treated as 1D path).
#' @param depth non-negative integer truncation level.
#' @param sep separator for Lyndon-word names (default `","`).
#'
#' @return Numeric matrix of shape `T x q` where
#'   `q = length(lyndonWords(d, depth))`. Column names are Lyndon
#'   words joined by `sep`; row names are `NULL`. For `depth = 0`
#'   returns a `T x 0` matrix.
#'
#' @section Implementation:
#' Computes the tensor-basis running log-signature via
#' `logSignaturePath()` (one C++ call), then batches the projection
#' as a single `qr.coef` against the cached QR. The projection step
#' is one matrix back-solve regardless of `T`.
#'
#' @examples
#' X  <- matrix(cumsum(rnorm(60)), 30, 2)
#' bp <- logSignaturePathLyndon(X, depth = 3)
#' dim(bp)                          # 30 x 5
#' bp[1, ]                          # zero (log of identity)
#' all.equal(unname(bp[30, ]),
#'           unname(logSignatureLyndon(X, depth = 3)))
#'
#' @seealso [logSignatureLyndon()], [logSignaturePath()]
#' @export
logSignaturePathLyndon <- function(X, depth, sep = ",") {
  X     <- .validatePath(X)
  depth <- .validateDepth(depth)
  if (!is.character(sep) || length(sep) != 1L || is.na(sep)) {
    stop("`sep` must be a single string")
  }

  d <- ncol(X); T <- nrow(X)
  if (depth == 0L) return(matrix(0, T, 0L))

  LP    <- logSignaturePath(X, depth = depth)
  cache <- .getLyndonProjection(d, depth)
  B     <- qr.coef(cache$qr, t(unname(LP)))           # q x T
  out   <- t(B)
  colnames(out) <- .lyndonWordNames(cache$words, sep = sep)
  rownames(out) <- NULL
  out
}
