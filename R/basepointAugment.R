#' Prepend a basepoint row to a path
#'
#' Adds a fixed reference point as the first row of `X`. This breaks the
#' translation invariance of the signature: with a basepoint at the
#' origin, level-1 signature terms become absolute terminal positions
#' `X[T, ]` rather than increments `X[T, ] - X[1, ]`, and higher-level
#' terms pick up cross-information between the starting position and
#' the rest of the path.
#'
#' @param X numeric matrix (T x d) or vector (treated as 1D path).
#' @param basepoint numeric vector of length `d`, or `NULL` (default)
#'   for `rep(0, d)`. Must be finite. No recycling: a length-`d`
#'   vector is required when supplied.
#'
#' @return a numeric matrix of shape `(T + 1) x d`. Dimnames are
#'   stripped, matching [timeAugment()].
#'
#' @section When to use this:
#' Basepoint augmentation is appropriate when the absolute level of the
#' path carries signal that translation invariance would discard:
#' bounded state variables (interest rates, volatility indices, credit
#' spreads), categorical embeddings where the integer encoding matters,
#' or spatial data in a fixed coordinate frame.
#'
#' For asset-price work it is usually unnecessary: practitioners
#' typically operate on log-prices or log-returns, where the level
#' carries no information the signature should care about, and
#' translation invariance is a feature rather than a bug.
#'
#' @section Composition with timeAugment:
#' If you want both, apply `basepointAugment()` first and `timeAugment()`
#' second. That way the basepoint row sits at `t = 0` and the original
#' first row at `t = 1 / T` (under default unit scaling), which is the
#' usual intended geometry. Reversing the order puts the basepoint at
#' a non-zero time, which is rarely what's wanted.
#'
#' @examples
#' X <- matrix(rnorm(20), 10, 2)
#' Xb <- basepointAugment(X)
#' dim(Xb)
#' #> [1] 11  2
#' Xb[1, ]
#' #> [1] 0 0
#'
#' # Custom basepoint
#' basepointAugment(X, basepoint = c(1, -1))[1, ]
#' #> [1]  1 -1
#'
#' # Recommended composition with timeAugment
#' signature(timeAugment(basepointAugment(X)), depth = 2)
#'
#' @seealso [timeAugment()], [signature()]
#' @export
basepointAugment <- function(X, basepoint = NULL) {
  X <- .validatePath(X)
  d <- ncol(X)

  if (is.null(basepoint)) {
    basepoint <- rep(0, d)
  } else {
    if (!is.numeric(basepoint) || !is.vector(basepoint)) {
      stop("`basepoint` must be a numeric vector")
    }
    if (length(basepoint) != d) {
      stop(sprintf("`basepoint` must have length %d (ncol(X)), got %d",
                   d, length(basepoint)))
    }
    if (any(!is.finite(basepoint))) {
      stop("`basepoint` must contain only finite values")
    }
    basepoint <- as.numeric(basepoint)
  }

  out <- rbind(basepoint, X, deparse.level = 0)
  dimnames(out) <- NULL
  out
}
