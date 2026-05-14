#' Truncated tensor inverse
#'
#' Computes `a^{-1}` in the truncated tensor algebra: the unique element
#' satisfying `a (x) a^{-1} = e`, where `e = (1, 0, 0, ...)` is the
#' identity tensor. Requires `a[1] = 1` (the standard condition under
#' which the inverse exists as a Neumann series).
#'
#' The implementation uses the truncated Neumann expansion
#' `a^{-1} = sum_{k = 0}^{N} (e - a)^{(x) k}`,
#' which is finite at truncation level `N` because `(e - a)` has zero
#' level-0 coefficient and thus `(e - a)^{(x) k}` has minimum level `k`.
#'
#' @param a numeric vector. Length must match `sum(d^(0:depth))` for
#'   some positive integer `d`. Level-0 coefficient `a[1]` must equal 1
#'   (within `1e-12`). Must be finite.
#' @param depth non-negative integer truncation level.
#'
#' @return a named numeric vector of the same length as `a`.
#'
#' @section Signature reversal:
#' If `X` is a path and `X_rev` is the same path traversed in reverse,
#' then `signature(X) (x) signature(X_rev) = e`, so
#' `signature(X_rev) = tensorInverse(signature(X))`. This is the
#' standard way to invert signatures under path reversal.
#'
#' @examples
#' X <- matrix(cumsum(rnorm(40)), 20, 2)
#' s <- signature(X, depth = 3)
#' sInv <- tensorInverse(s, depth = 3)
#' # Verify: s (x) sInv ≈ e
#' tensorProduct(s, sInv, depth = 3)
#'
#' @seealso [tensorProduct()], [signature()]
#' @export
tensorInverse <- function(a, depth) {
  depth <- .validateDepth(depth)
  if (!is.numeric(a)) stop("`a` must be numeric")
  if (any(!is.finite(a))) stop("`a` must contain only finite values")
  L <- length(a)
  d <- .inferDim(L, depth)

  if (abs(a[1] - 1) > 1e-12) {
    stop("`a[1]` (level-0 coefficient) must equal 1")
  }

  e <- c(1, rep(0, L - 1))
  if (depth == 0L) {
    names(e) <- ""
    return(e)
  }

  a       <- unname(as.numeric(a))
  c_vec   <- e - a                   # zero level-0
  result  <- e
  current <- e                       # (e - a)^{(x) 0}
  for (k in seq_len(depth)) {
    current <- unname(tensorProduct(current, c_vec, depth))
    result  <- result + current
  }
  names(result) <- .wordNames(d, depth)
  result
}
