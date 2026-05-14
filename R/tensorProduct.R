#' Truncated tensor product
#'
#' Computes the truncated tensor product of two elements `a` and `b`
#' of the truncated tensor algebra. For each output word `w`,
#' `(a (x) b)[w] = sum over splits w = u.v of a[u] * b[v]`.
#'
#' Both inputs must be numeric vectors of equal length matching
#' `sum(d^(0:depth))` for some integer `d`, which is inferred from
#' the length. The output sits in the same `enumerateWords` ordering
#' and carries canonical word names.
#'
#' @param a,b numeric vectors of equal length. Their length must equal
#'   `sum(d^(0:depth))` for some positive integer `d`. Must be finite.
#' @param depth non-negative integer truncation level. Must match the
#'   `depth` used to produce `a` and `b`.
#'
#' @return a named numeric vector of the same length as `a` and `b`,
#'   in canonical word order.
#'
#' @section Concatenation of signatures (Chen's identity):
#' If `X1` and `X2` are paths sharing a join point, then
#' `signature(X1.X2) = signature(X1) (x) signature(X2)`. This makes
#' `tensorProduct()` the natural way to combine signatures of
#' adjacent path segments.
#'
#' @examples
#' a <- signature(matrix(cumsum(rnorm(20)), 10, 2), depth = 2)
#' b <- signature(matrix(cumsum(rnorm(20)), 10, 2), depth = 2)
#' tensorProduct(a, b, depth = 2)
#'
#' @seealso [signature()], [enumerateWords()]
#' @export
tensorProduct <- function(a, b, depth) {
  depth <- .validateDepth(depth)
  if (!is.numeric(a) || !is.numeric(b)) {
    stop("`a` and `b` must be numeric")
  }
  if (length(a) != length(b)) {
    stop("`a` and `b` must have equal length")
  }
  if (any(!is.finite(a)) || any(!is.finite(b))) {
    stop("`a` and `b` must contain only finite values")
  }

  L <- length(a)
  d <- .inferDim(L, depth)

  if (depth == 0L) {
    out <- as.numeric(unname(a) * unname(b))
  } else {
    ws  <- .getWorkspace(d, depth)
    out <- tensor_product_cpp(unname(as.numeric(a)),
                              unname(as.numeric(b)), ws)
  }
  names(out) <- .wordNames(d, depth)
  out
}
