#' Chen product (alias for tensorProduct)
#'
#' The Chen product of two signatures is mathematically identical to the
#' truncated tensor product: for paths `X1`, `X2` sharing a join point,
#' `signature(X1.X2) = signature(X1) (x) signature(X2)`. This function
#' is provided as a named alias of [tensorProduct()] for code that
#' wants to emphasise the signature-concatenation semantics.
#'
#' Use `chenProduct()` when combining signatures of adjacent path
#' segments; use `tensorProduct()` when working with general tensor
#' algebra elements. Both call the same underlying kernel.
#'
#' @inheritParams tensorProduct
#'
#' @return a named numeric vector of the same length as `a` and `b`.
#'
#' @examples
#' X <- matrix(cumsum(rnorm(60)), 30, 2)
#' X1 <- X[1:15, , drop = FALSE]
#' X2 <- X[15:30, , drop = FALSE]
#' s1 <- signature(X1, depth = 2)
#' s2 <- signature(X2, depth = 2)
#' sFull <- chenProduct(s1, s2, depth = 2)
#' all.equal(unname(sFull), unname(signature(X, depth = 2)))
#'
#' @seealso [tensorProduct()], [signature()]
#' @export
chenProduct <- function(a, b, depth) tensorProduct(a, b, depth)
