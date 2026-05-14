#' Truncated shuffle product
#'
#' Computes the shuffle product of two elements `a` and `b` of the
#' truncated tensor algebra. For each output word `w` of length `k`,
#' `(a ⧢ b)[w] = sum over subsets S of {1..k} of a[w|S] * b[w|S^c]`.
#' Equivalently, sums over all `2^k` interleavings of subwords.
#'
#' The shuffle product is **commutative and associative**, unlike the
#' tensor product. The identity tensor `(1, 0, 0, ...)` is the unit
#' for both products.
#'
#' @param a,b numeric vectors of equal length matching `sum(d^(0:depth))`
#'   for some positive integer `d`. Must be finite.
#' @param depth non-negative integer truncation level.
#'
#' @return a named numeric vector of the same length as `a` and `b`.
#'
#' @section Shuffle identity for signatures:
#' For any path `X` with truncated signature `S`, and any two words
#' `u`, `v` with `|u| + |v| <= depth`,
#' `<S, u> * <S, v> = <S, u ⧢ v>` (scalar multiplication on the left
#' equals shuffle on the right). This means coordinates of a signature
#' aren't independent: products of low-level coordinates are linear
#' combinations of high-level coordinates.
#'
#' @examples
#' a <- c(1, 0.5, -0.2, 0.1, 0, 0.3, -0.4)
#' b <- c(1, -0.1, 0.2, 0, 0.5, 0, 0.1)
#' shuffleProduct(a, b, depth = 2)
#'
#' @seealso [tensorProduct()], [tensorInverse()], [signature()]
#' @export
shuffleProduct <- function(a, b, depth) {
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
    ws  <- .getShuffleWorkspace(d, depth)
    out <- shuffle_product_cpp(unname(as.numeric(a)),
                               unname(as.numeric(b)), ws)
  }
  names(out) <- .wordNames(d, depth)
  out
}
