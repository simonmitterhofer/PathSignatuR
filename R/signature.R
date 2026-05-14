#' Truncated path signature of a discrete path
#'
#' Computes the terminal value of the truncated path signature of `X` up to
#' level `depth`. Uses the exact piecewise-linear formulation (segment-wise
#' tensor exponentials and Chen products), so closed forms such as
#' `(b - a)^k / k!` for linear paths are recovered to machine precision.
#'
#' @param X a numeric matrix of shape `T x d` (rows = ordered observations,
#'   columns = channels). A numeric vector is treated as a 1D path
#'   (`T x 1`). Must be finite.
#' @param depth non-negative integer truncation level. `depth = 0` returns
#'   only the level-0 term (`1`).
#' @param includeLevelZero if `TRUE` (default), the returned vector starts
#'   with the empty-word entry (value `1`, name `""`). If `FALSE`, that
#'   entry is dropped.
#' @param sep string separating letters in word names (default `","`).
#'   For `d <= 9`, `sep = ""` gives compact `"11"`, `"12"`, ... names;
#'   for `d >= 10` the default is required to stay unambiguous.
#'
#' @return a named numeric vector. Length is `sum(d^(0:depth))` if
#'   `includeLevelZero = TRUE`, else `sum(d^(0:depth)) - 1`.
#'
#' @section Time augmentation:
#' `signature()` does not augment with time. To include time as a channel,
#' pre-process with [timeAugment()].
#'
#' @section Note on name clash:
#' `PathSignatuR::signature()` masks `methods::signature()` on attach.
#'
#' @examples
#' # 1D path: signature term at level k equals (x_T - x_0)^k / k!
#' x <- cumsum(rnorm(50))
#' s <- signature(x, depth = 4)
#' expected <- (x[50] - x[1])^(0:4) / factorial(0:4)
#' all.equal(unname(s), expected)
#'
#' @seealso [enumerateWords()], [signatureBatch()], [signaturePath()]
#' @export
signature <- function(X, depth, includeLevelZero = TRUE, sep = ",") {
  X     <- .validatePath(X)
  depth <- .validateDepth(depth)
  .validateSignatureOptions(sep, includeLevelZero)

  d   <- ncol(X)
  out <- sig_terminal_cpp(X, depth)
  names(out) <- .wordNames(d, depth, sep = sep)
  if (!includeLevelZero) out <- out[-1L]
  out
}
