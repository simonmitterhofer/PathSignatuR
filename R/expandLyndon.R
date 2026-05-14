#' Expand Lyndon coordinates to the tensor basis
#'
#' Inverse of [logSignatureLyndon()]: given Lyndon coordinates `beta`,
#' returns the corresponding free-Lie-algebra element in the tensor
#' basis as a named numeric vector of length `sum(dim^(0:depth))`.
#'
#' Computationally one matrix-vector product against the cached
#' bracket-to-tensor matrix.
#'
#' @param beta  numeric vector of length `length(lyndonWords(dim, depth))`.
#' @param dim   alphabet size, positive integer.
#' @param depth non-negative integer truncation level.
#' @param sep   separator for tensor-word names (default `","`).
#'
#' @return Named numeric vector of length `sum(dim^(0:depth))` in
#'   `enumerateWords` order. For `depth = 0`, returns `c("" = 0)`.
#'
#' @examples
#' X     <- matrix(cumsum(rnorm(60)), 30, 2)
#' beta  <- logSignatureLyndon(X, depth = 3)
#' lsT   <- expandLyndon(beta, dim = 2, depth = 3)
#' all.equal(unname(lsT), unname(logSignature(X, depth = 3)))
#'
#' @seealso [logSignatureLyndon()], [logSignature()]
#' @export
expandLyndon <- function(beta, dim, depth, sep = ",") {
  if (!is.numeric(beta) || !is.vector(beta)) {
    stop("`beta` must be a numeric vector")
  }
  if (any(!is.finite(beta))) {
    stop("`beta` must contain only finite values")
  }
  if (length(dim) != 1L || !is.numeric(dim) || dim < 1 || dim != trunc(dim)) {
    stop("`dim` must be a positive integer")
  }
  depth <- .validateDepth(depth)
  if (!is.character(sep) || length(sep) != 1L || is.na(sep)) {
    stop("`sep` must be a single string")
  }
  d <- as.integer(dim)

  if (depth == 0L) {
    if (length(beta) != 0L) {
      stop(sprintf("`beta` must have length 0 (Witt dim for d=%d, N=0), got %d",
                   d, length(beta)))
    }
    out <- 0
    names(out) <- ""
    return(out)
  }

  cache    <- .getLyndonProjection(d, depth)
  expected <- length(cache$words)
  if (length(beta) != expected) {
    stop(sprintf("`beta` must have length %d (Witt dim for d=%d, N=%d), got %d",
                 expected, d, depth, length(beta)))
  }

  out <- as.numeric(cache$M %*% unname(beta))
  names(out) <- .wordNames(d, depth, sep = sep)
  out
}
