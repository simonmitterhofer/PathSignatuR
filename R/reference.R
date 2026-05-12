#' Pure-R reference signature implementation
#'
#' Mirrors the parent-index recurrence literally: builds the full
#' path-valued signature `T x p` via the Stratonovich midpoint rule, then
#' returns its terminal row. Used in tests to verify the C++ kernel.
#'
#' Intentionally written for clarity (matches the parent-index formula
#' one-for-one), not speed. Not exported.
#'
#' @param X     numeric vector or matrix (T x d)
#' @param depth non-negative integer
#'
#' @return numeric vector of length `sum(d^(0:depth))`, unnamed.
#'
#' @keywords internal
.signatureRef <- function(X, depth) {
  # --- validate --------------------------------------------------------
  if (is.null(X)) stop("`X` must not be NULL")
  if (is.data.frame(X)) X <- as.matrix(X)
  if (is.vector(X) && !is.list(X)) X <- matrix(X, ncol = 1L)
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric vector or matrix (T x d)")
  }
  if (nrow(X) < 1L || ncol(X) < 1L) {
    stop("`X` must have at least one row and one column")
  }
  if (any(!is.finite(X))) stop("`X` must contain only finite values")
  if (length(depth) != 1L || !is.numeric(depth) || depth < 0 ||
      depth != trunc(depth) || !is.finite(depth)) {
    stop("`depth` must be a non-negative integer")
  }
  dimnames(X) <- NULL
  depth <- as.integer(depth)

  # --- compute ---------------------------------------------------------
  tLen <- nrow(X); d <- ncol(X)
  nWords <- as.integer(sum(as.double(d) ^ (0:depth)))
  sig <- matrix(0, tLen, nWords)
  sig[, 1L] <- 1
  if (depth == 0L || tLen < 2L) return(sig[tLen, ])

  inc <- diff(X)
  for (w in 2:nWords) {
    parentCol <- ((w - 2L) %/% d) + 1L
    letter    <- ((w - 2L) %%  d) + 1L
    mid <- (sig[-tLen, parentCol] + sig[-1L, parentCol]) / 2
    sig[, w] <- c(0, cumsum(mid * inc[, letter]))
  }
  sig[tLen, ]
}
