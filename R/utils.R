#' Validate and coerce a path argument.
#' @keywords internal
.validatePath <- function(X) {
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
  dimnames(X) <- NULL
  X
}

#' Validate the depth (truncation level).
#' @keywords internal
.validateDepth <- function(depth) {
  if (length(depth) != 1L) stop("`depth` must be a single value")
  if (!is.numeric(depth) || !is.finite(depth)) {
    stop("`depth` must be a finite non-negative integer")
  }
  if (depth < 0 || depth != trunc(depth)) {
    stop("`depth` must be a non-negative integer")
  }
  as.integer(depth)
}

#' Build the names vector for a (d, depth) signature.
#' Empty word -> "", otherwise letters joined by `sep`.
#' @keywords internal
.wordNames <- function(d, depth, sep = ",") {
  words <- enumerateWords(d, depth)
  vapply(words, function(w) {
    if (length(w) == 0L) "" else paste(w, collapse = sep)
  }, character(1L))
}
