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

#' Validate the shared `sep` / `includeLevelZero` options.
#' @keywords internal
.validateSignatureOptions <- function(sep, includeLevelZero) {
  if (!is.character(sep) || length(sep) != 1L || is.na(sep)) {
    stop("`sep` must be a single string")
  }
  if (!is.logical(includeLevelZero) || length(includeLevelZero) != 1L ||
      is.na(includeLevelZero)) {
    stop("`includeLevelZero` must be TRUE or FALSE")
  }
  invisible(NULL)
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

#' Infer alphabet size d from a (d, depth) signature length.
#' @keywords internal
.inferDim <- function(L, depth) {
  if (depth == 0L) {
    if (L == 1L) return(1L)
    stop(sprintf("length %d incompatible with depth 0 (expected 1)", L))
  }
  upper <- ceiling(L^(1 / depth)) + 1L
  for (d in 1L:upper) {
    if (sum(as.double(d)^(0:depth)) == L) return(d)
  }
  stop(sprintf("length %d is not sum(d^(0:%d)) for any integer d >= 1",
               L, depth))
}

# Package-level cache of SigWorkspace XPtrs, keyed by "d.depth".
# Per-session; refilled on demand by .getWorkspace().
.workspaceCache <- new.env(parent = emptyenv())

#' Get a cached SigWorkspace XPtr for the given (d, depth).
#' @keywords internal
.getWorkspace <- function(d, depth) {
  key <- paste(d, depth, sep = ".")
  ws  <- .workspaceCache[[key]]
  if (is.null(ws)) {
    ws <- build_sig_workspace(as.integer(d), as.integer(depth))
    .workspaceCache[[key]] <- ws
  }
  ws
}
# Package-level cache of ShuffleWorkspace XPtrs.
.shuffleCache <- new.env(parent = emptyenv())

#' Get a cached ShuffleWorkspace XPtr for the given (d, depth).
#' @keywords internal
.getShuffleWorkspace <- function(d, depth) {
  key <- paste(d, depth, sep = ".")
  ws  <- .shuffleCache[[key]]
  if (is.null(ws)) {
    ws <- build_shuffle_workspace(as.integer(d), as.integer(depth))
    .shuffleCache[[key]] <- ws
  }
  ws
}
