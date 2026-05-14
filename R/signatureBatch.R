#' Truncated path signatures of a batch of paths
#'
#' Computes terminal signatures for many paths in one call, sharing the
#' `(d, depth)` precompute across all paths. Equivalent to applying
#' [signature()] to each path individually but considerably faster for
#' large batches of short paths.
#'
#' @param X one of:
#'   * a list of `T_i x d` numeric matrices (variable `T_i` allowed; all
#'     must share the same `d`),
#'   * a 3D numeric array of shape `T x d x nPaths` (uniform `T`),
#'   * a single matrix or vector (treated as a 1-path batch).
#' @param depth non-negative integer truncation level.
#' @param includeLevelZero if `TRUE` (default), the first column is the
#'   empty-word entry. If `FALSE`, that column is dropped.
#' @param sep separator for word names (default `","`).
#'
#' @return a numeric matrix of shape `nPaths x p`, where `p = sum(d^(0:depth))`
#'   (or `p - 1` when `includeLevelZero = FALSE`). Row `i` is the terminal
#'   signature of the `i`-th path. Column names are word strings. Row
#'   names come from `names(X)` for a list input or `dimnames(X)[[3]]`
#'   for a 3D array.
#'
#' @examples
#' paths <- list(
#'   matrix(cumsum(rnorm(40)), 20, 2),
#'   matrix(cumsum(rnorm(60)), 30, 2),
#'   matrix(cumsum(rnorm(50)), 25, 2)
#' )
#' signatureBatch(paths, depth = 2)
#'
#' @seealso [signature()], [signaturePath()]
#' @export
signatureBatch <- function(X, depth, includeLevelZero = TRUE, sep = ",") {
  paths     <- .toPathList(X)
  depth     <- .validateDepth(depth)
  .validateSignatureOptions(sep, includeLevelZero)

  pathNames <- names(paths)
  if (length(paths) < 1L) stop("`X` must contain at least one path")

  d <- ncol(paths[[1]])
  for (i in seq_along(paths)) {
    if (ncol(paths[[i]]) != d) {
      stop(sprintf("path %d has ncol = %d, expected %d",
                   i, ncol(paths[[i]]), d))
    }
  }

  ws  <- .getWorkspace(d, depth)
  out <- sig_batch_cpp(paths, ws)
  colnames(out) <- .wordNames(d, depth, sep = sep)
  if (!is.null(pathNames)) rownames(out) <- pathNames
  if (!includeLevelZero) out <- out[, -1L, drop = FALSE]
  out
}

#' Coerce a batch input to a list of validated path matrices.
#' Dispatches on type: list, 3D array, or single path.
#' @keywords internal
.toPathList <- function(X) {
  if (is.null(X)) stop("`X` must not be NULL")

  # List of paths (data.frames excluded -- handled as single path below).
  if (is.list(X) && !is.data.frame(X)) {
    paths <- lapply(seq_along(X), function(i) {
      tryCatch(.validatePath(X[[i]]),
               error = function(e) {
                 stop(sprintf("path %d: %s", i, conditionMessage(e)),
                      call. = FALSE)
               })
    })
    names(paths) <- names(X)
    return(paths)
  }

  # 3D array T x d x nPaths.
  if (is.array(X) && length(dim(X)) == 3L) {
    if (!is.numeric(X)) stop("`X` must be numeric")
    if (any(!is.finite(X))) stop("`X` must contain only finite values")
    nP <- dim(X)[3L]
    paths <- lapply(seq_len(nP), function(i) {
      m <- X[, , i, drop = FALSE]
      dim(m) <- dim(m)[1:2]
      m
    })
    names(paths) <- dimnames(X)[[3L]]
    return(paths)
  }

  # Single matrix or vector.
  list(.validatePath(X))
}
