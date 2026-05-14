#' Prepend a time channel to a path
#'
#' Adds a time coordinate as the first column of the path. By default the
#' time grid is normalised to `[0, 1]` (Remark 3.4 in Cuchiero & Möller),
#' which keeps signature values on a stable scale across rolling windows
#' of differing length. Supply a custom grid via `time`, or use
#' `scale = "none"` to keep elapsed numeric time.
#'
#' @param X      numeric matrix (T x d) or vector (treated as 1D path)
#' @param time   optional numeric vector of length `T` giving the time
#'   coordinate of each row. If `NULL`, defaults to `seq(0, 1, length.out = T)`
#'   under `scale = "unit"`, or `0, 1, ..., T-1` under `scale = "none"`.
#' @param scale  one of `"unit"` (default, normalise to `[0, 1]`) or
#'   `"none"` (keep `time` as-is).
#'
#' @return a numeric matrix with `T` rows and `d + 1` columns. The first
#'   column is named `"time"`; user channel names on `X` are stripped.
#'
#' @examples
#' X <- matrix(rnorm(20), 10, 2)
#' Xt <- timeAugment(X)
#' colnames(Xt)
#' #> "time" NULL NULL  (only "time" is named)
#'
#' @export
timeAugment <- function(X, time = NULL, scale = c("unit", "none")) {
  X     <- .validatePath(X)
  scale <- match.arg(scale)
  tLen  <- nrow(X)

  if (is.null(time)) {
    tCol <- if (scale == "unit") {
      if (tLen == 1L) 0 else seq(0, 1, length.out = tLen)
    } else {
      seq_len(tLen) - 1L
    }
  } else {
    if (length(time) != tLen) {
      stop("`time` must have length equal to nrow(X)")
    }
    if (any(!is.finite(time))) stop("`time` must be finite")
    time <- as.numeric(time)
    tCol <- if (scale == "unit") {
      rng <- range(time)
      if (diff(rng) == 0) rep(0, tLen) else (time - rng[1L]) / diff(rng)
    } else {
      time
    }
  }

  out <- cbind(tCol, X, deparse.level = 0)
  colnames(out) <- c("time", rep("", ncol(X)))
  out
}
