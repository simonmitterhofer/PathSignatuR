#' Resample a path to a new time grid
#'
#' Interpolates `X` onto a new grid, preserving piecewise-linearity so
#' the result composes cleanly with [signature()] and friends. The
#' source path is treated as living on `seq(0, 1, length.out = T)`;
#' supply either `n` for a uniform regrid or `time` for an explicit
#' target grid.
#'
#' @param X numeric matrix (T x d, with T >= 2) or vector (treated as
#'   a 1D path). Single-row paths are rejected: there is nothing to
#'   interpolate between.
#' @param n positive integer. Number of equispaced points to produce
#'   on `[0, 1]`. Exactly one of `n` or `time` must be supplied.
#' @param time numeric vector of target times. Must be finite,
#'   non-decreasing, and contained in `[0, 1]` (up to a tolerance of
#'   `1e-10`).
#' @param method one of `"linear"` (default) or `"step"`. `"step"` is
#'   piecewise-constant left-continuous: the value at time `t` is
#'   `X[s, ]` for the largest source time `s <= t`.
#'
#' @return a numeric matrix of shape `n x d` (or `length(time) x d`).
#'   Dimnames are stripped.
#'
#' @section No splines:
#' Only linear and step interpolation are supported by design. Spline
#' or higher-order interpolation would produce a path that is not
#' piecewise linear, breaking the assumption the signature kernel
#' rests on. If you need smoother resampling for a different purpose,
#' do it outside this package before calling [signature()].
#'
#' @section Use cases:
#' Two common motivations: aligning multiple paths to a common grid
#' before [signatureBatch()], and densifying or thinning a single
#' path for rolling-window applications. Step interpolation is the
#' natural choice for state variables observed at irregular times
#' (last-known value carries forward).
#'
#' @examples
#' # Uniform regrid to 50 points
#' X <- matrix(cumsum(rnorm(60)), 30, 2)
#' pathInterpolate(X, n = 50)
#'
#' # Custom target grid
#' pathInterpolate(X, time = c(0, 0.25, 0.5, 0.75, 1))
#'
#' # Step interpolation
#' pathInterpolate(X, n = 100, method = "step")
#'
#' @seealso [signature()], [signatureBatch()]
#' @export
pathInterpolate <- function(X, n = NULL, time = NULL,
                            method = c("linear", "step")) {
  X      <- .validatePath(X)
  method <- match.arg(method)
  tLen   <- nrow(X); d <- ncol(X)

  if (tLen < 2L) {
    stop("`X` must have at least 2 rows; nothing to interpolate")
  }
  if (is.null(n) && is.null(time)) {
    stop("supply exactly one of `n` or `time`")
  }
  if (!is.null(n) && !is.null(time)) {
    stop("supply exactly one of `n` or `time`, not both")
  }

  # --- Build target grid ---
  if (!is.null(n)) {
    if (length(n) != 1L || !is.numeric(n) || !is.finite(n) ||
        n != trunc(n) || n < 1) {
      stop("`n` must be a positive integer")
    }
    n <- as.integer(n)
    target <- seq(0, 1, length.out = n)
  } else {
    if (!is.numeric(time) || !is.vector(time)) {
      stop("`time` must be a numeric vector")
    }
    if (length(time) < 1L) {
      stop("`time` must have length at least 1")
    }
    if (any(!is.finite(time))) {
      stop("`time` must contain only finite values")
    }
    tol <- 1e-10
    if (min(time) < -tol || max(time) > 1 + tol) {
      stop("`time` must lie in [0, 1]")
    }
    if (any(diff(time) < 0)) {
      stop("`time` must be non-decreasing")
    }
    target <- pmin(pmax(as.numeric(time), 0), 1)  # clamp tiny FP excess
  }

  # --- Interpolate ---
  src    <- seq(0, 1, length.out = tLen)
  approxMethod <- if (method == "linear") "linear" else "constant"
  out <- matrix(0, nrow = length(target), ncol = d)
  for (i in seq_len(d)) {
    out[, i] <- stats::approx(x = src, y = X[, i], xout = target,
                              method = approxMethod, f = 0,
                              rule = 2)$y
  }
  out
}
