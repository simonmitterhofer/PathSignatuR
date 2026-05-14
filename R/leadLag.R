#' Lead-lag transform of a path
#'
#' Applies the Flint-Hambly-Lyons lead-lag transform: each channel of
#' `X` is duplicated into a "lead" and a "lag" copy that advance in
#' alternation, producing a path of length `2T - 1` and dimension `2d`.
#'
#' The point of the transform is that the level-2 cross-signature term
#' `S^{lead_i, lag_i}` of the transformed path recovers the quadratic
#' variation of channel `i` of the original path (up to a factor of
#' 1/2). This lets the signature pick up volatility-like information
#' that a piecewise-linear path would otherwise miss.
#'
#' For a 1D path `(x_1, x_2, x_3)` the output is
#' ```
#'         lead_1   lag_1
#'    t=1:  x_1     x_1
#'    t=2:  x_2     x_1     (lead jumped)
#'    t=3:  x_2     x_2     (lag caught up)
#'    t=4:  x_3     x_2     (lead jumped)
#'    t=5:  x_3     x_3     (lag caught up)
#' ```
#'
#' @param X numeric matrix (T x d) or vector (treated as 1D path).
#'
#' @return a numeric matrix of shape `(2T - 1) x (2d)`. Columns are
#'   interleaved `lead_1, lag_1, lead_2, lag_2, ...` and named
#'   accordingly. Row names are stripped.
#'
#' @section Quadratic variation:
#' For `X` a 1D path and `LL = leadLag(X)`, the level-2 cross term
#' `S^{lead, lag}(LL) - S^{lag, lead}(LL)` equals the signed area
#' enclosed by the lead-lag staircase, which in turn equals (up to
#' sign and a factor) the quadratic variation of `X`. See Flint,
#' Hambly & Lyons (2016) for the precise statement.
#'
#' @section Sizing:
#' Lead-lag doubles the channel count, which inflates the truncated
#' signature size from `sum(d^(0:N))` to `sum((2d)^(0:N))`. For
#' `(d, N) = (3, 4)` that's 121 vs 1365 terms. Plan depth accordingly.
#'
#' @examples
#' # 1D path: classic staircase
#' leadLag(c(1, 2, 3))
#' #>      lead_1 lag_1
#' #> [1,]    1     1
#' #> [2,]    2     1
#' #> [3,]    2     2
#' #> [4,]    3     2
#' #> [5,]    3     3
#'
#' # 2D path: columns are lead_1, lag_1, lead_2, lag_2
#' X <- matrix(rnorm(20), 10, 2)
#' LL <- leadLag(X)
#' dim(LL)
#' #> [1] 19  4
#' colnames(LL)
#' #> "lead_1" "lag_1" "lead_2" "lag_2"
#'
#' @seealso [signature()], [basepointAugment()], [timeAugment()]
#' @export
leadLag <- function(X) {
  X <- .validatePath(X)
  tLen <- nrow(X); d <- ncol(X)

  nOut <- 2L * tLen - 1L
  out  <- matrix(0, nrow = nOut, ncol = 2L * d)

  # Lead column for channel i lives at output column 2i - 1; lag at 2i.
  # Lead at output row r:
  #   r odd  (r = 2t - 1, t = 1..T):  X[t, ]
  #   r even (r = 2t,     t = 1..T-1): X[t + 1, ]   (just jumped)
  # Lag at output row r:
  #   r = 1:                     X[1, ]
  #   r odd  (r = 2t - 1, t > 1): X[t, ]            (just caught up)
  #   r even (r = 2t,     t > 0): X[t, ]            (still holding)
  #
  # Compactly: lead row r holds X[ceiling((r + 1) / 2), ],
  #           lag  row r holds X[ceiling(r / 2),       ].
  leadIdx <- ceiling((seq_len(nOut) + 1L) / 2L)
  lagIdx  <- ceiling( seq_len(nOut)       / 2L)

  for (i in seq_len(d)) {
    out[, 2L * i - 1L] <- X[leadIdx, i]
    out[, 2L * i]      <- X[lagIdx,  i]
  }

  colnames(out) <- as.vector(rbind(paste0("lead_", seq_len(d)),
                                   paste0("lag_",  seq_len(d))))
  rownames(out) <- NULL
  out
}
