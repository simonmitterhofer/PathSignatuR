#' Log-signature in the Lyndon basis
#'
#' Projects the log-signature onto the Lyndon basis of the free Lie
#' algebra, returning canonical independent coordinates. The Lyndon
#' basis has dimension given by Witt's formula
#' `dim L_N(d) = sum_{k=1}^{N} (1/k) sum_{j | k} mu(j) d^(k/j)`,
#' strictly smaller than the tensor-algebra dimension `sum(d^(0:N))`
#' for `d >= 2, N >= 2`. For `(d, N) = (3, 4)`: 32 vs 121.
#'
#' @param X numeric matrix (T x d) or vector (treated as 1D path).
#' @param depth non-negative integer truncation level.
#' @param sep separator for Lyndon-word names (default `","`).
#'
#' @return Named numeric vector of length `length(lyndonWords(d, depth))`.
#'   Names are Lyndon words joined by `sep`. For `depth = 0`, returns
#'   `numeric(0)`.
#'
#' @section Why this exists:
#' The log-signature lies in the free Lie algebra, a strict subspace
#' of the truncated tensor algebra. The tensor-basis representation
#' carries linear redundancy from the Jacobi identity and bracket
#' antisymmetry. Lyndon coordinates remove this algebraic
#' collinearity, yielding a canonical minimal representation —
#' orthogonal to and complementary with statistical reductions like
#' PCA or JL.
#'
#' @section Implementation:
#' The `(d, depth)` projection matrix `M` (tensor-basis expansions of
#' standard Lyndon brackets) and its QR decomposition are built once
#' per `(d, depth)` and cached. Each call computes the tensor-basis
#' log-signature via the C++ kernel, then back-solves against the
#' cached QR.
#'
#' @examples
#' X    <- matrix(cumsum(rnorm(60)), 30, 2)
#' beta <- logSignatureLyndon(X, depth = 3)
#' length(beta)                   # 5 = length(lyndonWords(2, 3))
#' names(beta)                    # "1" "2" "1,2" "1,1,2" "1,2,2"
#'
#' # Round-trip
#' all.equal(unname(expandLyndon(beta, 2, 3)),
#'           unname(logSignature(X, 3)))
#'
#' @seealso [lyndonWords()], [logSignature()], [expandLyndon()],
#'   [logSignaturePathLyndon()]
#' @export
logSignatureLyndon <- function(X, depth, sep = ",") {
  X     <- .validatePath(X)
  depth <- .validateDepth(depth)
  if (!is.character(sep) || length(sep) != 1L || is.na(sep)) {
    stop("`sep` must be a single string")
  }

  d <- ncol(X)
  if (depth == 0L) return(numeric(0L))

  ls    <- logSignature(X, depth = depth)
  cache <- .getLyndonProjection(d, depth)
  beta  <- qr.coef(cache$qr, unname(ls))
  names(beta) <- .lyndonWordNames(cache$words, sep = sep)
  beta
}
