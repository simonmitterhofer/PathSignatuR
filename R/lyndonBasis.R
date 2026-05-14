#' Standard bracketing of a Lyndon word.
#'
#' For `|w| = 1`, returns the letter (atomic). For `|w| >= 2`, factors
#' `w = u . v` where `v` is the longest proper Lyndon suffix of `w`,
#' and returns `list(left = bracket(u), right = bracket(v))`. The
#' recursion gives the canonical free-Lie representation of `w`.
#'
#' @keywords internal
.standardBracket <- function(w) {
  k <- length(w)
  if (k == 1L) return(w[1L])
  for (i in 2:k) {
    suffix <- w[i:k]
    if (.isLyndon(suffix)) {
      return(list(left  = .standardBracket(w[1:(i - 1L)]),
                  right = .standardBracket(suffix)))
    }
  }
  stop(".standardBracket: input is not a Lyndon word")
}

#' Expand a bracket to tensor algebra coordinates.
#'
#' `[a, b] = a (x) b - b (x) a`, linear in both arguments. Singletons
#' map to the corresponding level-1 basis vector. Output is a numeric
#' vector of length `sum(d^(0:N))` in `enumerateWords` order.
#'
#' @keywords internal
.bracketToTensor <- function(bracket, d, N) {
  L <- as.integer(sum(as.double(d)^(0:N)))
  if (is.atomic(bracket) && length(bracket) == 1L) {
    v <- numeric(L)
    v[1L + as.integer(bracket)] <- 1
    return(v)
  }
  a <- .bracketToTensor(bracket$left,  d, N)
  b <- .bracketToTensor(bracket$right, d, N)
  unname(tensorProduct(a, b, N) - tensorProduct(b, a, N))
}

#' Build the Lyndon -> tensor projection matrix.
#'
#' Column `i` is the standard bracket of Lyndon word `i` expanded in
#' tensor basis. Shape: `sum(d^(0:N)) x length(words)`.
#'
#' @keywords internal
.lyndonProjectionMatrix <- function(d, N, words = NULL) {
  if (is.null(words)) words <- lyndonWords(d, N)
  L <- as.integer(sum(as.double(d)^(0:N)))
  M <- matrix(0, L, length(words))
  for (i in seq_along(words)) {
    M[, i] <- .bracketToTensor(.standardBracket(words[[i]]), d, N)
  }
  M
}
