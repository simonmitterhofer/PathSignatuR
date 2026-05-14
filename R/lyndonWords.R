#' Lyndon words over an alphabet of size d
#'
#' Enumerates all Lyndon words of length 1 to `depth` over the alphabet
#' `{1, ..., dim}`. A Lyndon word is a non-empty word that is strictly
#' lexicographically less than all its proper rotations. Equivalently,
#' a Lyndon word is the unique lex-smallest representative of its
#' rotation class, restricted to aperiodic words.
#'
#' Lyndon words of length up to `N` index a basis of the free Lie
#' algebra of rank `d` truncated at level `N`. The count at length `k`
#' is given by Witt's formula
#' `L(d, k) = (1/k) sum_{j | k} mu(j) d^(k/j)`,
#' which is strictly smaller than `d^k` for `d >= 2, k >= 2`.
#'
#' @param dim   alphabet size, positive integer.
#' @param depth maximum word length, non-negative integer.
#'
#' @return A list of integer vectors, in canonical order: by length,
#'   then lex within length. Empty list for `depth = 0`.
#'
#' @examples
#' lyndonWords(2, 3)
#' #> [[1]] 1
#' #> [[2]] 2
#' #> [[3]] 1 2
#' #> [[4]] 1 1 2
#' #> [[5]] 1 2 2
#'
#' # Counts match Witt's formula
#' lengths(lyndonWords(3, 4))
#' length(lyndonWords(3, 4))    # 32
#'
#' @seealso [enumerateWords()], [logSignatureLyndon()]
#' @export
lyndonWords <- function(dim, depth) {
  if (length(dim) != 1L || !is.numeric(dim) || dim < 1 || dim != trunc(dim)) {
    stop("`dim` must be a positive integer")
  }
  if (length(depth) != 1L || !is.numeric(depth) || depth < 0 ||
      depth != trunc(depth) || !is.finite(depth)) {
    stop("`depth` must be a non-negative integer")
  }
  d <- as.integer(dim); N <- as.integer(depth)
  if (N == 0L) return(list())

  words <- enumerateWords(d, N)[-1L]            # drop empty word
  words[vapply(words, .isLyndon, logical(1L))]
}

#' Check whether a non-empty integer word is Lyndon.
#'
#' A word is Lyndon iff strictly lex-less than all its proper rotations.
#' Brute O(n^2) check; n is small (truncation depth) so fine.
#'
#' @keywords internal
.isLyndon <- function(w) {
  n <- length(w)
  if (n == 0L) return(FALSE)
  if (n == 1L) return(TRUE)
  for (i in 1:(n - 1L)) {
    rot <- c(w[(i + 1L):n], w[1:i])
    if (!.lexLess(w, rot)) return(FALSE)
  }
  TRUE
}

#' Strict lex-less comparison of equal-length integer vectors.
#' @keywords internal
.lexLess <- function(a, b) {
  ne <- a != b
  if (!any(ne)) return(FALSE)
  i <- which.max(ne)
  a[i] < b[i]
}
