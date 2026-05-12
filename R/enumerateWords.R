#' Enumerate all multi-indices over an alphabet of size d
#'
#' Generates the list of words used to index the truncated signature, in the
#' canonical level-major, lexicographic-within-level ordering: empty word
#' first (position 1), then all level-1 words (positions 2 to d+1), then
#' all level-2 words, and so on. Within each level the last letter varies
#' fastest.
#'
#' For a level-k word at position `w >= 2`, its parent (level k-1) is at
#' position `floor((w-2)/d) + 1` and its last letter is `((w-2) %% d) + 1`.
#' This identity is what the signature kernel uses to build level-k columns
#' from level-(k-1) columns.
#'
#' @param dim   alphabet size (number of path channels), a positive integer
#' @param depth truncation level, a non-negative integer
#'
#' @return A list of integer vectors of length `sum(dim^(0:depth))`. The
#'   empty word is `integer(0)`.
#'
#' @examples
#' enumerateWords(2, 2)
#' #> [[1]] integer(0)
#' #> [[2]] 1
#' #> [[3]] 2
#' #> [[4]] 1 1
#' #> [[5]] 1 2
#' #> [[6]] 2 1
#' #> [[7]] 2 2
#'
#' @export
enumerateWords <- function(dim, depth) {
  if (length(dim) != 1L || !is.numeric(dim) || dim < 1 || dim != trunc(dim)) {
    stop("`dim` must be a positive integer")
  }
  if (length(depth) != 1L || !is.numeric(depth) || depth < 0 ||
      depth != trunc(depth) || !is.finite(depth)) {
    stop("`depth` must be a non-negative integer")
  }
  dim   <- as.integer(dim)
  depth <- as.integer(depth)

  nWords <- as.integer(sum(as.double(dim) ^ (0:depth)))
  words  <- vector("list", nWords)
  words[[1]] <- integer(0)
  if (depth == 0L) return(words)

  idx <- 1L; prevStart <- 1L; prevEnd <- 1L
  for (k in seq_len(depth)) {
    for (j in prevStart:prevEnd) {
      for (letter in seq_len(dim)) {
        idx <- idx + 1L
        words[[idx]] <- c(words[[j]], letter)
      }
    }
    prevStart <- prevEnd + 1L
    prevEnd   <- idx
  }
  words
}
