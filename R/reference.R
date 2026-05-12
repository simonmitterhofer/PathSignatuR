#' Pure-R reference signature implementation
#'
#' Exact piecewise-linear signature via segment-wise tensor exponentials
#' and Chen products. For a path X with increments delta_1, ..., delta_(T-1):
#'
#'   S(X) = exp(delta_1) (x) exp(delta_2) (x) ... (x) exp(delta_(T-1))
#'
#' where exp() is the truncated tensor exponential and (x) is the truncated
#' tensor product. Bitwise exact for piecewise-linear paths; closed forms
#' like (b - a)^k / k! are recovered to machine precision.
#'
#' Used in tests as the comparison baseline for the C++ kernel.
#'
#' @param X     numeric vector or matrix (T x d)
#' @param depth non-negative integer
#'
#' @return numeric vector of length `sum(d^(0:depth))`, unnamed.
#'
#' @keywords internal
.signatureRef <- function(X, depth) {
  # --- validate ---
  X <- .validatePath(X)
  N <- .validateDepth(depth)
  tLen <- nrow(X); d <- ncol(X)
  nWords <- as.integer(sum(as.double(d) ^ (0:N)))

  # Identity tensor (1, 0, 0, ...).
  S <- numeric(nWords); S[1L] <- 1
  if (N == 0L || tLen < 2L) return(S)

  # --- precompute split table for the Chen product ----------------------
  # splitTable[[w]] is a (|w|+1) x 2 integer matrix; row (j+1) holds the
  # flat positions of (prefix of length j, suffix of length |w|-j).
  words <- enumerateWords(d, N)

  wordPos <- function(letters) {
    k <- length(letters)
    if (k == 0L) return(1L)
    if (d == 1L) return(k + 1L)
    cumLevel <- (d^k - 1L) %/% (d - 1L)      # 1 + d + ... + d^{k-1}
    as.integer(cumLevel + 1L + sum((letters - 1L) * d^((k - 1L):0L)))
  }

  splitTable <- lapply(words, function(word) {
    k <- length(word)
    if (k == 0L) return(matrix(c(1L, 1L), nrow = 1L))
    out <- matrix(0L, nrow = k + 1L, ncol = 2L)
    for (j in 0:k) {
      pre <- if (j == 0L) integer(0) else word[1:j]
      suf <- if (j == k) integer(0) else word[(j + 1L):k]
      out[j + 1L, ] <- c(wordPos(pre), wordPos(suf))
    }
    out
  })

  # --- main loop: S := S (x) exp(delta_i) for each segment --------------
  inc <- diff(X)
  for (i in seq_len(tLen - 1L)) {
    delta <- inc[i, ]

    # exp(delta) via parent recurrence: E[w of length k] = E[parent] * delta[last] / k
    # (so the level-k term is delta^{(x)k} / k!).
    E <- numeric(nWords); E[1L] <- 1
    for (w in 2:nWords) {
      parentCol <- ((w - 2L) %/% d) + 1L
      letter    <- ((w - 2L) %%  d) + 1L
      k         <- length(words[[w]])
      E[w] <- E[parentCol] * delta[letter] / k
    }

    # Chen product: newS[w] = sum over splits S[prefix] * E[suffix].
    newS <- numeric(nWords)
    for (w in seq_len(nWords)) {
      st <- splitTable[[w]]
      newS[w] <- sum(S[st[, 1L]] * E[st[, 2L]])
    }
    S <- newS
  }
  S
}
