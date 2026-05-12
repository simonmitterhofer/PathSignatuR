# Shared test fixtures and tolerance constants.
# helper-* files are auto-loaded by testthat before any test runs.

TOL_STRICT  <- 1e-12   # closed-form agreement under tensor-exp
TOL_DEFAULT <- 1e-10   # C++ vs R reference, Chen identity, etc.

# 1D straight line: x(t) = a + (b - a) t.
makeLine1D <- function(a = 0, b = 1, n = 100) {
  matrix(seq(a, b, length.out = n), ncol = 1)
}

# d-dim straight line from a to b, with n equally spaced points.
makeLineND <- function(a, b, n = 100) {
  stopifnot(length(a) == length(b))
  d <- length(a)
  t <- seq(0, 1, length.out = n)
  outer(t, rep(1, d)) * matrix(b - a, n, d, byrow = TRUE) +
    matrix(a, n, d, byrow = TRUE)
}

# Closed-form signature of a linear path from a to b, in enumerateWords order.
# For word (i_1, ..., i_k): prod_j (b[i_j] - a[i_j]) / k!.
linearPathSignature <- function(a, b, depth) {
  d <- length(a)
  words <- enumerateWords(d, depth)
  delta <- b - a
  vapply(words, function(w) {
    k <- length(w)
    if (k == 0L) 1 else prod(delta[w]) / factorial(k)
  }, numeric(1L))
}

# Reproducible random walk (cumsum of N(0, 1) increments).
makeRandomPath <- function(T, d, seed = 1L) {
  set.seed(seed)
  matrix(cumsum(rnorm(T * d)), T, d)
}
