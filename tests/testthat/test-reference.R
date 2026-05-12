test_that("depth = 0 returns 1 regardless of input", {
  expect_equal(PathSignatuR:::.signatureRef(matrix(rnorm(10), 5, 2), 0), 1)
  expect_equal(PathSignatuR:::.signatureRef(matrix(0, 1, 3), 0), 1)
  expect_equal(PathSignatuR:::.signatureRef(rnorm(8), 0), 1)
})

test_that("T = 1 yields (1, 0, 0, ...) at any depth", {
  X <- matrix(c(0.3, -0.1), 1, 2)
  for (N in 0:3) {
    s <- PathSignatuR:::.signatureRef(X, N)
    expect_equal(s[1], 1)
    if (length(s) > 1) expect_equal(s[-1], rep(0, length(s) - 1))
  }
})

test_that("1D linear path: bitwise (b-a)^k / k!", {
  x <- makeLine1D(a = 0.5, b = 2.7, n = 200)
  s <- PathSignatuR:::.signatureRef(x, 5)
  expect_equal(s, (2.7 - 0.5)^(0:5) / factorial(0:5), tolerance = TOL_STRICT)
})

test_that("1D non-linear path also gives (x_T - x_0)^k / k!", {
  # 1D signature depends only on endpoints; tensor-exp recovers this exactly.
  set.seed(7)
  x <- matrix(cumsum(rnorm(300)), ncol = 1)
  s <- PathSignatuR:::.signatureRef(x, 5)
  delta <- x[300] - x[1]
  expect_equal(s, delta^(0:5) / factorial(0:5), tolerance = TOL_STRICT)
})

test_that("multi-d linear path matches closed form", {
  set.seed(13)
  for (d in 2:4) {
    a <- runif(d, -1, 1); b <- runif(d, -1, 1)
    X <- makeLineND(a, b, n = 150)
    for (N in 1:3) {
      s <- PathSignatuR:::.signatureRef(X, N)
      expect_equal(s, linearPathSignature(a, b, N),
                   tolerance = TOL_STRICT,
                   info = sprintf("d=%d N=%d", d, N))
    }
  }
})

test_that("Chen identity holds for the reference itself", {
  # Splitting a path and recombining via tensor product must equal the whole.
  set.seed(23)
  X    <- makeRandomPath(100, 3, seed = 42)
  X1   <- X[1:60, , drop = FALSE]
  X2   <- X[60:100, , drop = FALSE]
  N    <- 2

  sFull <- PathSignatuR:::.signatureRef(X,  N)
  s1    <- PathSignatuR:::.signatureRef(X1, N)
  s2    <- PathSignatuR:::.signatureRef(X2, N)

  words <- enumerateWords(3, N)
  # Manual Chen product reconstruction
  sReconstructed <- vapply(seq_along(words), function(wi) {
    w <- words[[wi]]; k <- length(w)
    if (k == 0L) return(s1[1] * s2[1])
    total <- 0
    for (j in 0:k) {
      pre <- if (j == 0L) integer(0) else w[1:j]
      suf <- if (j == k) integer(0) else w[(j+1):k]
      posPre <- if (length(pre) == 0L) 1L else {
        cum <- (3^length(pre) - 1L) %/% 2L
        cum + 1L + sum((pre - 1L) * 3^((length(pre)-1L):0L))
      }
      posSuf <- if (length(suf) == 0L) 1L else {
        cum <- (3^length(suf) - 1L) %/% 2L
        cum + 1L + sum((suf - 1L) * 3^((length(suf)-1L):0L))
      }
      total <- total + s1[posPre] * s2[posSuf]
    }
    total
  }, numeric(1L))

  expect_equal(sFull, sReconstructed, tolerance = TOL_DEFAULT)
})
