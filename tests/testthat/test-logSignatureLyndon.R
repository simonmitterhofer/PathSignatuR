test_that("output length equals the Witt dimension", {
  for (d in 2:3) for (N in 1:4) {
    X <- makeRandomPath(20, d, seed = d * 10 + N)
    expect_equal(length(logSignatureLyndon(X, depth = N)),
                 length(lyndonWords(d, N)),
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("output is named with Lyndon words joined by sep", {
  X    <- makeRandomPath(20, 2, seed = 1)
  beta <- logSignatureLyndon(X, depth = 3)
  expect_equal(names(beta), c("1", "2", "1,2", "1,1,2", "1,2,2"))
})

test_that("sep controls name formatting", {
  X    <- makeRandomPath(20, 2, seed = 1)
  beta <- logSignatureLyndon(X, depth = 3, sep = "")
  expect_equal(names(beta), c("1", "2", "12", "112", "122"))
})

test_that("round-trip: expandLyndon(logSignatureLyndon(X)) = logSignature(X)", {
  set.seed(7)
  for (d in 2:3) for (N in 1:4) {
    X    <- makeRandomPath(20, d, seed = d * 100 + N)
    beta <- logSignatureLyndon(X, depth = N)
    expect_equal(unname(expandLyndon(beta, d, N)),
                 unname(logSignature(X, depth = N)),
                 tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("linear path: only level-1 Lyndon coords nonzero, equal to increments", {
  set.seed(3)
  for (d in 2:3) {
    a <- runif(d, -1, 1); b <- runif(d, -1, 1)
    X    <- makeLineND(a, b, n = 100)
    beta <- logSignatureLyndon(X, depth = 4)
    expect_equal(unname(beta[1:d]), b - a, tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d level-1", d))
    expect_equal(unname(beta[-(1:d)]), rep(0, length(beta) - d),
                 tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d higher levels", d))
  }
})

test_that("1D path: single Lyndon coord equals the increment", {
  x    <- cumsum(rnorm(50))
  beta <- logSignatureLyndon(x, depth = 4)
  expect_length(beta, 1L)
  expect_equal(unname(beta), x[length(x)] - x[1], tolerance = TOL_DEFAULT)
})

test_that("depth = 0 returns numeric(0)", {
  X <- matrix(rnorm(10), 5, 2)
  expect_identical(logSignatureLyndon(X, depth = 0), numeric(0))
})

test_that("vector input is treated as a 1D path", {
  v    <- cumsum(rnorm(20))
  beta <- logSignatureLyndon(v, depth = 3)
  expect_length(beta, 1L)
  expect_equal(unname(beta), v[20] - v[1], tolerance = TOL_DEFAULT)
})

test_that("invalid input is rejected", {
  expect_error(logSignatureLyndon(NULL, 2),                 "NULL")
  expect_error(logSignatureLyndon("not numeric", 2),        "numeric")
  expect_error(logSignatureLyndon(matrix(NA_real_, 3, 2), 2), "finite")
  expect_error(logSignatureLyndon(matrix(rnorm(6), 3, 2), -1), "non-negative")
  expect_error(logSignatureLyndon(matrix(rnorm(6), 3, 2), 2,
                                  sep = c(",", ";")), "single string")
})

test_that("Lyndon cache is populated and reused", {
  rm(list = ls(envir = PathSignatuR:::.lyndonCache, all.names = TRUE),
     envir = PathSignatuR:::.lyndonCache)
  X <- makeRandomPath(20, 2, seed = 1)
  logSignatureLyndon(X, depth = 3)
  expect_setequal(ls(envir = PathSignatuR:::.lyndonCache), "2.3")
  logSignatureLyndon(X, depth = 3)
  expect_setequal(ls(envir = PathSignatuR:::.lyndonCache), "2.3")
})
