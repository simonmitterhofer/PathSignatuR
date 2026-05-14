test_that("1D path: log-signature has only level-1 entry equal to the increment", {
  set.seed(1)
  x <- cumsum(rnorm(50))
  delta <- x[length(x)] - x[1]
  ls <- logSignature(x, depth = 5)
  expect_equal(unname(ls[1]), 0,     tolerance = TOL_STRICT)  # empty word
  expect_equal(unname(ls[2]), delta, tolerance = TOL_DEFAULT) # level-1
  expect_equal(unname(ls[-(1:2)]), rep(0, length(ls) - 2L),
               tolerance = TOL_DEFAULT)
})

test_that("multi-D linear path: log-signature is (0, delta_1, ..., delta_d, 0, ...)", {
  set.seed(3)
  for (d in 2:4) {
    a <- runif(d, -1, 1); b <- runif(d, -1, 1)
    X <- makeLineND(a, b, n = 100)
    ls <- logSignature(X, depth = 3)
    L <- length(ls)
    expect_equal(unname(ls[1]), 0, tolerance = TOL_STRICT,
                 info = sprintf("d=%d empty", d))
    expect_equal(unname(ls[2:(d + 1L)]), b - a, tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d level-1", d))
    expect_equal(unname(ls[(d + 2L):L]), rep(0, L - d - 1L),
                 tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d level >= 2", d))
  }
})

test_that("log/exp roundtrip: exp(log(S)) == S for actual signatures", {
  set.seed(7)
  for (d in 2:3) for (N in 1:3) {
    X  <- makeRandomPath(40, d, seed = d * 100 + N)
    s  <- signature(X, depth = N)
    ls <- logSignature(X, depth = N)
    sRecon <- PathSignatuR:::.tensorExpRef(unname(ls), N)
    expect_equal(unname(s), sRecon, tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("depth = 0 returns a single zero", {
  ls <- logSignature(matrix(rnorm(10), 5, 2), depth = 0)
  expect_length(ls, 1L)
  expect_equal(unname(ls), 0)
  expect_equal(names(ls), "")
})

test_that("includeLevelZero = FALSE drops the empty-word entry", {
  X <- makeRandomPath(20, 2, seed = 5)
  full <- logSignature(X, depth = 3, includeLevelZero = TRUE)
  drop <- logSignature(X, depth = 3, includeLevelZero = FALSE)
  expect_length(drop, length(full) - 1L)
  expect_equal(unname(drop), unname(full[-1L]), tolerance = TOL_STRICT)
})

test_that("output is named in canonical word order", {
  X <- makeRandomPath(20, 2, seed = 9)
  ls <- logSignature(X, depth = 2)
  expect_equal(names(ls), c("", "1", "2", "1,1", "1,2", "2,1", "2,2"))
})

test_that("sep controls name formatting", {
  X <- makeRandomPath(20, 2, seed = 11)
  expect_equal(names(logSignature(X, 2, sep = ""))[4:6],
               c("11", "12", "21"))
})

test_that("invalid input is rejected", {
  expect_error(logSignature(NULL, 2),                 "NULL")
  expect_error(logSignature("a", 2),                  "numeric")
  expect_error(logSignature(matrix(NA_real_, 3, 2), 2), "finite")
  expect_error(logSignature(matrix(rnorm(6), 3, 2), -1), "non-negative")
})

test_that("log-signature antisymmetry: level-2 (i,j) and (j,i) sum to zero in linear paths", {
  # Algebraic check: for linear paths, level-2 vanishes entirely.
  # More general: level-2 logSig satisfies ls[ij] + ls[ji] = 0 only
  # for closed paths. Here we just confirm the linear-path zero.
  set.seed(13)
  a <- c(0.3, -0.7); b <- c(1.1, 0.4)
  X <- makeLineND(a, b, n = 50)
  ls <- logSignature(X, depth = 2)
  expect_equal(unname(ls["1,1"]), 0, tolerance = TOL_DEFAULT)
  expect_equal(unname(ls["1,2"] + ls["2,1"]), 0, tolerance = TOL_DEFAULT)
  expect_equal(unname(ls["2,2"]), 0, tolerance = TOL_DEFAULT)
})
