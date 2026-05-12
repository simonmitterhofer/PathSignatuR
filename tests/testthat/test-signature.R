test_that("output shape and naming, includeLevelZero = TRUE", {
  X <- makeRandomPath(20, 3)
  s <- signature(X, depth = 2)
  expect_length(s, 1 + 3 + 9)
  expect_equal(s[[1]], 1)
  expect_equal(names(s)[1], "")
  expect_equal(names(s)[2:4], c("1", "2", "3"))
  expect_equal(names(s)[5],  "1,1")
  expect_equal(names(s)[13], "3,3")
})

test_that("includeLevelZero = FALSE drops the empty word entirely", {
  X  <- makeRandomPath(20, 2)
  s0 <- signature(X, depth = 3, includeLevelZero = FALSE)
  s1 <- signature(X, depth = 3, includeLevelZero = TRUE)
  expect_length(s0, length(s1) - 1L)
  expect_equal(names(s0)[1], "1")
  expect_false("" %in% names(s0))
  expect_equal(unname(s0), unname(s1[-1L]))
})

test_that("sep argument controls name formatting", {
  X <- makeRandomPath(15, 2)
  expect_equal(names(signature(X, depth = 2, sep = "" ))[4:6],
               c("11", "12", "21"))
  expect_equal(names(signature(X, depth = 2, sep = "-"))[4:6],
               c("1-1", "1-2", "2-1"))
})

test_that("depth = 0 returns c(empty = 1) with sensible name", {
  X <- makeRandomPath(10, 2)
  s <- signature(X, depth = 0)
  expect_length(s, 1L)
  expect_equal(unname(s), 1)
  expect_equal(names(s), "")
  expect_length(signature(X, depth = 0, includeLevelZero = FALSE), 0L)
})

test_that("vector input is coerced to a 1D path", {
  v  <- rnorm(50)
  sV <- signature(v, depth = 3)
  sM <- signature(matrix(v, ncol = 1), depth = 3)
  expect_equal(unname(sV), unname(sM), tolerance = TOL_STRICT)
})

test_that("T = 2, d = 1 reduces to (1, dx, dx^2/2, ...)", {
  X <- matrix(c(0, 0.7), 2, 1)
  s <- signature(X, depth = 4)
  expect_equal(unname(s), 0.7^(0:4) / factorial(0:4),
               tolerance = TOL_STRICT)
})

test_that("invalid input raises informative errors", {
  expect_error(signature("not numeric", 2),  "numeric")
  expect_error(signature(matrix(NA_real_, 3, 2), 2), "finite")
  expect_error(signature(NULL, 2),           "NULL")
  expect_error(signature(matrix(rnorm(6), 3, 2), -1),  "non-negative")
  expect_error(signature(matrix(rnorm(6), 3, 2), 1.5), "non-negative")
  expect_error(signature(matrix(rnorm(6), 3, 2), c(1, 2)), "single value")
  expect_error(signature(matrix(rnorm(6), 3, 2), 2,
                         includeLevelZero = NA), "TRUE or FALSE")
  expect_error(signature(matrix(rnorm(6), 3, 2), 2, sep = c(",", ";")),
               "single string")
})

test_that("1D path: signature term at level k is (x_T - x_0)^k / k!", {
  set.seed(11); x <- cumsum(rnorm(200))
  s <- signature(x, depth = 5)
  expect_equal(unname(s), (x[200] - x[1])^(0:5) / factorial(0:5),
               tolerance = TOL_DEFAULT)
})

test_that("multi-d linear path matches closed form at levels 1..3", {
  set.seed(31)
  for (d in 2:4) {
    a <- runif(d, -1, 1); b <- runif(d, -1, 1)
    X <- makeLineND(a, b, n = 200)
    for (N in 1:3) {
      expect_equal(unname(signature(X, depth = N)),
                   linearPathSignature(a, b, N),
                   tolerance = TOL_STRICT,
                   info = sprintf("d=%d N=%d", d, N))
    }
  }
})

test_that("Chen identity at levels 1 and 2 (multi-d)", {
  for (d in 2:3) {
    X    <- makeRandomPath(200, d, seed = d)
    mid  <- 100L
    X1   <- X[1:mid, , drop = FALSE]
    X2   <- X[mid:nrow(X), , drop = FALSE]   # share the join point

    sFull <- signature(X,  depth = 2)
    s1    <- signature(X1, depth = 2)
    s2    <- signature(X2, depth = 2)

    # Level 1
    for (i in seq_len(d)) {
      key <- as.character(i)
      expect_equal(sFull[[key]], s1[[key]] + s2[[key]],
                   tolerance = TOL_DEFAULT,
                   info = sprintf("level-1 d=%d i=%d", d, i))
    }
    # Level 2: deconcatenation S^{ij} = S^{ij}(1) + S^i(1) S^j(2) + S^{ij}(2)
    for (i in seq_len(d)) for (j in seq_len(d)) {
      key  <- paste(i, j, sep = ",")
      keyI <- as.character(i); keyJ <- as.character(j)
      lhs <- sFull[[key]]
      rhs <- s1[[key]] + s1[[keyI]] * s2[[keyJ]] + s2[[key]]
      expect_equal(lhs, rhs, tolerance = TOL_DEFAULT,
                   info = sprintf("level-2 d=%d (%d,%d)", d, i, j))
    }
  }
})

test_that("C++ kernel matches R reference across (d, N) grid", {
  for (d in 1:4) for (N in 0:3) {
    X    <- makeRandomPath(80, d, seed = 100L + 10L * d + N)
    sCpp <- signature(X, depth = N)
    sRef <- PathSignatuR:::.signatureRef(X, N)
    expect_equal(unname(sCpp), unname(sRef),
                 tolerance = TOL_STRICT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("kernel handles depth = 0 and T = 1 edge cases", {
  expect_equal(unname(signature(matrix(rnorm(10), 5, 2), 0)), 1)
  s <- signature(matrix(c(0.3, -0.1), 1, 2), 2)
  expect_equal(s[[1]], 1)
  expect_equal(unname(s[-1]), rep(0, length(s) - 1))
})

test_that("translation invariance: signature ignores path basepoint", {
  X <- makeRandomPath(100, 3, seed = 99)
  shift <- matrix(rnorm(3), nrow = nrow(X), ncol = 3, byrow = TRUE)
  expect_equal(unname(signature(X,         depth = 3)),
               unname(signature(X + shift, depth = 3)),
               tolerance = TOL_STRICT)
})
