test_that("output shape is T x p with column names", {
  X <- matrix(rnorm(20), 10, 2)
  sp <- signaturePath(X, depth = 2)
  expect_equal(dim(sp), c(10L, 7L))
  expect_equal(colnames(sp), c("", "1", "2", "1,1", "1,2", "2,1", "2,2"))
  expect_null(rownames(sp))
})

test_that("row 1 is the identity (1, 0, ..., 0)", {
  X <- matrix(rnorm(40), 20, 2)
  sp <- signaturePath(X, depth = 3)
  expect_equal(unname(sp[1, 1]), 1)
  expect_equal(unname(sp[1, -1L]), rep(0, ncol(sp) - 1L))
})

test_that("last row matches signature() on the full path", {
  set.seed(11)
  for (d in 1:3) {
    X <- makeRandomPath(30, d, seed = d)
    for (N in 0:3) {
      sp  <- signaturePath(X, depth = N)
      sig <- signature(X, depth = N)
      expect_equal(unname(sp[nrow(sp), ]), unname(sig),
                   tolerance = TOL_STRICT,
                   info = sprintf("d=%d N=%d", d, N))
    }
  }
})

test_that("row t matches signature() on the prefix X[1:t, ]", {
  set.seed(7)
  X <- makeRandomPath(20, 2, seed = 5)
  sp <- signaturePath(X, depth = 2)
  for (t in 2:nrow(X)) {
    expect_equal(unname(sp[t, ]),
                 unname(signature(X[1:t, , drop = FALSE], depth = 2)),
                 tolerance = TOL_DEFAULT,
                 info = sprintf("t=%d", t))
  }
})

test_that("1D linear path: row t equals (x_t - x_1)^k / k!", {
  x  <- seq(0, 2, length.out = 50)
  sp <- signaturePath(matrix(x, ncol = 1), depth = 4)
  for (t in 1:50) {
    expect_equal(unname(sp[t, ]),
                 (x[t] - x[1])^(0:4) / factorial(0:4),
                 tolerance = TOL_STRICT,
                 info = sprintf("t=%d", t))
  }
})

test_that("vector input gives T x p", {
  v <- rnorm(15)
  expect_equal(dim(signaturePath(v, depth = 2)), c(15L, 3L))
})

test_that("includeLevelZero = FALSE drops the first column", {
  X <- matrix(rnorm(20), 10, 2)
  full <- signaturePath(X, depth = 2, includeLevelZero = TRUE)
  drop <- signaturePath(X, depth = 2, includeLevelZero = FALSE)
  expect_equal(ncol(drop), ncol(full) - 1L)
  expect_equal(unname(drop), unname(full[, -1L, drop = FALSE]))
})

test_that("sep controls column-name format", {
  X <- matrix(rnorm(20), 10, 2)
  expect_equal(colnames(signaturePath(X, 2, sep = ""))[4:6],
               c("11", "12", "21"))
})

test_that("depth = 0 returns a T x 1 column of ones", {
  X  <- matrix(rnorm(20), 10, 2)
  sp <- signaturePath(X, depth = 0)
  expect_equal(dim(sp), c(10L, 1L))
  expect_equal(unname(sp), matrix(1, 10, 1))
})

test_that("T = 1 returns a single identity row", {
  X  <- matrix(c(0.3, -0.1), 1, 2)
  sp <- signaturePath(X, depth = 3)
  expect_equal(dim(sp), c(1L, 15L))
  expect_equal(unname(sp[1, 1]), 1)
  expect_equal(unname(sp[1, -1L]), rep(0, 14))
})

test_that("invalid input errors cleanly", {
  expect_error(signaturePath(NULL, 2),                 "NULL")
  expect_error(signaturePath("not numeric", 2),        "numeric")
  expect_error(signaturePath(matrix(NA_real_, 3, 2), 2), "finite")
  expect_error(signaturePath(matrix(rnorm(6), 3, 2), -1), "non-negative")
})
