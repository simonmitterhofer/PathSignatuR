test_that("default unit-scaled grid is [0, 1] with T points", {
  X <- matrix(rnorm(20), 10, 2)
  Xt <- timeAugment(X)
  expect_equal(dim(Xt), c(10L, 3L))
  expect_equal(colnames(Xt)[1], "time")
  expect_equal(Xt[, 1], seq(0, 1, length.out = 10))
  expect_equal(unname(Xt[, -1]), unname(X))
})

test_that("scale = 'none' with NULL time gives 0..T-1", {
  X <- matrix(rnorm(15), 5, 3)
  expect_equal(timeAugment(X, time = NULL, scale = "none")[, 1], 0:4)
})

test_that("explicit time vector is accepted and normalised under 'unit'", {
  X  <- matrix(rnorm(10), 5, 2)
  Xt <- timeAugment(X, time = c(10, 12, 14, 16, 18), scale = "unit")
  expect_equal(Xt[, 1], seq(0, 1, length.out = 5))
})

test_that("explicit time vector is preserved under 'none'", {
  X  <- matrix(rnorm(10), 5, 2)
  Xt <- timeAugment(X, time = c(10, 12, 14, 16, 18), scale = "none")
  expect_equal(Xt[, 1], c(10, 12, 14, 16, 18))
})

test_that("vector input is treated as a 1D path and yields T x 2", {
  Xt <- timeAugment(rnorm(8))
  expect_equal(dim(Xt), c(8L, 2L))
  expect_equal(colnames(Xt), c("time", ""))
})

test_that("invalid input is rejected", {
  expect_error(timeAugment("not numeric"), "numeric")
  expect_error(timeAugment(matrix(NA_real_, 3, 2)), "finite")
  expect_error(timeAugment(matrix(rnorm(6), 3, 2), time = c(1, 2)),
               "length equal to nrow")
  expect_error(timeAugment(matrix(rnorm(6), 3, 2), time = c(1, NA, 3)),
               "finite")
  expect_error(timeAugment(matrix(rnorm(6), 3, 2), scale = "bogus"))
})

test_that("constant time vector under 'unit' degenerates safely", {
  X <- matrix(rnorm(6), 3, 2)
  expect_equal(timeAugment(X, time = c(5, 5, 5), scale = "unit")[, 1],
               rep(0, 3))
})

test_that("output strips user channel names", {
  X  <- matrix(rnorm(6), 3, 2, dimnames = list(NULL, c("ret1", "ret2")))
  Xt <- timeAugment(X)
  expect_equal(colnames(Xt), c("time", "", ""))
})

test_that("integration: signature on augmented path has (d+1)-channel size", {
  X <- makeRandomPath(50, 2, seed = 5)
  expect_length(signature(timeAugment(X), depth = 2), 1 + 3 + 9)
})

test_that("T = 1 input is handled (unit scale gives time = 0)", {
  X  <- matrix(c(0.3, -0.1), 1, 2)
  Xt <- timeAugment(X)
  expect_equal(dim(Xt), c(1L, 3L))
  expect_equal(unname(Xt[1, 1]), 0)
})
