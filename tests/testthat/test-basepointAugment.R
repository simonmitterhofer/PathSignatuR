test_that("default basepoint is the origin", {
  X  <- makeRandomPath(10, 2, seed = 1)
  Xb <- basepointAugment(X)
  expect_equal(dim(Xb), c(11L, 2L))
  expect_equal(Xb[1, ], c(0, 0))
  expect_equal(unname(Xb[-1, ]), unname(X))
})

test_that("custom basepoint is used as the first row", {
  X  <- makeRandomPath(10, 3, seed = 2)
  bp <- c(1, -1, 0.5)
  Xb <- basepointAugment(X, basepoint = bp)
  expect_equal(Xb[1, ], bp)
  expect_equal(unname(Xb[-1, ]), unname(X))
})

test_that("vector input is treated as a 1D path", {
  v  <- rnorm(8)
  Xb <- basepointAugment(v)
  expect_equal(dim(Xb), c(9L, 1L))
  expect_equal(Xb[1, 1], 0)
  expect_equal(unname(Xb[-1, 1]), v)
})

test_that("output has stripped dimnames", {
  X <- matrix(rnorm(6), 3, 2, dimnames = list(c("a", "b", "c"),
                                              c("x", "y")))
  expect_null(dimnames(basepointAugment(X)))
})

test_that("T = 1 input gives a 2-row path", {
  X  <- matrix(c(0.3, -0.1), 1, 2)
  Xb <- basepointAugment(X)
  expect_equal(dim(Xb), c(2L, 2L))
  expect_equal(Xb[1, ], c(0, 0))
  expect_equal(Xb[2, ], c(0.3, -0.1))
})

test_that("invalid input is rejected", {
  X <- matrix(rnorm(6), 3, 2)
  expect_error(basepointAugment(NULL),                    "NULL")
  expect_error(basepointAugment("not numeric"),           "numeric")
  expect_error(basepointAugment(matrix(NA_real_, 3, 2)),  "finite")
  expect_error(basepointAugment(X, basepoint = "a"),      "numeric")
  expect_error(basepointAugment(X, basepoint = c(1, 2, 3)),
               "length 2")
  expect_error(basepointAugment(X, basepoint = c(1, NA)), "finite")
  expect_error(basepointAugment(X, basepoint = matrix(c(1, 2), 1, 2)),
               "vector")
})

test_that("breaks translation invariance: level-1 terms shift", {
  X     <- makeRandomPath(50, 2, seed = 7)
  shift <- c(3, -2)

  # Without basepoint: signature is translation invariant.
  sPlain        <- signature(X, depth = 1)
  sPlainShifted <- signature(X + matrix(shift, nrow(X), 2, byrow = TRUE),
                             depth = 1)
  expect_equal(unname(sPlain), unname(sPlainShifted),
               tolerance = TOL_STRICT)

  # With basepoint at origin: shifting changes level-1.
  sBp        <- signature(basepointAugment(X), depth = 1)
  sBpShifted <- signature(basepointAugment(
    X + matrix(shift, nrow(X), 2, byrow = TRUE)),
    depth = 1)
  expect_false(isTRUE(all.equal(unname(sBp), unname(sBpShifted),
                                tolerance = TOL_DEFAULT)))
})

test_that("level-1 terms equal the terminal absolute position", {
  # With basepoint at 0, S^i = X[T, i] - 0 = X[T, i].
  X <- makeRandomPath(40, 3, seed = 11)
  s <- signature(basepointAugment(X), depth = 1,
                 includeLevelZero = FALSE)
  expect_equal(unname(s), X[nrow(X), ], tolerance = TOL_STRICT)
})

test_that("custom basepoint b: level-1 terms equal X[T, ] - b", {
  X  <- makeRandomPath(40, 2, seed = 13)
  bp <- c(0.7, -0.2)
  s  <- signature(basepointAugment(X, basepoint = bp), depth = 1,
                  includeLevelZero = FALSE)
  expect_equal(unname(s), X[nrow(X), ] - bp, tolerance = TOL_STRICT)
})

test_that("composition with timeAugment puts basepoint at t = 0", {
  X  <- makeRandomPath(10, 2, seed = 17)
  Xt <- timeAugment(basepointAugment(X))
  expect_equal(dim(Xt), c(11L, 3L))
  expect_equal(unname(Xt[1, 1]), 0)         # basepoint sits at t = 0
  expect_equal(unname(Xt[nrow(Xt), 1]), 1)  # last row at t = 1
})
