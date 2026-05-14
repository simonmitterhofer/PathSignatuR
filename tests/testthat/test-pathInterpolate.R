test_that("uniform regrid produces n rows and preserves endpoints", {
  X <- makeRandomPath(30, 2, seed = 1)
  Y <- pathInterpolate(X, n = 50)
  expect_equal(dim(Y), c(50L, 2L))
  expect_equal(Y[1, ],  X[1, ],        tolerance = TOL_STRICT)
  expect_equal(Y[50, ], X[nrow(X), ],  tolerance = TOL_STRICT)
})

test_that("regrid to original n returns the input unchanged", {
  X <- makeRandomPath(30, 3, seed = 2)
  Y <- pathInterpolate(X, n = 30)
  expect_equal(unname(Y), unname(X), tolerance = TOL_STRICT)
})

test_that("linear path is interpolated exactly", {
  x <- seq(0, 2, length.out = 10)
  y <- pathInterpolate(matrix(x, ncol = 1), n = 100)
  expect_equal(y[, 1], seq(0, 2, length.out = 100), tolerance = TOL_STRICT)
})

test_that("multi-d linear path is interpolated exactly", {
  X <- makeLineND(a = c(0, -1, 2), b = c(1, 3, -2), n = 25)
  Y <- pathInterpolate(X, n = 100)
  expect_equal(Y, makeLineND(a = c(0, -1, 2), b = c(1, 3, -2), n = 100),
               tolerance = TOL_STRICT)
})

test_that("step interpolation is piecewise-constant left-continuous", {
  xs <- pathInterpolate(matrix(c(1, 2, 3, 4), ncol = 1), n = 7,
                        method = "step")
  expect_equal(xs[, 1], c(1, 1, 2, 2, 3, 3, 4))
})

test_that("step interpolation preserves the endpoint", {
  X  <- makeRandomPath(20, 2, seed = 3)
  Y  <- pathInterpolate(X, n = 100, method = "step")
  expect_equal(Y[1, ],         X[1, ],        tolerance = TOL_STRICT)
  expect_equal(Y[nrow(Y), ],   X[nrow(X), ],  tolerance = TOL_STRICT)
})

test_that("explicit time grid is honored", {
  X <- makeRandomPath(10, 2, seed = 4)
  Y <- pathInterpolate(X, time = c(0, 0.5, 1))
  expect_equal(dim(Y), c(3L, 2L))
  expect_equal(Y[1, ], X[1, ],       tolerance = TOL_STRICT)
  expect_equal(Y[3, ], X[nrow(X), ], tolerance = TOL_STRICT)
})

test_that("explicit time at source-grid points reproduces source rows", {
  X   <- makeRandomPath(11, 2, seed = 5)
  src <- seq(0, 1, length.out = 11)
  Y   <- pathInterpolate(X, time = src)
  expect_equal(unname(Y), unname(X), tolerance = TOL_STRICT)
})

test_that("vector input is treated as a 1D path", {
  v <- cumsum(rnorm(20))
  Y <- pathInterpolate(v, n = 50)
  expect_equal(dim(Y), c(50L, 1L))
})

test_that("dimnames are stripped", {
  X <- matrix(rnorm(20), 10, 2,
              dimnames = list(paste0("r", 1:10), c("a", "b")))
  Y <- pathInterpolate(X, n = 25)
  expect_null(dimnames(Y))
})

test_that("signature is invariant under regridding a linear path", {
  # For piecewise-linear (in fact linear) paths the signature is a
  # function of the geometry, not the parametrization. Regridding
  # changes the parametrization but not the geometry.
  for (d in 1:3) {
    a <- runif(d, -1, 1); b <- runif(d, -1, 1)
    X <- makeLineND(a, b, n = 30)
    Y <- pathInterpolate(X, n = 73)
    for (N in 1:3) {
      expect_equal(unname(signature(X, depth = N)),
                   unname(signature(Y, depth = N)),
                   tolerance = TOL_STRICT,
                   info = sprintf("d=%d N=%d", d, N))
    }
  }
})

test_that("signature is invariant under regridding that preserves source knots", {
  # Signature invariance under reparametrization requires the target grid
  # to CONTAIN the source knots (otherwise the regridded path cuts
  # corners and the polygon changes). A finer grid alone is not enough.
  X   <- makeRandomPath(20, 2, seed = 6)
  src <- seq(0, 1, length.out = nrow(X))
  # Insert 9 extra points between each pair of source knots: 20 + 19*9 = 191.
  insert <- seq(0, 1, length.out = 10)[-c(1, 10)]
  target <- sort(unique(c(src, as.vector(outer(insert / (length(src) - 1),
                                               head(src, -1), `+`)))))
  Y <- pathInterpolate(X, time = target)
  expect_equal(unname(signature(X, depth = 3)),
               unname(signature(Y, depth = 3)),
               tolerance = TOL_DEFAULT)
})

test_that("regridding to a grid that misses source knots changes the signature", {
  # Documents the geometric reality: subdividing onto a target grid that
  # doesn't contain the source knots cuts corners off the polygon and
  # produces a different signature at level >= 2. Level 1 is unaffected
  # because it depends only on endpoints.
  X <- makeRandomPath(20, 2, seed = 6)
  Y <- pathInterpolate(X, n = 200)   # grid doesn't align with source knots
  sX <- signature(X, depth = 2)
  sY <- signature(Y, depth = 2)
  expect_equal(unname(sX[2:3]), unname(sY[2:3]), tolerance = TOL_STRICT) # level 1
  expect_false(isTRUE(all.equal(unname(sX[4:7]), unname(sY[4:7]),        # level 2
                                tolerance = TOL_DEFAULT)))
})

test_that("n = 1 returns the source basepoint", {
  X <- makeRandomPath(10, 2, seed = 7)
  Y <- pathInterpolate(X, n = 1)
  expect_equal(dim(Y), c(1L, 2L))
  expect_equal(unname(Y[1, ]), unname(X[1, ]), tolerance = TOL_STRICT)
})

test_that("invalid inputs are rejected", {
  X <- makeRandomPath(10, 2, seed = 8)
  expect_error(pathInterpolate(matrix(1, 1, 2), n = 10), "at least 2")
  expect_error(pathInterpolate(X),                       "exactly one")
  expect_error(pathInterpolate(X, n = 10, time = c(0, 1)), "not both")
  expect_error(pathInterpolate(X, n = 0),                "positive")
  expect_error(pathInterpolate(X, n = -3),               "positive")
  expect_error(pathInterpolate(X, n = 1.5),              "positive")
  expect_error(pathInterpolate(X, time = c(0, -0.1)),    "\\[0, 1\\]")
  expect_error(pathInterpolate(X, time = c(0, 1.5)),     "\\[0, 1\\]")
  expect_error(pathInterpolate(X, time = c(0.5, 0.2)),   "non-decreasing")
  expect_error(pathInterpolate(X, time = c(0, NA, 1)),   "finite")
  expect_error(pathInterpolate(X, time = numeric(0)),    "length at least 1")
  expect_error(pathInterpolate(X, n = 10, method = "spline"))
})

test_that("non-strict monotonic time is allowed (repeated times OK)", {
  # We required non-decreasing, not strictly increasing. Repeats produce
  # duplicate rows, which is harmless for the signature.
  X <- makeRandomPath(10, 2, seed = 9)
  Y <- pathInterpolate(X, time = c(0, 0.5, 0.5, 1))
  expect_equal(dim(Y), c(4L, 2L))
  expect_equal(Y[2, ], Y[3, ], tolerance = TOL_STRICT)
})
