test_that(".inferDim solves (d, N) -> length uniquely when constrained by depth", {
  expect_identical(PathSignatuR:::.inferDim(1L,  0L),  1L)
  expect_identical(PathSignatuR:::.inferDim(7L,  2L),  2L)
  expect_identical(PathSignatuR:::.inferDim(13L, 2L),  3L)
  expect_identical(PathSignatuR:::.inferDim(15L, 3L),  2L)
  expect_identical(PathSignatuR:::.inferDim(40L, 3L),  3L)
  expect_identical(PathSignatuR:::.inferDim(341L, 4L), 4L)
})

test_that(".inferDim rejects impossible lengths", {
  expect_error(PathSignatuR:::.inferDim(2L,  0L), "incompatible with depth 0")
  expect_error(PathSignatuR:::.inferDim(8L,  2L), "is not sum")
  expect_error(PathSignatuR:::.inferDim(14L, 2L), "is not sum")
})

test_that("depth = 0: tensor product reduces to scalar multiplication", {
  expect_equal(PathSignatuR:::.tensorProductRef(c(3), c(4), 0), 12)
})

test_that("identity tensor (1, 0, ...) is left and right unit", {
  for (d in 1:3) for (N in 1:3) {
    L <- sum(d^(0:N))
    e <- c(1, rep(0, L - 1))
    set.seed(d * 10 + N)
    a <- rnorm(L)
    expect_equal(PathSignatuR:::.tensorProductRef(a, e, N), a,
                 tolerance = TOL_STRICT,
                 info = sprintf("d=%d N=%d (right unit)", d, N))
    expect_equal(PathSignatuR:::.tensorProductRef(e, a, N), a,
                 tolerance = TOL_STRICT,
                 info = sprintf("d=%d N=%d (left unit)", d, N))
  }
})

test_that("tensor product is associative up to truncation", {
  set.seed(7)
  for (d in 2:3) for (N in 2:3) {
    L <- sum(d^(0:N))
    a <- rnorm(L); b <- rnorm(L); c <- rnorm(L)
    ab_c <- PathSignatuR:::.tensorProductRef(
      PathSignatuR:::.tensorProductRef(a, b, N), c, N)
    a_bc <- PathSignatuR:::.tensorProductRef(
      a, PathSignatuR:::.tensorProductRef(b, c, N), N)
    expect_equal(ab_c, a_bc, tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("tensor product is bilinear", {
  set.seed(11)
  d <- 2; N <- 3
  L <- sum(d^(0:N))
  a1 <- rnorm(L); a2 <- rnorm(L); b <- rnorm(L)
  alpha <- 0.7; beta <- -1.3
  lhs <- PathSignatuR:::.tensorProductRef(alpha * a1 + beta * a2, b, N)
  rhs <- alpha * PathSignatuR:::.tensorProductRef(a1, b, N) +
    beta  * PathSignatuR:::.tensorProductRef(a2, b, N)
  expect_equal(lhs, rhs, tolerance = TOL_DEFAULT)
})

test_that("Chen identity: signatures concatenate via tensor product", {
  # S(X1 . X2) = S(X1) (x) S(X2), where . is path concatenation.
  set.seed(13)
  for (d in 2:3) for (N in 2:3) {
    X  <- makeRandomPath(50, d, seed = d * 10 + N)
    mid <- 25L
    X1 <- X[1:mid, , drop = FALSE]
    X2 <- X[mid:nrow(X), , drop = FALSE]   # share join point
    s1    <- signature(X1, depth = N)
    s2    <- signature(X2, depth = N)
    sFull <- signature(X,  depth = N)
    sProd <- PathSignatuR:::.tensorProductRef(unname(s1), unname(s2), N)
    expect_equal(unname(sFull), sProd, tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that(".tensorProductRef validates inputs", {
  expect_error(PathSignatuR:::.tensorProductRef("a", c(1), 0),
               "must be numeric")
  expect_error(PathSignatuR:::.tensorProductRef(c(1, 2), c(1), 0),
               "equal length")
  expect_error(PathSignatuR:::.tensorProductRef(c(1, NA), c(1, 2), 1),
               "finite")
  expect_error(PathSignatuR:::.tensorProductRef(c(1, 2), c(1, 2), 0),
               "incompatible with depth 0")
})

test_that("exported tensorProduct matches the R reference across (d, depth)", {
  set.seed(101)
  for (d in 1:3) for (N in 0:3) {
    L <- sum(d^(0:N))
    a <- rnorm(L); b <- rnorm(L)
    expect_equal(unname(tensorProduct(a, b, N)),
                 PathSignatuR:::.tensorProductRef(a, b, N),
                 tolerance = TOL_STRICT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("tensorProduct output is named in canonical word order", {
  a <- c(1, 0.5, -0.2, 0.1, 0, 0.3, -0.4)
  b <- c(1, -0.1, 0.2, 0,  0.5, 0,  0.1)
  out <- tensorProduct(a, b, 2)
  expect_equal(names(out), c("", "1", "2", "1,1", "1,2", "2,1", "2,2"))
})

test_that("tensorProduct(s1, s2) recovers signature of concatenated path", {
  set.seed(17)
  for (d in 2:3) for (N in 2:3) {
    X   <- makeRandomPath(60, d, seed = d * 100 + N)
    mid <- 30L
    X1  <- X[1:mid,        , drop = FALSE]
    X2  <- X[mid:nrow(X),  , drop = FALSE]
    s1  <- signature(X1, depth = N)
    s2  <- signature(X2, depth = N)
    sP  <- tensorProduct(s1, s2, depth = N)
    sF  <- signature(X, depth = N)
    expect_equal(unname(sP), unname(sF), tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("workspace cache populates and is reused", {
  rm(list = ls(envir = PathSignatuR:::.workspaceCache, all.names = TRUE),
     envir = PathSignatuR:::.workspaceCache)
  expect_length(ls(envir = PathSignatuR:::.workspaceCache), 0L)

  a <- rnorm(7); b <- rnorm(7)
  tensorProduct(a, b, depth = 2)
  expect_setequal(ls(envir = PathSignatuR:::.workspaceCache), "2.2")

  # Second call same (d, depth): no new entry.
  tensorProduct(a + 1, b - 1, depth = 2)
  expect_setequal(ls(envir = PathSignatuR:::.workspaceCache), "2.2")

  # Different (d, depth): new entry.
  a3 <- rnorm(13); b3 <- rnorm(13)
  tensorProduct(a3, b3, depth = 2)
  expect_setequal(ls(envir = PathSignatuR:::.workspaceCache), c("2.2", "3.2"))
})

test_that("tensorProduct validates inputs", {
  expect_error(tensorProduct("a", c(1), 0),       "must be numeric")
  expect_error(tensorProduct(c(1, 2), c(1), 0),   "equal length")
  expect_error(tensorProduct(c(1, NA), c(1, 2), 1), "finite")
  expect_error(tensorProduct(c(1, 2), c(1, 2), 0), "incompatible with depth 0")
  expect_error(tensorProduct(c(1, 2, 3, 4, 5, 6, 7, 8), c(1:8), 2),
               "is not sum")
})
