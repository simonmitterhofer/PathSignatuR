test_that("a (x) tensorInverse(a) = e (identity)", {
  set.seed(1)
  for (d in 1:3) for (N in 0:3) {
    L <- sum(d^(0:N))
    a <- c(1, rnorm(L - 1))   # level-0 = 1
    if (L == 1L) a <- 1
    aInv <- tensorInverse(a, N)
    prod <- tensorProduct(a, aInv, N)
    e    <- c(1, rep(0, L - 1))
    expect_equal(unname(prod), e, tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d N=%d (right)", d, N))
    prodL <- tensorProduct(aInv, a, N)
    expect_equal(unname(prodL), e, tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d N=%d (left)", d, N))
  }
})

test_that("tensorInverse is an involution: (a^{-1})^{-1} = a", {
  set.seed(2)
  for (d in 2:3) for (N in 1:3) {
    L <- sum(d^(0:N))
    a <- c(1, rnorm(L - 1))
    expect_equal(unname(tensorInverse(tensorInverse(a, N), N)),
                 a, tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("identity tensor inverts to itself", {
  for (d in 1:3) for (N in 0:3) {
    L <- sum(d^(0:N))
    e <- c(1, rep(0, L - 1))
    expect_equal(unname(tensorInverse(e, N)), e,
                 tolerance = TOL_STRICT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("depth = 0 returns the trivial inverse (1)", {
  expect_equal(unname(tensorInverse(c(1), 0)), 1)
})

test_that("signature reversal identity: sig(X_rev) = sig(X)^{-1}", {
  set.seed(7)
  for (d in 2:3) for (N in 2:3) {
    X    <- makeRandomPath(40, d, seed = d * 10 + N)
    Xrev <- X[nrow(X):1, , drop = FALSE]
    s    <- signature(X,    depth = N)
    sRev <- signature(Xrev, depth = N)
    expect_equal(unname(sRev),
                 unname(tensorInverse(s, N)),
                 tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("output is named in canonical word order", {
  a <- c(1, 0.5, -0.2, 0.1, 0, 0.3, -0.4)
  out <- tensorInverse(a, 2)
  expect_equal(names(out), c("", "1", "2", "1,1", "1,2", "2,1", "2,2"))
})

test_that("invalid inputs are rejected", {
  expect_error(tensorInverse("a", 0),               "must be numeric")
  expect_error(tensorInverse(c(1, NA), 1),          "finite")
  expect_error(tensorInverse(c(0.5, 1, 2), 1),      "level-0 coefficient")
  expect_error(tensorInverse(c(1, 0, 0, 0, 0), 2),  "is not sum")
})
