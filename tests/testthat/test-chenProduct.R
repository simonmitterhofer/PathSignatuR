test_that("chenProduct is identical to tensorProduct", {
  set.seed(1)
  for (d in 1:3) for (N in 0:3) {
    L <- sum(d^(0:N))
    a <- rnorm(L); b <- rnorm(L)
    expect_identical(chenProduct(a, b, N), tensorProduct(a, b, N))
  }
})

test_that("chenProduct(s1, s2) recovers signature of concatenated path", {
  set.seed(7)
  X   <- makeRandomPath(40, 2, seed = 1)
  mid <- 20L
  s1  <- signature(X[1:mid, , drop = FALSE], depth = 3)
  s2  <- signature(X[mid:nrow(X), , drop = FALSE], depth = 3)
  sP  <- chenProduct(s1, s2, depth = 3)
  expect_equal(unname(sP), unname(signature(X, depth = 3)),
               tolerance = TOL_DEFAULT)
})
