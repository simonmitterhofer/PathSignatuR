test_that("output shape is T x q with column names", {
  X  <- matrix(rnorm(20), 10, 2)
  bp <- logSignaturePathLyndon(X, depth = 3)
  expect_equal(dim(bp), c(10L, 5L))
  expect_equal(colnames(bp), c("1", "2", "1,2", "1,1,2", "1,2,2"))
  expect_null(rownames(bp))
})

test_that("row 1 is zero (log of identity)", {
  X  <- matrix(rnorm(40), 20, 2)
  bp <- logSignaturePathLyndon(X, depth = 3)
  expect_equal(unname(bp[1, ]), rep(0, ncol(bp)))
})

test_that("last row matches logSignatureLyndon on the full path", {
  set.seed(11)
  for (d in 2:3) for (N in 1:3) {
    X    <- makeRandomPath(30, d, seed = d * 10 + N)
    bp   <- logSignaturePathLyndon(X, depth = N)
    beta <- logSignatureLyndon(X, depth = N)
    expect_equal(unname(bp[nrow(bp), ]), unname(beta),
                 tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("row t matches logSignatureLyndon on the prefix X[1:t, ]", {
  set.seed(7)
  X  <- makeRandomPath(15, 2, seed = 5)
  bp <- logSignaturePathLyndon(X, depth = 3)
  for (t in 2:nrow(X)) {
    expect_equal(unname(bp[t, ]),
                 unname(logSignatureLyndon(X[1:t, , drop = FALSE], depth = 3)),
                 tolerance = TOL_DEFAULT,
                 info = sprintf("t=%d", t))
  }
})

test_that("depth = 0 returns T x 0 matrix", {
  X  <- matrix(rnorm(10), 5, 2)
  bp <- logSignaturePathLyndon(X, depth = 0)
  expect_equal(dim(bp), c(5L, 0L))
})

test_that("sep controls column-name format", {
  X <- matrix(rnorm(20), 10, 2)
  expect_equal(colnames(logSignaturePathLyndon(X, 3, sep = ""))[4:5],
               c("112", "122"))
})

test_that("vector input gives T x q", {
  v  <- cumsum(rnorm(20))
  bp <- logSignaturePathLyndon(v, depth = 3)
  expect_equal(dim(bp), c(20L, 1L))
})

test_that("invalid input errors cleanly", {
  expect_error(logSignaturePathLyndon(NULL, 2),                 "NULL")
  expect_error(logSignaturePathLyndon("not numeric", 2),        "numeric")
  expect_error(logSignaturePathLyndon(matrix(NA_real_, 3, 2), 2), "finite")
  expect_error(logSignaturePathLyndon(matrix(rnorm(6), 3, 2), -1), "non-negative")
})
