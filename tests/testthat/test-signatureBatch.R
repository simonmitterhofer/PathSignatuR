test_that("list input produces nPaths x p matrix with column names", {
  paths <- list(matrix(rnorm(20), 10, 2),
                matrix(rnorm(30), 15, 2))
  out <- signatureBatch(paths, depth = 2)
  expect_equal(dim(out), c(2L, 7L))
  expect_equal(colnames(out), c("", "1", "2", "1,1", "1,2", "2,1", "2,2"))
})

test_that("list names propagate to row names", {
  paths <- list(a = matrix(rnorm(20), 10, 2),
                b = matrix(rnorm(30), 15, 2))
  expect_equal(rownames(signatureBatch(paths, 2)), c("a", "b"))
})

test_that("unnamed list yields no row names", {
  paths <- list(matrix(rnorm(20), 10, 2),
                matrix(rnorm(30), 15, 2))
  expect_null(rownames(signatureBatch(paths, 2)))
})

test_that("3D array dispatches as a batch with uniform T", {
  arr <- array(rnorm(60), dim = c(10, 2, 3))
  out <- signatureBatch(arr, depth = 2)
  expect_equal(dim(out), c(3L, 7L))
})

test_that("3D array dimnames[[3]] propagate to row names", {
  arr <- array(rnorm(60), dim = c(10, 2, 3),
               dimnames = list(NULL, NULL, c("p1", "p2", "p3")))
  expect_equal(rownames(signatureBatch(arr, 2)), c("p1", "p2", "p3"))
})

test_that("single matrix input gives a 1 x p matrix", {
  X <- matrix(rnorm(20), 10, 2)
  expect_equal(dim(signatureBatch(X, depth = 2)), c(1L, 7L))
})

test_that("single vector input gives a 1 x p matrix (1D path)", {
  v <- rnorm(10)
  out <- signatureBatch(v, depth = 3)
  expect_equal(dim(out), c(1L, 4L))
})

test_that("each row equals signature() on the corresponding path", {
  set.seed(1)
  paths <- lapply(1:5, function(i) makeRandomPath(20 + i, 2, seed = i))
  batch <- signatureBatch(paths, depth = 3)
  for (i in seq_along(paths)) {
    expect_equal(unname(batch[i, ]),
                 unname(signature(paths[[i]], depth = 3)),
                 tolerance = TOL_STRICT,
                 info = sprintf("path %d", i))
  }
})

test_that("ragged T is allowed", {
  paths <- list(makeRandomPath(10, 2, seed = 1),
                makeRandomPath(50, 2, seed = 2),
                makeRandomPath(3,  2, seed = 3))
  expect_silent(signatureBatch(paths, depth = 2))
})

test_that("mismatched d in a list errors", {
  paths <- list(matrix(rnorm(20), 10, 2),
                matrix(rnorm(30), 10, 3))
  expect_error(signatureBatch(paths, depth = 2), "expected 2")
})

test_that("includeLevelZero = FALSE drops the first column", {
  paths <- list(matrix(rnorm(20), 10, 2),
                matrix(rnorm(30), 15, 2))
  full <- signatureBatch(paths, depth = 2, includeLevelZero = TRUE)
  drop <- signatureBatch(paths, depth = 2, includeLevelZero = FALSE)
  expect_equal(ncol(drop), ncol(full) - 1L)
  expect_equal(unname(drop), unname(full[, -1L, drop = FALSE]))
  expect_false("" %in% colnames(drop))
})

test_that("sep controls column-name format", {
  paths <- list(matrix(rnorm(20), 10, 2))
  expect_equal(colnames(signatureBatch(paths, 2, sep = ""))[4:6],
               c("11", "12", "21"))
})

test_that("depth = 0 gives nPaths x 1 column of ones", {
  paths <- list(matrix(rnorm(20), 10, 2),
                matrix(rnorm(30), 15, 2))
  out <- signatureBatch(paths, depth = 0)
  expect_equal(dim(out), c(2L, 1L))
  expect_equal(unname(out), matrix(1, 2, 1))
})

test_that("invalid input errors cleanly", {
  expect_error(signatureBatch(NULL, 2), "NULL")
  expect_error(signatureBatch(list(), 2), "at least one")
  expect_error(signatureBatch(list("not numeric"), 2), "path 1:")
  expect_error(signatureBatch(list(matrix(NA_real_, 3, 2)), 2),
               "path 1:.*finite")
})
