test_that("output shape is T x p with column names", {
  X  <- matrix(rnorm(20), 10, 2)
  lp <- logSignaturePath(X, depth = 2)
  expect_equal(dim(lp), c(10L, 7L))
  expect_equal(colnames(lp), c("", "1", "2", "1,1", "1,2", "2,1", "2,2"))
  expect_null(rownames(lp))
})

test_that("row 1 is zero (log of identity)", {
  X  <- matrix(rnorm(40), 20, 2)
  lp <- logSignaturePath(X, depth = 3)
  expect_equal(unname(lp[1, ]), rep(0, ncol(lp)))
})

test_that("last row matches logSignature() on the full path", {
  set.seed(11)
  for (d in 1:3) for (N in 0:3) {
    X  <- makeRandomPath(30, d, seed = d * 10 + N)
    lp <- logSignaturePath(X, depth = N)
    ls <- logSignature(X, depth = N)
    expect_equal(unname(lp[nrow(lp), ]), unname(ls),
                 tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("row t matches logSignature() on the prefix X[1:t, ]", {
  set.seed(7)
  X  <- makeRandomPath(20, 2, seed = 5)
  lp <- logSignaturePath(X, depth = 2)
  for (t in 2:nrow(X)) {
    expect_equal(unname(lp[t, ]),
                 unname(logSignature(X[1:t, , drop = FALSE], depth = 2)),
                 tolerance = TOL_DEFAULT,
                 info = sprintf("t=%d", t))
  }
})

test_that("1D linear path: only level-1 entry is nonzero, equals x_t - x_1", {
  x  <- seq(0, 2, length.out = 50)
  lp <- logSignaturePath(matrix(x, ncol = 1), depth = 4)
  for (t in 1:50) {
    expect_equal(unname(lp[t, 1]), 0, tolerance = TOL_STRICT,
                 info = sprintf("t=%d empty", t))
    expect_equal(unname(lp[t, 2]), x[t] - x[1], tolerance = TOL_DEFAULT,
                 info = sprintf("t=%d level-1", t))
    if (ncol(lp) > 2) {
      expect_equal(unname(lp[t, -(1:2)]), rep(0, ncol(lp) - 2L),
                   tolerance = TOL_DEFAULT,
                   info = sprintf("t=%d higher", t))
    }
  }
})

test_that("exp/log roundtrip across the running signature", {
  set.seed(3)
  X  <- makeRandomPath(20, 2, seed = 5)
  lp <- logSignaturePath(X, depth = 3)
  sp <- signaturePath(X, depth = 3)
  for (t in seq_len(nrow(X))) {
    sRecon <- PathSignatuR:::.tensorExpRef(unname(lp[t, ]), 3)
    expect_equal(unname(sp[t, ]), sRecon, tolerance = TOL_DEFAULT,
                 info = sprintf("t=%d", t))
  }
})

test_that("vector input gives T x p", {
  v <- rnorm(15)
  expect_equal(dim(logSignaturePath(v, depth = 2)), c(15L, 3L))
})

test_that("includeLevelZero = FALSE drops the first column", {
  X    <- matrix(rnorm(20), 10, 2)
  full <- logSignaturePath(X, depth = 2, includeLevelZero = TRUE)
  drop <- logSignaturePath(X, depth = 2, includeLevelZero = FALSE)
  expect_equal(ncol(drop), ncol(full) - 1L)
  expect_equal(unname(drop), unname(full[, -1L, drop = FALSE]))
  expect_false("" %in% colnames(drop))
})

test_that("sep controls column-name format", {
  X <- matrix(rnorm(20), 10, 2)
  expect_equal(colnames(logSignaturePath(X, 2, sep = ""))[4:6],
               c("11", "12", "21"))
})

test_that("depth = 0 returns a T x 1 column of zeros", {
  X  <- matrix(rnorm(20), 10, 2)
  lp <- logSignaturePath(X, depth = 0)
  expect_equal(dim(lp), c(10L, 1L))
  expect_equal(unname(lp), matrix(0, 10, 1))
})

test_that("T = 1 returns a single zero row at any depth", {
  X  <- matrix(c(0.3, -0.1), 1, 2)
  lp <- logSignaturePath(X, depth = 3)
  expect_equal(dim(lp), c(1L, 15L))
  expect_equal(unname(lp[1, ]), rep(0, 15))
})

test_that("workspace cache is reused across calls", {
  rm(list = ls(envir = PathSignatuR:::.workspaceCache, all.names = TRUE),
     envir = PathSignatuR:::.workspaceCache)
  X <- makeRandomPath(15, 2, seed = 1)
  logSignaturePath(X, depth = 2)
  expect_setequal(ls(envir = PathSignatuR:::.workspaceCache), "2.2")
  logSignaturePath(X, depth = 2)
  expect_setequal(ls(envir = PathSignatuR:::.workspaceCache), "2.2")
})

test_that("invalid input errors cleanly", {
  expect_error(logSignaturePath(NULL, 2),                   "NULL")
  expect_error(logSignaturePath("not numeric", 2),          "numeric")
  expect_error(logSignaturePath(matrix(NA_real_, 3, 2), 2), "finite")
  expect_error(logSignaturePath(matrix(rnorm(6), 3, 2), -1),  "non-negative")
  expect_error(logSignaturePath(matrix(rnorm(6), 3, 2), 1.5), "non-negative")
})
