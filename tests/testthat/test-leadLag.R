test_that("1D path: output matches the FHL staircase exactly", {
  expected <- matrix(c(
    1, 1,
    2, 1,
    2, 2,
    3, 2,
    3, 3
  ), 5, 2, byrow = TRUE)
  colnames(expected) <- c("lead_1", "lag_1")
  expect_equal(leadLag(c(1, 2, 3)), expected)
})

test_that("output shape is (2T - 1) x (2d)", {
  for (T in c(1, 2, 5, 20)) for (d in 1:4) {
    X <- makeRandomPath(T, d, seed = 10L * T + d)
    expect_equal(dim(leadLag(X)), c(2L * T - 1L, 2L * d),
                 info = sprintf("T=%d d=%d", T, d))
  }
})

test_that("column names are interleaved lead_i, lag_i", {
  X <- matrix(rnorm(30), 10, 3)
  expect_equal(colnames(leadLag(X)),
               c("lead_1", "lag_1", "lead_2", "lag_2", "lead_3", "lag_3"))
})

test_that("row names are stripped", {
  X <- matrix(rnorm(6), 3, 2, dimnames = list(c("a", "b", "c"), NULL))
  expect_null(rownames(leadLag(X)))
})

test_that("vector input is treated as a 1D path", {
  v  <- c(0.5, -0.2, 1.1)
  LL <- leadLag(v)
  expect_equal(dim(LL), c(5L, 2L))
  expect_equal(colnames(LL), c("lead_1", "lag_1"))
})

test_that("T = 1 returns a single row with all columns equal to the input", {
  X  <- matrix(c(0.3, -0.1), 1, 2)
  LL <- leadLag(X)
  expect_equal(dim(LL), c(1L, 4L))
  expect_equal(unname(LL[1, ]), c(0.3, 0.3, -0.1, -0.1))
})

test_that("lead column for channel i contains exactly the values of X[, i]", {
  X  <- makeRandomPath(15, 3, seed = 5)
  LL <- leadLag(X)
  for (i in 1:3) {
    expect_setequal(unique(LL[, paste0("lead_", i)]), X[, i])
    expect_setequal(unique(LL[, paste0("lag_",  i)]), X[, i])
  }
})

test_that("each channel's lead and lag start at X[1, i]", {
  X  <- makeRandomPath(10, 2, seed = 8)
  LL <- leadLag(X)
  expect_equal(unname(LL[1, "lead_1"]), X[1, 1])
  expect_equal(unname(LL[1, "lag_1"]),  X[1, 1])
  expect_equal(unname(LL[1, "lead_2"]), X[1, 2])
  expect_equal(unname(LL[1, "lag_2"]),  X[1, 2])
})

test_that("each channel's lead and lag end at X[T, i]", {
  X  <- makeRandomPath(10, 2, seed = 8)
  LL <- leadLag(X)
  expect_equal(unname(LL[nrow(LL), "lead_1"]), X[nrow(X), 1])
  expect_equal(unname(LL[nrow(LL), "lag_1"]),  X[nrow(X), 1])
  expect_equal(unname(LL[nrow(LL), "lead_2"]), X[nrow(X), 2])
  expect_equal(unname(LL[nrow(LL), "lag_2"]),  X[nrow(X), 2])
})

test_that("FHL identity: S^{lead, lag} - S^{lag, lead} = QV (1D)", {
  set.seed(1)
  x  <- cumsum(rnorm(100))
  qv <- sum(diff(x)^2)
  s  <- signature(leadLag(x), depth = 2)
  expect_equal(unname(s["1,2"] - s["2,1"]), qv, tolerance = TOL_STRICT)
})

test_that("FHL identity holds across multiple 1D paths and seeds", {
  for (seed in 1:5) {
    set.seed(seed)
    x  <- cumsum(rnorm(50 + 10 * seed))
    qv <- sum(diff(x)^2)
    s  <- signature(leadLag(x), depth = 2)
    expect_equal(unname(s["1,2"] - s["2,1"]), qv,
                 tolerance = TOL_STRICT,
                 info = sprintf("seed=%d", seed))
  }
})

test_that("FHL identity holds per-channel in a multi-d path", {
  # For each channel i, S^{lead_i, lag_i} - S^{lag_i, lead_i} = QV(X[, i]).
  # In integer-letter naming after the transform, channel i's lead is
  # letter (2i - 1) and lag is letter (2i).
  set.seed(2)
  d <- 3
  X <- matrix(cumsum(rnorm(50 * d)), 50, d)
  s <- signature(leadLag(X), depth = 2)
  for (i in seq_len(d)) {
    leadLetter <- 2L * i - 1L
    lagLetter  <- 2L * i
    keyLL <- paste(leadLetter, lagLetter, sep = ",")
    keyLLrev <- paste(lagLetter, leadLetter, sep = ",")
    qv <- sum(diff(X[, i])^2)
    expect_equal(unname(s[keyLL] - s[keyLLrev]), qv,
                 tolerance = TOL_STRICT,
                 info = sprintf("channel %d", i))
  }
})

test_that("invalid input is rejected", {
  expect_error(leadLag(NULL),                   "NULL")
  expect_error(leadLag("not numeric"),          "numeric")
  expect_error(leadLag(matrix(NA_real_, 3, 2)), "finite")
})
