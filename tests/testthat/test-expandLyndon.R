test_that("expandLyndon: at depth >= 1 returns named tensor-basis vector", {
  X    <- makeRandomPath(20, 2, seed = 5)
  beta <- logSignatureLyndon(X, depth = 3)
  out  <- expandLyndon(beta, dim = 2, depth = 3)
  expect_length(out, 15L)
  expect_equal(names(out),
               c("", "1", "2", "1,1", "1,2", "2,1", "2,2",
                 "1,1,1", "1,1,2", "1,2,1", "1,2,2",
                 "2,1,1", "2,1,2", "2,2,1", "2,2,2"))
})

test_that("expandLyndon: level-0 entry is always zero", {
  beta <- runif(length(lyndonWords(3, 4)), -1, 1)
  out  <- expandLyndon(beta, dim = 3, depth = 4)
  expect_equal(unname(out[1]), 0, tolerance = TOL_STRICT)
})

test_that("expandLyndon: depth = 0 with empty beta returns c('' = 0)", {
  out <- expandLyndon(numeric(0), dim = 3, depth = 0)
  expect_length(out, 1L)
  expect_equal(unname(out), 0)
  expect_equal(names(out), "")
})

test_that("expandLyndon: round-trip from logSignatureLyndon to logSignature", {
  set.seed(13)
  for (d in 2:3) for (N in 1:4) {
    X    <- makeRandomPath(20, d, seed = d * 50 + N)
    beta <- logSignatureLyndon(X, depth = N)
    expect_equal(unname(expandLyndon(beta, d, N)),
                 unname(logSignature(X, depth = N)),
                 tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("expandLyndon: sep controls tensor-word naming", {
  beta <- c(0.5, -0.3, 1.2, 0.7, -0.4)
  out  <- expandLyndon(beta, dim = 2, depth = 3, sep = "-")
  expect_true("1-1-2" %in% names(out))
})

test_that("expandLyndon: length mismatch errors", {
  expect_error(expandLyndon(c(1, 2), dim = 3, depth = 2),
               "length 6")
  expect_error(expandLyndon(c(1), dim = 2, depth = 0),
               "length 0")
})

test_that("expandLyndon: invalid input errors", {
  expect_error(expandLyndon("a", 2, 1),         "numeric")
  expect_error(expandLyndon(c(1, NA), 2, 1),    "finite")
  expect_error(expandLyndon(c(1, 2), 0, 1),     "positive integer")
  expect_error(expandLyndon(c(1, 2), 2, -1),    "non-negative")
  expect_error(expandLyndon(c(1, 2), 2, 1, sep = c(",", ";")),
               "single string")
})
