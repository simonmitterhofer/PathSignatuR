test_that("standardBracket: singleton returns the letter", {
  expect_identical(PathSignatuR:::.standardBracket(c(1L)), 1L)
  expect_identical(PathSignatuR:::.standardBracket(c(3L)), 3L)
})

test_that("standardBracket: [1, 2] = list(1, 2)", {
  expect_equal(PathSignatuR:::.standardBracket(c(1L, 2L)),
               list(left = 1L, right = 2L))
})

test_that("standardBracket: [1, 1, 2] = [1, [1, 2]]", {
  expect_equal(PathSignatuR:::.standardBracket(c(1L, 1L, 2L)),
               list(left = 1L, right = list(left = 1L, right = 2L)))
})

test_that("standardBracket: [1, 2, 2] = [[1, 2], 2]", {
  expect_equal(PathSignatuR:::.standardBracket(c(1L, 2L, 2L)),
               list(left = list(left = 1L, right = 2L), right = 2L))
})

test_that("bracketToTensor: singleton i is the level-1 basis vector e_i", {
  v <- PathSignatuR:::.bracketToTensor(1L, d = 2, N = 2)
  expect_equal(v, c(0, 1, 0, 0, 0, 0, 0))  # position 2 = letter "1"
  v <- PathSignatuR:::.bracketToTensor(2L, d = 2, N = 2)
  expect_equal(v, c(0, 0, 1, 0, 0, 0, 0))
})

test_that("bracketToTensor: [1, 2] expands to 1(x)2 - 2(x)1", {
  v <- PathSignatuR:::.bracketToTensor(list(left = 1L, right = 2L),
                                       d = 2, N = 2)
  expect_equal(v, c(0, 0, 0, 0, 1, -1, 0))  # "1,2" = +1, "2,1" = -1
})

test_that("bracketToTensor: [1, [1, 2]] = 1(x)1(x)2 - 2*1(x)2(x)1 + 2(x)1(x)1", {
  b <- list(left = 1L, right = list(left = 1L, right = 2L))
  v <- PathSignatuR:::.bracketToTensor(b, d = 2, N = 3)
  expected <- numeric(15)
  expected[9]  <- 1                       # "1,1,2"
  expected[10] <- -2                      # "1,2,1"
  expected[12] <- 1                       # "2,1,1"
  expect_equal(v, expected, tolerance = TOL_STRICT)
})

test_that("lyndonProjectionMatrix has shape and rank as expected", {
  for (d in 2:3) for (N in 1:4) {
    M <- PathSignatuR:::.lyndonProjectionMatrix(d, N)
    L <- sum(d^(0:N))
    q <- length(lyndonWords(d, N))
    expect_equal(dim(M), c(L, q),
                 info = sprintf("d=%d N=%d shape", d, N))
    expect_equal(qr(M)$rank, q,
                 info = sprintf("d=%d N=%d rank", d, N))
  }
})

test_that("projection matrix has zero level-0 rows (brackets produce level >= 1)", {
  for (d in 2:3) for (N in 1:3) {
    M <- PathSignatuR:::.lyndonProjectionMatrix(d, N)
    expect_equal(M[1, ], rep(0, ncol(M)), tolerance = TOL_STRICT)
  }
})

test_that("projection matrix is block-diagonal in levels", {
  # Column for a Lyndon word of length k has support only at level k of
  # the tensor algebra. Check this for d = 2, N = 3.
  d <- 2; N <- 3
  M <- PathSignatuR:::.lyndonProjectionMatrix(d, N)
  lyn <- lyndonWords(d, N)
  # Level boundaries in enumerateWords order: levels 0, 1, 2, 3 sit at
  # rows 1, 2:3, 4:7, 8:15.
  levRows <- list(`0` = 1, `1` = 2:3, `2` = 4:7, `3` = 8:15)
  for (i in seq_along(lyn)) {
    k <- length(lyn[[i]])
    own <- levRows[[as.character(k)]]
    other <- setdiff(seq_len(nrow(M)), own)
    expect_equal(M[other, i], rep(0, length(other)),
                 tolerance = TOL_STRICT,
                 info = sprintf("Lyndon word #%d (length %d)", i, k))
  }
})
