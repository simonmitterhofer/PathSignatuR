test_that("counts match Witt's formula across (d, N) grid", {
  cases <- list(
    c(1, 1, 1),  c(1, 5, 1),
    c(2, 1, 2),  c(2, 2, 3),  c(2, 3, 5),  c(2, 4, 8),  c(2, 5, 14),
    c(3, 1, 3),  c(3, 2, 6),  c(3, 3, 14), c(3, 4, 32),
    c(4, 1, 4),  c(4, 2, 10), c(4, 3, 30), c(4, 4, 90)
  )
  for (cs in cases) {
    expect_equal(length(lyndonWords(cs[1], cs[2])), cs[3],
                 info = sprintf("d=%d N=%d", cs[1], cs[2]))
  }
})

test_that("depth = 0 returns empty list", {
  expect_identical(lyndonWords(3, 0), list())
})

test_that("output is in canonical order (by length, then lex)", {
  w <- lyndonWords(2, 3)
  expect_identical(w, list(1L, 2L, c(1L, 2L), c(1L, 1L, 2L), c(1L, 2L, 2L)))
})

test_that("d = 2 N = 4 enumeration is exactly the 8 expected words", {
  w <- lyndonWords(2, 4)
  expect_identical(w, list(
    1L, 2L,
    c(1L, 2L),
    c(1L, 1L, 2L), c(1L, 2L, 2L),
    c(1L, 1L, 1L, 2L), c(1L, 1L, 2L, 2L), c(1L, 2L, 2L, 2L)
  ))
})

test_that("every Lyndon word is strictly less than its rotations", {
  for (d in 2:3) for (N in 1:4) {
    w <- lyndonWords(d, N)
    for (word in w) {
      n <- length(word)
      if (n == 1L) next
      for (i in 1:(n - 1L)) {
        rot <- c(word[(i + 1L):n], word[1:i])
        expect_true(PathSignatuR:::.lexLess(word, rot),
                    info = sprintf("d=%d word=(%s) rot=(%s)",
                                   d,
                                   paste(word, collapse = ","),
                                   paste(rot,  collapse = ",")))
      }
    }
  }
})

test_that(".isLyndon and .lexLess basics", {
  expect_true(PathSignatuR:::.isLyndon(c(1L, 2L)))
  expect_false(PathSignatuR:::.isLyndon(c(2L, 1L)))
  expect_false(PathSignatuR:::.isLyndon(c(1L, 1L)))      # periodic
  expect_true(PathSignatuR:::.isLyndon(c(1L, 1L, 2L)))
  expect_false(PathSignatuR:::.isLyndon(c(1L, 2L, 1L, 2L)))  # period 2

  expect_true(PathSignatuR:::.lexLess(c(1L, 1L), c(1L, 2L)))
  expect_false(PathSignatuR:::.lexLess(c(1L, 2L), c(1L, 2L)))  # equal
  expect_false(PathSignatuR:::.lexLess(c(1L, 2L), c(1L, 1L)))
})

test_that("invalid input is rejected", {
  expect_error(lyndonWords(0,  1), "positive integer")
  expect_error(lyndonWords(-1, 1), "positive integer")
  expect_error(lyndonWords(2, -1), "non-negative integer")
  expect_error(lyndonWords(2, 1.5), "non-negative integer")
})
