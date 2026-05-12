test_that("empty case: enumerateWords(d, 0) returns just the empty word", {
  for (d in 1:5) {
    w <- enumerateWords(d, 0)
    expect_length(w, 1L)
    expect_identical(w[[1]], integer(0))
  }
})

test_that("size matches sum(d^(0:N))", {
  for (d in 1:5) for (N in 0:4) {
    expect_length(enumerateWords(d, N), as.integer(sum(d^(0:N))))
  }
})

test_that("d = 2, depth = 2 enumerates in level-major lex order", {
  w <- enumerateWords(2, 2)
  expect_identical(w[[1]], integer(0))
  expect_identical(w[[2]], 1L)
  expect_identical(w[[3]], 2L)
  expect_identical(w[[4]], c(1L, 1L))
  expect_identical(w[[5]], c(1L, 2L))
  expect_identical(w[[6]], c(2L, 1L))
  expect_identical(w[[7]], c(2L, 2L))
})

test_that("parent-index identity holds for every level-k word", {
  for (d in 2:4) for (N in 1:3) {
    w <- enumerateWords(d, N)
    for (i in 2:length(w)) {
      parentIdx <- ((i - 2L) %/% d) + 1L
      letter    <- ((i - 2L) %%  d) + 1L
      expect_identical(w[[i]], c(w[[parentIdx]], as.integer(letter)),
                       info = sprintf("d=%d N=%d i=%d", d, N, i))
    }
  }
})

test_that("level-1 words are exactly the singletons in order", {
  for (d in 1:6) {
    w <- enumerateWords(d, 1)
    expect_identical(lapply(2:length(w), function(i) w[[i]]),
                     lapply(seq_len(d), function(j) as.integer(j)))
  }
})

test_that("invalid input is rejected", {
  expect_error(enumerateWords(0, 1),   "positive integer")
  expect_error(enumerateWords(-1, 1),  "positive integer")
  expect_error(enumerateWords(1.5, 1), "positive integer")
  expect_error(enumerateWords(2, -1),  "non-negative integer")
  expect_error(enumerateWords(2, 1.5), "non-negative integer")
  expect_error(enumerateWords(c(2, 3), 1), "positive integer")
})
