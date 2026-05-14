test_that("depth = 0: shuffle reduces to scalar multiplication", {
  expect_equal(PathSignatuR:::.shuffleProductRef(c(3), c(4), 0), 12)
})

test_that("identity tensor is the shuffle unit", {
  for (d in 1:3) for (N in 1:3) {
    L <- sum(d^(0:N))
    e <- c(1, rep(0, L - 1))
    set.seed(d * 10 + N)
    a <- rnorm(L)
    expect_equal(PathSignatuR:::.shuffleProductRef(a, e, N), a,
                 tolerance = TOL_STRICT,
                 info = sprintf("d=%d N=%d (right)", d, N))
    expect_equal(PathSignatuR:::.shuffleProductRef(e, a, N), a,
                 tolerance = TOL_STRICT,
                 info = sprintf("d=%d N=%d (left)", d, N))
  }
})

test_that("shuffle product is commutative", {
  set.seed(11)
  for (d in 2:3) for (N in 1:3) {
    L <- sum(d^(0:N))
    a <- rnorm(L); b <- rnorm(L)
    expect_equal(PathSignatuR:::.shuffleProductRef(a, b, N),
                 PathSignatuR:::.shuffleProductRef(b, a, N),
                 tolerance = TOL_STRICT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("shuffle product is associative", {
  set.seed(13)
  for (d in 2:3) for (N in 2:3) {
    L <- sum(d^(0:N))
    a <- rnorm(L); b <- rnorm(L); c <- rnorm(L)
    ab_c <- PathSignatuR:::.shuffleProductRef(
      PathSignatuR:::.shuffleProductRef(a, b, N), c, N)
    a_bc <- PathSignatuR:::.shuffleProductRef(
      a, PathSignatuR:::.shuffleProductRef(b, c, N), N)
    expect_equal(ab_c, a_bc, tolerance = TOL_DEFAULT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("shuffle product is bilinear", {
  set.seed(17)
  d <- 2; N <- 3
  L <- sum(d^(0:N))
  a1 <- rnorm(L); a2 <- rnorm(L); b <- rnorm(L)
  alpha <- 0.7; beta <- -1.3
  lhs <- PathSignatuR:::.shuffleProductRef(alpha * a1 + beta * a2, b, N)
  rhs <- alpha * PathSignatuR:::.shuffleProductRef(a1, b, N) +
    beta  * PathSignatuR:::.shuffleProductRef(a2, b, N)
  expect_equal(lhs, rhs, tolerance = TOL_DEFAULT)
})

test_that("shuffle of e_u and e_v: explicit cases", {
  # d = 2, depth = 2. Words: "", "1", "2", "1,1", "1,2", "2,1", "2,2".
  L <- 7
  e <- function(i) { v <- numeric(L); v[i] <- 1; v }

  # e_"1" ⧢ e_"2" = e_"1,2" + e_"2,1"
  sh <- PathSignatuR:::.shuffleProductRef(e(2), e(3), 2)
  expect_equal(sh, c(0, 0, 0, 0, 1, 1, 0), tolerance = TOL_STRICT)

  # e_"1" ⧢ e_"1" = 2 * e_"1,1"
  sh <- PathSignatuR:::.shuffleProductRef(e(2), e(2), 2)
  expect_equal(sh, c(0, 0, 0, 2, 0, 0, 0), tolerance = TOL_STRICT)

  # e_∅ ⧢ e_"1,2" = e_"1,2"
  sh <- PathSignatuR:::.shuffleProductRef(e(1), e(5), 2)
  expect_equal(sh, c(0, 0, 0, 0, 1, 0, 0), tolerance = TOL_STRICT)
})

test_that("shuffle identity: <S, u> * <S, v> = <S, u ⧢ v> for signatures", {
  # The headline algebraic property. Requires depth >= |u| + |v|.
  set.seed(31)
  for (d in 2:3) {
    X <- makeRandomPath(40, d, seed = d * 100)
    N <- 4
    s <- signature(X, depth = N)
    L <- length(s)

    # Test several (u, v) pairs with |u| + |v| <= N.
    words <- enumerateWords(d, N)
    pairs <- list(c(2, 3),          # u="1", v="2": sums to length 2
                  c(2, 5),          # u="1", v="1,1": sums to length 3
                  c(3, 4),          # u="2", v="1,1": sums to length 3
                  c(5, 6))          # u="1,1", v="1,2": sums to length 4
    for (p in pairs) {
      uPos <- p[1]; vPos <- p[2]
      uLen <- length(words[[uPos]])
      vLen <- length(words[[vPos]])
      if (uLen + vLen > N) next
      eu <- numeric(L); eu[uPos] <- 1
      ev <- numeric(L); ev[vPos] <- 1
      sh <- PathSignatuR:::.shuffleProductRef(eu, ev, N)
      lhs <- sum(sh * unname(s))
      rhs <- unname(s[uPos]) * unname(s[vPos])
      expect_equal(lhs, rhs, tolerance = TOL_DEFAULT,
                   info = sprintf("d=%d uPos=%d vPos=%d", d, uPos, vPos))
    }
  }
})

test_that(".shuffleProductRef validates inputs", {
  expect_error(PathSignatuR:::.shuffleProductRef("a", c(1), 0),
               "must be numeric")
  expect_error(PathSignatuR:::.shuffleProductRef(c(1, 2), c(1), 0),
               "equal length")
  expect_error(PathSignatuR:::.shuffleProductRef(c(1, NA), c(1, 2), 1),
               "finite")
})

test_that("exported shuffleProduct matches the R reference", {
  set.seed(101)
  for (d in 1:3) for (N in 0:3) {
    L <- sum(d^(0:N))
    a <- rnorm(L); b <- rnorm(L)
    expect_equal(unname(shuffleProduct(a, b, N)),
                 PathSignatuR:::.shuffleProductRef(a, b, N),
                 tolerance = TOL_STRICT,
                 info = sprintf("d=%d N=%d", d, N))
  }
})

test_that("shuffleProduct output is named in canonical word order", {
  a <- c(1, 0.5, -0.2, 0.1, 0, 0.3, -0.4)
  b <- c(1, -0.1, 0.2, 0, 0.5, 0, 0.1)
  out <- shuffleProduct(a, b, 2)
  expect_equal(names(out), c("", "1", "2", "1,1", "1,2", "2,1", "2,2"))
})

test_that("shuffle workspace cache populates and is reused", {
  rm(list = ls(envir = PathSignatuR:::.shuffleCache, all.names = TRUE),
     envir = PathSignatuR:::.shuffleCache)
  a <- rnorm(7); b <- rnorm(7)
  shuffleProduct(a, b, 2)
  expect_setequal(ls(envir = PathSignatuR:::.shuffleCache), "2.2")
  shuffleProduct(a + 1, b - 1, 2)
  expect_setequal(ls(envir = PathSignatuR:::.shuffleCache), "2.2")
})
