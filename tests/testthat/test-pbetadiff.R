# Reference values: Python_reference/ref_distributions_mpmath.py (mpmath, 40
# digits, two integral representations), cross-checked with scipy.
ref_pbetadiff <- data.frame(
  q   = c(0, 0.1, -0.2, 0.05, 0.2, 0.001, 0.95, -0.95),
  a_t = c(5, 6, 2.5, 0.5, 7.5, 30.5, 20, 1),
  b_t = c(5, 3, 10.5, 12.5, 5.5, 9970.5, 1, 20),
  a_c = c(5, 5, 4.5, 0.5, 4.5, 10.5, 1, 20),
  b_c = c(5, 8, 8.5, 12.5, 8.5, 9990.5, 20, 1),
  upper = c(0.5, 0.81832837711327888, 0.60911710328926547,
            0.16845257184576472, 0.57601763350543618, 0.94619967208045943,
            0.27051374743310446, 0.72948625256689554)
)

test_that("pbetadiff() matches high-precision upper tail references", {
  r <- ref_pbetadiff
  for (i in seq_len(nrow(r))) {
    expect_equal(pbetadiff(r$q[i], r$a_t[i], r$b_t[i], r$a_c[i], r$b_c[i],
                           lower.tail = FALSE),
                 r$upper[i], tolerance = 1e-9, info = paste("case", i))
  }
})

test_that("pbetadiff() lower tail is the complement of the upper tail", {
  r <- ref_pbetadiff
  lower <- pbetadiff(r$q, r$a_t, r$b_t, r$a_c, r$b_c)
  expect_equal(lower, 1 - r$upper, tolerance = 1e-9)
})

test_that("pbetadiff() is vectorised with recycling", {
  q <- c(-0.1, 0, 0.1)
  vec <- pbetadiff(q, 6, 3, 5, 8, lower.tail = FALSE)
  one <- vapply(q, pbetadiff, numeric(1), a_t = 6, b_t = 3, a_c = 5, b_c = 8,
                lower.tail = FALSE)
  expect_identical(vec, one)
  expect_true(all(diff(vec) < 0))
})

test_that("pbetadiff() returns 0 and 1 outside the support", {
  expect_identical(pbetadiff(1, 2, 3, 4, 5, lower.tail = FALSE), 0)
  expect_identical(pbetadiff(-1, 2, 3, 4, 5, lower.tail = FALSE), 1)
})

test_that("pbetadiff() validates its arguments", {
  expect_error(pbetadiff(0, 0, 1, 1, 1), "'a_t' must be greater than 0.")
  expect_error(pbetadiff(c(0, 0.1), c(1, 2, 3), 1, 1, 1),
               "'q' must have length 1 or 3.")
  expect_error(pbetadiff(0, 1, 1, 1, 1, lower.tail = NA),
               "'lower.tail' must be TRUE or FALSE.")
})

test_that("pbetadiff() matches the finite-sum formula for integer a_t at q = 0", {
  # For integer a_t, P(pi_t > pi_c) is a finite sum (Python_reference/
  # check_closed_forms.py confirms it against numerical integration).
  grid <- expand.grid(a_t = 1:8, b_t = c(0.5, 1, 3.5, 12),
                      a_c = c(0.5, 1, 3.5, 12), b_c = c(0.5, 1, 3.5, 12))
  exact <- vapply(seq_len(nrow(grid)), function(k) {
    g <- grid[k, ]
    i <- seq(0, g$a_t - 1)
    sum(exp(lbeta(g$a_c + i, g$b_c + g$b_t) - log(g$b_t + i) -
              lbeta(1 + i, g$b_t) - lbeta(g$a_c, g$b_c)))
  }, numeric(1))
  got <- pbetadiff(0, grid$a_t, grid$b_t, grid$a_c, grid$b_c,
                   lower.tail = FALSE)
  expect_equal(got, exact, tolerance = 1e-9)
})

test_that("pbetadiff() agrees with integration over the treatment variable", {
  # Alternative representation: P(pi_t - pi_c > q) =
  #   int f_t(y) F_c(y - q) dy + [q < 0] P(pi_t > 1 + q)
  grid <- expand.grid(q = c(-0.6, -0.2, 0.05, 0.3, 0.7),
                      shapes = 1:3)
  shapes <- list(c(2.5, 4.5, 3.5, 6.5), c(0.5, 9.5, 1.5, 8.5),
                 c(15.5, 5.5, 10.5, 10.5))
  for (k in seq_len(nrow(grid))) {
    q <- grid$q[k]
    s <- shapes[[grid$shapes[k]]]
    alt <- stats::integrate(
      function(y) stats::dbeta(y, s[1], s[2]) * stats::pbeta(y - q, s[3], s[4]),
      lower = max(0, q), upper = min(1, 1 + q), rel.tol = 1e-11
    )$value
    if (q < 0) {
      alt <- alt + stats::pbeta(1 + q, s[1], s[2], lower.tail = FALSE)
    }
    expect_equal(pbetadiff(q, s[1], s[2], s[3], s[4], lower.tail = FALSE),
                 alt, tolerance = 1e-7, info = paste("case", k))
  }
})

test_that("pbetadiff() is a non-increasing function of q in [0, 1]", {
  q <- seq(-1.2, 1.2, length.out = 241)
  p <- pbetadiff(q, 3.5, 9.5, 6.5, 6.5, lower.tail = FALSE)
  expect_true(all(diff(p) <= 1e-10))
  expect_true(all(p >= 0 & p <= 1))
  expect_identical(p[q <= -1], rep(1, sum(q <= -1)))
  expect_identical(p[q >= 1], rep(0, sum(q >= 1)))
})

test_that("pbetadiff() handles long vectors", {
  q <- seq(-0.5, 0.5, length.out = 1200)
  p <- pbetadiff(q, 7.5, 5.5, 4.5, 8.5)
  expect_length(p, 1200L)
  expect_identical(p[c(1, 600, 1200)],
                   pbetadiff(q[c(1, 600, 1200)], 7.5, 5.5, 4.5, 8.5))
})

test_that("pbetadiff() agrees with Monte Carlo simulation", {
  skip_on_cran()
  set.seed(11)
  n <- 1e6
  d <- stats::rbeta(n, 7.5, 5.5) - stats::rbeta(n, 4.5, 8.5)
  for (q in c(-0.1, 0, 0.2)) {
    p <- pbetadiff(q, 7.5, 5.5, 4.5, 8.5, lower.tail = FALSE)
    expect_lt(abs(mean(d > q) - p), 4 * sqrt(p * (1 - p) / n))
  }
})
