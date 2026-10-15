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

# Stress cases: Python_reference/ref_pbetadiff_stress_mpmath.py (mpmath, 30
# digits, tanh-sinh quadrature with the singular points of the densities at
# the ends of the range of integration). Python_reference/
# select_pbetadiff_stress.py selects every tenth case of the grid and all
# cases in which the integration without the change of variables near the
# singular endpoints had an absolute error above 1e-9 (17 cases).
stress_pbetadiff <- data.frame(
  case = c(0, 5, 10, 14, 20, 30, 32, 37, 40, 41, 50, 60,
           70, 73, 76, 80, 90, 98, 100, 110, 113, 117, 119, 120,
           123, 124, 130, 139, 140, 150, 160, 170, 180, 190, 200, 210,
           220, 230, 240, 250),
  q = c(-0.9, -0.9, -0.9, -0.9, -0.9, -0.3, -0.3, -0.3, -0.3, -0.3,
        -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -1e-06, -1e-06, -1e-06, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-06,
        1e-06, 1e-06, 0.05, 0.05, 0.05, 0.05, 0.05, 0.3, 0.3, 0.9),
  a_t = c(0.3, 0.8, 1, 3.5, 12, 0.3, 0.5, 1.5, 3.5, 3.5,
          0.3, 0.5, 1.5, 3.5, 12, 50.5, 0.8, 1.5, 3.5, 0.3,
          0.5, 1, 1, 1, 1, 1.5, 3.5, 1000.5, 1000.5, 0.5,
          0.8, 3.5, 0.3, 0.8, 1.5, 12, 1000.5, 12, 1000.5, 1.5),
  b_t = c(0.5, 0.3, 0.8, 0.3, 12, 1, 0.3, 12, 0.3, 1,
          0.5, 2.5, 1000.5, 0.5, 1.5, 3.5, 1.5, 12, 0.3, 2.5,
          2.5, 0.5, 2.5, 2.5, 12, 0.8, 0.8, 0.8, 1000.5, 1.5,
          12, 3.5, 2.5, 0.5, 0.3, 0.3, 3.5, 1, 1.5, 0.3),
  a_c = c(12, 0.3, 1.5, 0.3, 0.8, 50.5, 3.5, 0.5, 1.5, 0.8,
          1000.5, 0.3, 0.3, 12, 0.3, 50.5, 1.5, 0.8, 0.5, 0.8,
          0.3, 0.3, 0.3, 0.5, 0.5, 3.5, 1.5, 3.5, 0.5, 1.5,
          12, 0.3, 0.3, 50.5, 1, 0.8, 1000.5, 12, 12, 0.8),
  b_c = c(1, 0.8, 1000.5, 0.8, 1, 3.5, 0.8, 0.8, 0.8, 0.8,
          1, 4, 50.5, 0.8, 0.8, 50.5, 1.5, 0.8, 0.5, 4,
          4, 0.3, 4, 4, 3.5, 0.8, 0.5, 0.8, 1, 3.5,
          50.5, 0.8, 4, 1000.5, 1, 1000.5, 0.8, 12, 0.5, 0.5),
  upper = c(0.78147756325279410379, 0.99814248309371167287,
            1.0, 0.99999930346318561614,
            0.99999999687305614588, 0.12754816221221219706,
            0.63185169360613470981, 0.57394455466253133721,
            0.98222197132984332322, 0.9481908891291601313,
            0.10037577211254767805, 0.82947269524096216324,
            0.99011640073675305789, 0.59433463828456951789,
            0.96336946952969003811, 0.99999999996782670512,
            0.33150884141014086879, 0.13816863642211098422,
            0.86371317408546832994, 0.31860889966307426188,
            0.69800257255950772626, 0.60918576977096304155,
            0.85558188487959095894, 0.77513966687846307973,
            0.46228623390197753798, 0.31785262915599758083,
            0.53637529537897118811, 0.99128280000881698558,
            0.70706261048164388536, 0.3847177119395777698,
            0.068325499164829594361, 0.75125177720507794449,
            0.33295432566138725666, 0.91326733770048470407,
            0.78338650747441059305, 0.99999999999999998181,
            1.6216022603770513545e-24, 0.84579210478621164998,
            0.0036725899305870677651, 0.040901440062194422561)
)

test_that("pbetadiff() is accurate when shape parameters are below 1", {
  r <- stress_pbetadiff
  upper <- pbetadiff(r$q, r$a_t, r$b_t, r$a_c, r$b_c, lower.tail = FALSE)
  lower <- pbetadiff(r$q, r$a_t, r$b_t, r$a_c, r$b_c)
  for (i in seq_len(nrow(r))) {
    expect_lt(abs(upper[i] - r$upper[i]), 1e-9,
              label = paste("upper tail error, stress case", r$case[i]))
    expect_lt(abs(lower[i] - (1 - r$upper[i])), 1e-9,
              label = paste("lower tail error, stress case", r$case[i]))
  }
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
