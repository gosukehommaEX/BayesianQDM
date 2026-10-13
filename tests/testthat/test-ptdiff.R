# Reference values: Python_reference/ref_distributions_mpmath.py (mpmath, 40
# digits, two integral representations), cross-checked with scipy.
s1 <- sqrt(15.62 / 36)
s2 <- sqrt(13.28 / 31)
ref_ptdiff <- data.frame(
  q       = c(1, 0, 1, 2, 0.5, 0, -2, 1),
  mu_t    = c(3, 1, 5, 0, 0, -5.2, -5.2, 2),
  scale_t = c(1, 1, 2, 1, 0.01, s1, s1, 1),
  df_t    = c(15, 5, 10, 3, 10, 35, 35, Inf),
  mu_c    = c(1, 1, 3, 0, 0, -2.7, -2.7, 0),
  scale_c = c(1, 1, 1.5, 0.2, 1, s2, s2, 1.5),
  df_c    = c(15, 20, 15, 3, 4.5, 30, 30, Inf),
  ni = c(0.7502661261453839, 0.5, 0.64781149700404485, 0.074147617800499862,
         0.3202865229079607, 0.0051019590255298554, 0.2989120884753876,
         0.71045012902304058),
  mm = c(0.75034682914089527, 0.5, 0.64835811602628018, NA,
         0.32028123131924755, 0.0051031869369207273, 0.29890417235967946,
         0.71045012902304058)
)

test_that("ptdiff() with NI matches high-precision references", {
  r <- ref_ptdiff
  for (i in seq_len(nrow(r))) {
    expect_equal(ptdiff(r$q[i], r$mu_t[i], r$scale_t[i], r$df_t[i],
                        r$mu_c[i], r$scale_c[i], r$df_c[i],
                        lower.tail = FALSE),
                 r$ni[i], tolerance = 1e-9, info = paste("case", i))
  }
})

test_that("ptdiff() with MM matches the closed-form references", {
  r <- ref_ptdiff[!is.na(ref_ptdiff$mm), ]
  for (i in seq_len(nrow(r))) {
    expect_equal(ptdiff(r$q[i], r$mu_t[i], r$scale_t[i], r$df_t[i],
                        r$mu_c[i], r$scale_c[i], r$df_c[i], method = "MM",
                        lower.tail = FALSE),
                 r$mm[i], tolerance = 1e-10, info = paste("case", i))
  }
})

test_that("ptdiff() reproduces the glaucoma example of Yamaguchi et al.", {
  # Moment-matching paper, Section 5.2: P(Delta < 0) = 0.9949 and
  # P(Delta <= -2) = 0.7011
  for (m in c("NI", "MM")) {
    p0 <- ptdiff(0, -5.2, s1, 35, -2.7, s2, 30, method = m)
    p2 <- ptdiff(-2, -5.2, s1, 35, -2.7, s2, 30, method = m)
    expect_equal(round(p0, 4), 0.9949, info = m)
    expect_equal(round(p2, 4), 0.7011, info = m)
  }
})

test_that("ptdiff() lower tail is the complement of the upper tail", {
  r <- ref_ptdiff
  lower <- ptdiff(r$q, r$mu_t, r$scale_t, r$df_t, r$mu_c, r$scale_c, r$df_c)
  expect_equal(lower, 1 - r$ni, tolerance = 1e-9)
})

test_that("ptdiff() with MC agrees with NI within Monte Carlo error", {
  r <- ref_ptdiff[c(1, 3, 6), ]
  set.seed(2026)
  est <- ptdiff(r$q, r$mu_t, r$scale_t, r$df_t, r$mu_c, r$scale_c, r$df_c,
                method = "MC", nMC = 1e5, lower.tail = FALSE)
  se <- sqrt(r$ni * (1 - r$ni) / 1e5)
  expect_true(all(abs(est - r$ni) < 4 * se))
})

test_that("ptdiff() with MC is reproducible with a seed", {
  set.seed(1)
  a <- ptdiff(1, 3, 1, 15, 1, 1, 15, method = "MC", nMC = 1000)
  set.seed(1)
  b <- ptdiff(1, 3, 1, 15, 1, 1, 15, method = "MC", nMC = 1000)
  expect_identical(a, b)
})

test_that("ptdiff() validates its arguments", {
  expect_error(ptdiff(0, 0, 1, 3, 0, 1, 10, method = "MM"),
               "requires 'df_t' and 'df_c' greater than 4")
  expect_error(ptdiff(0, 0, 0, 3, 0, 1, 10), "'scale_t' must be greater")
  expect_error(ptdiff(0, 0, 1, 3, 0, 1, 10, method = "MC", nMC = 0),
               "'nMC' must be at least 1.")
  expect_error(ptdiff(0, 0, 1, 3, 0, 1, 10, method = "XX"),
               "should be one of")
})

test_that("ptdiff() is exact for Cauchy variables (one degree of freedom)", {
  # T_t - T_c is Cauchy with location mu_t - mu_c and scale s_t + s_c
  grid <- expand.grid(q = c(-5, -1, 0, 0.3, 2, 20),
                      s = 1:3)
  sc <- list(c(1, 1), c(0.2, 3), c(5, 0.5))
  for (k in seq_len(nrow(grid))) {
    s <- sc[[grid$s[k]]]
    exact <- stats::pcauchy(grid$q[k], 1.5 - 0.5, s[1] + s[2],
                            lower.tail = FALSE)
    expect_equal(ptdiff(grid$q[k], 1.5, s[1], 1, 0.5, s[2], 1,
                        lower.tail = FALSE),
                 exact, tolerance = 1e-8, info = paste("case", k))
  }
})

test_that("ptdiff() is exact for normal variables (infinite df)", {
  q <- c(-3, -0.5, 0, 0.7, 4)
  exact <- stats::pnorm(q, 2 - 0.5, sqrt(1.2 ^ 2 + 0.8 ^ 2),
                        lower.tail = FALSE)
  for (m in c("NI", "MM")) {
    expect_equal(ptdiff(q, 2, 1.2, Inf, 0.5, 0.8, Inf, method = m,
                        lower.tail = FALSE),
                 exact, tolerance = 1e-9, info = m)
  }
  # Large degrees of freedom are close to the normal case
  expect_equal(ptdiff(q, 2, 1.2, 1e7, 0.5, 0.8, 1e7, lower.tail = FALSE),
               exact, tolerance = 1e-6)
})

test_that("ptdiff() with MM is close to NI for moderate degrees of freedom", {
  # Accuracy of the moment-matching approximation (Yamaguchi et al.)
  grid <- expand.grid(q = c(-2, 0, 1, 3), df_t = c(10, 30), df_c = c(12, 60),
                      scale_c = c(0.5, 1, 2))
  ni <- ptdiff(grid$q, 0.5, 1, grid$df_t, 0, grid$scale_c, grid$df_c)
  mm <- ptdiff(grid$q, 0.5, 1, grid$df_t, 0, grid$scale_c, grid$df_c,
               method = "MM")
  expect_lt(max(abs(ni - mm)), 2e-3)
})

test_that("ptdiff() is a non-increasing function of q", {
  q <- seq(-10, 10, length.out = 201)
  p <- ptdiff(q, 1, 1, 3, 0, 2, 7, lower.tail = FALSE)
  expect_true(all(diff(p) <= 1e-10))
  expect_true(all(p >= 0 & p <= 1))
})

test_that("ptdiff() handles long vectors", {
  q <- seq(-3, 3, length.out = 1200)
  p <- ptdiff(q, 1, 1, 8, 0, 1.5, 12)
  expect_length(p, 1200L)
  expect_identical(p[c(1, 600, 1200)],
                   ptdiff(q[c(1, 600, 1200)], 1, 1, 8, 0, 1.5, 12))
})

test_that("ptdiff() with MC uses the default nMC and both tails", {
  set.seed(5)
  up <- ptdiff(1, 3, 1, 15, 1, 1, 15, method = "MC", lower.tail = FALSE)
  lo <- ptdiff(1, 3, 1, 15, 1, 1, 15, method = "MC")
  ref <- 0.7502661261453839
  se <- sqrt(ref * (1 - ref) / 10000)
  expect_lt(abs(up - ref), 4 * se)
  expect_lt(abs(lo - (1 - ref)), 4 * se)
})

test_that("ptdiff() with MC agrees with NI for many draws", {
  skip_on_cran()
  set.seed(6)
  ref <- ptdiff(c(-1, 0.5, 2), 1, 1, 6, 0, 1.5, 9, lower.tail = FALSE)
  est <- ptdiff(c(-1, 0.5, 2), 1, 1, 6, 0, 1.5, 9, method = "MC",
                nMC = 1e6, lower.tail = FALSE)
  expect_true(all(abs(est - ref) < 4 * sqrt(ref * (1 - ref) / 1e6)))
})

test_that("ptdiff() checks the lengths of its arguments", {
  expect_error(ptdiff(c(0, 1), c(0, 1, 2), 1, 5, 0, 1, 5),
               "'q' must have length 1 or 3.")
  expect_error(ptdiff(0, 0, 1, 0, 0, 1, 5), "'df_t' must be greater than 0.")
  expect_error(ptdiff(0, 0, 1, 5, 0, 1, 5, lower.tail = "yes"),
               "'lower.tail' must be TRUE or FALSE.")
})
