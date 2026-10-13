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
