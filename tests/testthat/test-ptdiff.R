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

# Stress cases: Python_reference/ref_ptdiff_stress_mpmath.py (mpmath, 30
# digits, integrals over the control and over the treatment variable, which
# agree). Locations mu_t = 0.3 and mu_c = -0.2; scale ratios from 1e-9 to
# 1e3 and degrees of freedom from 0.5 to Inf. The integration scheme used
# before the change of variables had absolute errors above 1e-9 in 15 of
# these cases (Python_reference/compare_ptdiff_stress.py).
stress_ptdiff <- data.frame(
  case = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
           10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
           20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
           30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
           40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
           50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
           60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
           70, 71),
  q = c(-3, 0, 0.7, 5, 40, 10000, -3, 0, 0, 0.7,
        5, 40, 10000, -3, 0, 0.7, 0.7, 5, 40, 10000,
        -3, 0, 0.7, 5, 5, 40, 10000, -3, 0, 0.7,
        5, 40, 40, 10000, -3, 0, 0.7, 5, 40, 10000,
        10000, -3, 0, 0.7, 5, 40, 10000, -3, -3, 0,
        0.7, 5, 40, 10000, -3, 0, 0, 0.7, 5, 40,
        10000, -3, 0, 0.7, 0.7, 5, 40, 10000, -3, 0,
        0.7, 5),
  scale_t = c(1, 1, 1, 1, 1, 1, 1, 1, 0.001, 0.001,
              0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 1, 1, 1, 1,
              1, 1, 1, 1, 50, 50, 50, 50, 50, 50,
              50, 50, 1, 1, 1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1, 1, 1, 1e-06, 1e-06,
              1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 2, 2, 2, 2,
              2, 2, 2, 2, 1, 1, 1, 1, 1, 1,
              1, 1),
  df_t = c(1, 2.5, Inf, 0.5, 3, 100, Inf, 0.5, 1, 2.5,
           Inf, 0.5, 3, 100, Inf, 0.5, 1, 2.5, Inf, 0.5,
           3, 100, Inf, 0.5, 1, 2.5, Inf, 0.5, 3, 100,
           Inf, 0.5, 1, 2.5, Inf, 0.5, 3, 100, Inf, 0.5,
           1, 2.5, Inf, 0.5, 3, 100, Inf, 0.5, 1, 2.5,
           Inf, 0.5, 3, 100, Inf, 0.5, 1, 2.5, Inf, 0.5,
           3, 100, Inf, 0.5, 1, 2.5, Inf, 0.5, 3, 100,
           Inf, 0.5),
  scale_c = c(0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 1, 1,
              1, 1, 1, 1, 1, 1, 1000, 1000, 1000, 1000,
              1000, 1000, 1000, 1000, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
              0.02, 0.02, 1, 1, 1, 1, 1, 1, 1, 1,
              1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1, 1,
              1, 1, 1, 1, 1, 1, 3, 3, 3, 3,
              3, 3, 3, 3, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09,
              1e-09, 1e-09),
  df_c = c(1, 30, 3, 0.5, Inf, 1, Inf, 50, 1, 30,
           3, 0.5, Inf, 1, Inf, 50, 1, 30, 3, 0.5,
           Inf, 1, Inf, 50, 1, 30, 3, 0.5, Inf, 1,
           Inf, 50, 1, 30, 3, 0.5, Inf, 1, Inf, 50,
           1, 30, 3, 0.5, Inf, 1, Inf, 50, 1, 30,
           3, 0.5, Inf, 1, Inf, 50, 1, 30, 3, 0.5,
           Inf, 1, Inf, 50, 1, 30, 3, 0.5, Inf, 1,
           Inf, 50),
  upper = c(0.91133039170764425418, 0.67115093567010934752,
            0.42074040765924459801, 0.15394679082069016453,
            0.000017850414896238118547, 3.1832580572246072843e-8,
            0.99976736939376622699, 0.62134085990839123404,
            0.64745639548047548244, 0.42141527088928027295,
            0.010245208541066846977, 0.0524735833335331264,
            1.1028232722294309817e-21, 0.91141446074178731388,
            0.69146237325774192374, 0.422790303642966117,
            0.49993640162198838471, 0.49821966580044579444,
            0.48548671355610357016, 0.10389508478043332841,
            0.50139629304005562656, 0.50015915476742657095,
            0.49992021158434583382, 0.49825110004024779717,
            0.4714404483635505827, 0.24877584569437109426,
            8.8239094167393032372e-18, 0.51854475341116702796,
            0.50367544390155720288, 0.4984087287039126337,
            0.46414361027560233783, 0.32883148281684476407,
            0.016103204440453942834, 7.194296606551776707e-11,
            0.9709385519672941852, 0.55020871728690150942,
            0.44944380850956707758, 0.073394727911349497345,
            1.2196097806026371956e-171, 0.0032070899425161695899,
            0.000031832611973852153915, 0.9738272265198326242,
            0.69146246127374905511, 0.4471399630492135827,
            0.010245206172229073533, 8.063758991918860477e-9,
            8.3925074411958435885e-21712558, 0.82960734002820512989,
            0.91141438313524486913, 0.68963849755699125182,
            0.42713516462318657448, 0.15073157555078742611,
            1.7960747258925735792e-23, 0.000031832580141272224882,
            0.99976737092096294777, 0.69024244472609524235,
            0.53172551743055356951, 0.48037994515799900416,
            0.14540739116563212103, 0.14767836587184530183,
            8.8225891425226291726e-12, 0.74572201641115103476,
            0.5551465322334494115, 0.48676599817665109517,
            0.43716704187221225578, 0.014947029504366996746,
            1.7960814320930068024e-32, 0.0032071906965681227405,
            0.98025948119035860002, 0.6909132169382523281,
            0.420740290560896977, 0.15062658127776870293)
)

test_that("ptdiff() with NI is accurate for extreme scales and tails", {
  r <- stress_ptdiff
  upper <- ptdiff(r$q, 0.3, r$scale_t, r$df_t, -0.2, r$scale_c, r$df_c,
                  lower.tail = FALSE)
  lower <- ptdiff(r$q, 0.3, r$scale_t, r$df_t, -0.2, r$scale_c, r$df_c)
  for (i in seq_len(nrow(r))) {
    expect_lt(abs(upper[i] - r$upper[i]), 1e-9,
              label = paste("upper tail error, stress case", r$case[i]))
    expect_lt(abs(lower[i] - (1 - r$upper[i])), 1e-9,
              label = paste("lower tail error, stress case", r$case[i]))
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
               "'nMC' must be at least 1 and at most 2147483647.")
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
