test_that("qdm_scenarios() gives the columns of the design", {
  des <- ref_design("D1")
  expect_identical(
    qdm_scenarios(des, pi_t = c(0.2, 0.4), pi_c = c(0.3, 0.4)),
    data.frame(pi_t = c(0.2, 0.4, 0.2, 0.4), pi_c = c(0.3, 0.3, 0.4, 0.4))
  )
  # The columns are in the order of qdm_oc() whatever the order given
  expect_identical(qdm_scenarios(des, pi_c = 0.4, pi_t = c(0.2, 0.4)),
                   data.frame(pi_t = c(0.2, 0.4), pi_c = 0.4))
  expect_identical(
    qdm_scenarios(des, pi_t = c(0.2, 0.4), pi_c = c(0.3, 0.4), grid = FALSE),
    data.frame(pi_t = c(0.2, 0.4), pi_c = c(0.3, 0.4))
  )
  # Without a concurrent control group, a specification
  expect_identical(qdm_scenarios(ref_design("D3"), pi_t = c(0.2, 0.5)),
                   data.frame(pi_t = c(0.2, 0.5)))
  expect_identical(qdm_scenarios(spec_d1(), pi_t = 0.3, pi_c = 0.3),
                   data.frame(pi_t = 0.3, pi_c = 0.3))
  # The scenarios are accepted by qdm_oc()
  scen <- qdm_scenarios(des, pi_t = c(0.2, 0.4, 0.6, 0.8), pi_c = 0.4)
  expect_identical(qdm_oc(des, scen), qdm_oc(des, scen_c))
})

test_that("qdm_scenarios() gives a value without the group to both groups", {
  expect_identical(
    qdm_scenarios(ref_design_cont("C1"), mu_t = c(0, 20), mu_c = 0,
                  sd = c(20, 25)),
    data.frame(mu_t = c(0, 20, 0, 20), sd_t = c(20, 20, 25, 25), mu_c = 0,
               sd_c = c(20, 20, 25, 25))
  )
  expect_identical(qdm_scenarios(ref_design_cont("C3"), mu_t = 10, sd = 20),
                   data.frame(mu_t = 10, sd_t = 20))
  s2 <- qdm_scenarios(ref_design_cont2("K1"), mu_t1 = 10, mu_t2 = 2,
                      mu_c1 = 0, mu_c2 = 0, sd1 = 70, sd2 = 6,
                      rho = c(0, 0.5))
  expect_identical(names(s2), scenario_columns(ref_design_cont2("K1"))$names)
  expect_identical(s2$rho_t, c(0, 0.5))
  expect_identical(s2$rho_c, c(0, 0.5))
  expect_identical(s2$sd_c2, c(6, 6))
  s3 <- qdm_scenarios(ref_design_bin2("B5"), pi_t1 = 0.3, pi_t2 = 0.4,
                      pi_c1 = 0.2, pi_c2 = 0.2, rho = 0.1)
  expect_identical(c(s3$rho_t, s3$rho_c), c(0.1, 0.1))
  expect_identical(names(qdm_scenarios(ref_design_bin2("B7"), pi_t1 = 0.3,
                                       pi_t2 = 0.4, rho = 0)),
                   c("pi_t1", "pi_t2", "rho_t"))
})

test_that("qdm_scenarios() computes the treatment group from theta", {
  th <- c(-0.1, 0, 0.15, 0.3)
  s <- qdm_scenarios(ref_design("D1"), theta = th, pi_c = 0.4)
  expect_identical(names(s), c("theta", "pi_t", "pi_c"))
  expect_identical(s$theta, th)
  expect_equal(s$pi_t, 0.4 + th)
  expect_equal(oc_axis(ref_design("D1"), s)$x, th)
  # Without a concurrent control: the mean of Beta(5, 9) (hypothetical
  # control), of Beta(4, 8) (external control data) or the fixed value
  expect_equal(qdm_scenarios(ref_design("D3"), theta = th)$pi_t, 5 / 14 + th)
  expect_equal(qdm_scenarios(ref_design("D4"), theta = th)$pi_t, 1 / 3 + th)
  expect_equal(qdm_scenarios(ref_design("D8"), theta = th)$pi_t, 0.2 + th)
  # The thresholds of the plots are shifted by the same location
  des <- ref_design("D3")
  ax <- oc_axis(des, qdm_scenarios(des, theta = th))
  expect_equal(ax$x - (ax$lines[1L] - des$criterion$theta_mav), th)
  # Continuous endpoint: concurrent, hypothetical (mean -10), fixed (-10)
  # and external control data only (posterior mean -12)
  s <- qdm_scenarios(ref_design_cont("C1"), theta = c(0, 15, 30), mu_c = 5,
                     sd = 25)
  expect_equal(s$mu_t, 5 + c(0, 15, 30))
  expect_identical(names(s), c("theta", "mu_t", "sd_t", "mu_c", "sd_c"))
  expect_equal(qdm_scenarios(ref_design_cont("C3"), theta = 15, sd = 20)$mu_t,
               5)
  expect_equal(qdm_scenarios(ref_design_cont("C7"), theta = 15, sd = 20)$mu_t,
               5)
  expect_equal(qdm_scenarios(ref_design_cont("C4"), theta = 15, sd = 20)$mu_t,
               3)
})

test_that("qdm_scenarios() gives the median survival times for theta", {
  # References: Python_reference/ref_scenarios.out
  s <- qdm_scenarios(ref_design_tte("T2"), theta = c(0, 0.15), median_c = 6)
  expect_equal(s$median_t, c(6, 13.4859860856831), tolerance = 1e-12)
  expect_equal(oc_axis(ref_design_tte("T2"), s)$x, c(0, 0.15))
  expect_equal(qdm_scenarios(ref_design_tte("T1"), theta = c(0, 0.15))$median_t,
               c(6.33003286944648, 14.9318589153208), tolerance = 1e-12)
  expect_equal(qdm_scenarios(ref_design_tte("T5"), theta = c(0, 0.15))$median_t,
               c(5.3473308322019, 11.0046366165297), tolerance = 1e-12)
  s <- qdm_scenarios(ref_design_hr("H1"), theta = c(1, 0.7), median_c = 8)
  expect_equal(s$median_t, c(8, 11.4285714285714), tolerance = 1e-12)
  expect_equal(oc_axis(ref_design_hr("H1"), s)$x, c(1, 0.7))
  expect_equal(qdm_scenarios(ref_design_hr("H2"), theta = c(1, 0.7))$median_t,
               c(6, 8.57142857142857), tolerance = 1e-12)
  expect_equal(qdm_scenarios(ref_design_hr("H3"), theta = c(1, 0.7))$median_t,
               c(4.63628246599699, 6.62326066570998), tolerance = 1e-12)
  expect_equal(oc_axis(ref_design_hr("H3"),
                       qdm_scenarios(ref_design_hr("H3"),
                                     theta = c(1, 0.7)))$x, c(1, 0.7))
})

test_that("qdm_scenarios() gives the effects of two endpoints", {
  s <- qdm_scenarios(ref_design_cont2("K1"), theta = list(c(0, 5), c(0, 1)),
                     mu_c1 = -10, mu_c2 = -1, sd1 = 70, sd2 = 6, rho = 0.4)
  expect_identical(names(s), c("theta1", "theta2",
                               scenario_columns(ref_design_cont2("K1"))$names))
  expect_identical(s$theta1, c(0, 5, 0, 5))
  expect_identical(s$theta2, c(0, 0, 1, 1))
  expect_equal(s$mu_t1 - s$mu_c1, s$theta1)
  expect_equal(s$mu_t2 - s$mu_c2, s$theta2)
  # Without a concurrent control: hypothetical mean, external data and fixed
  # values (control_location())
  s <- qdm_scenarios(ref_design_cont2("K4"), theta = list(5, 1), sd1 = 30,
                     sd2 = 3, rho = 0.5)
  expect_equal(c(s$mu_t1, s$mu_t2), c(-5, 0))
  s <- qdm_scenarios(ref_design_bin2("B8"), theta = list(0.1, 0.2), rho = 0)
  expect_equal(c(s$pi_t1, s$pi_t2), c(2.5 / 11 + 0.1, 2 / 11 + 0.2))
  s <- qdm_scenarios(ref_design_bin2("BF1"), theta = "thresholds", rho = 0)
  expect_identical(s$theta1, c(0, 0.1))
  expect_identical(s$theta2, c(0, 0))
  expect_equal(s$pi_t1, c(0.25, 0.35))
  expect_equal(s$pi_t2, c(0.5, 0.5))
})

test_that("qdm_scenarios() gives the thresholds of the criterion", {
  s <- qdm_scenarios(ref_design("D1"), theta = "thresholds", pi_c = 0.3)
  expect_identical(s$theta, c(0, 0.1, 0.2))
  expect_equal(s$pi_t, c(0.3, 0.4, 0.5))
  expect_identical(qdm_scenarios(ref_design("D5"), theta = "thresholds",
                                 pi_c = 0.3)$theta, c(0, 0.05))
  expect_identical(qdm_scenarios(ref_design_hr("H1"), theta = "thresholds",
                                 median_c = 8)$theta, c(1, 0.7))
  s <- qdm_scenarios(ref_design_cont2("K1"), theta = "thresholds",
                     mu_c1 = 0, mu_c2 = 0, sd1 = 70, sd2 = 6, rho = 0)
  expect_identical(nrow(s), 9L)
  expect_identical(unique(s$theta1), c(0, 5, 10))
  expect_identical(unique(s$theta2), c(0, 1, 2))
})

test_that("qdm_scenarios() rejects missing, unknown and repeated columns", {
  des <- ref_design("D1")
  expect_error(qdm_scenarios(des, pi_t = 0.3),
               "The column(s) 'pi_c' must be given in '...'.", fixed = TRUE)
  expect_error(qdm_scenarios(des, pi_t = 0.3, pi_c = 0.3, pi = 1),
               paste("The scenarios of this design have no column 'pi'; they",
                     "have the column(s) pi_t, pi_c."), fixed = TRUE)
  expect_error(qdm_scenarios(ref_design("D3"), pi_t = 0.3, pi_c = 0.3),
               "no column 'pi_c'", fixed = TRUE)
  expect_error(qdm_scenarios(ref_design_cont("C1"), mu_t = 0, mu_c = 0,
                             sd = 20, sd_t = 20),
               "The column(s) 'sd_t' are given more than once in '...'.",
               fixed = TRUE)
  expect_error(qdm_scenarios(des, 0.3, pi_c = 0.3),
               "The values of the columns in '...' must be named.",
               fixed = TRUE)
  expect_error(qdm_scenarios(des, pi_t = 0.3, pi_t = 0.4, pi_c = 0.3),
               "Each name in '...' must be given only once.", fixed = TRUE)
  expect_error(qdm_scenarios(des, theta = 0.1, pi_t = 0.3, pi_c = 0.3),
               "'pi_t' is computed from 'theta' and must not be given in '...'.",
               fixed = TRUE)
  expect_error(qdm_scenarios(des, pi_t = "a", pi_c = 0.3),
               "'pi_t' must be numeric.", fixed = TRUE)
  expect_error(qdm_scenarios(des, pi_t = 1.2, pi_c = 0.3),
               "'scenarios$pi_t' must be at least 0 and at most 1.",
               fixed = TRUE)
  expect_error(qdm_scenarios(des, pi_t = c(0.2, 0.3, 0.4), pi_c = c(0.3, 0.4),
                             grid = FALSE),
               "With 'grid = FALSE', each value must have length 1",
               fixed = TRUE)
  expect_error(qdm_scenarios(des, pi_t = 0.3, pi_c = 0.3, grid = NA),
               "'grid' must be TRUE or FALSE.", fixed = TRUE)
  expect_error(qdm_scenarios(ref_design_cont2("K1"), theta = list(0, 1),
                             mu_c1 = 0, mu_c2 = 0, sd1 = 1, sd2 = 1),
               "The column(s) 'rho_t', 'rho_c' must be given in '...'.",
               fixed = TRUE)
  expect_error(qdm_scenarios(list(), pi_t = 0.3),
               "'design' must be created by a design_*() function",
               fixed = TRUE)
})
