test_that("boundary_tte() specifies the grid of events", {
  des <- ref_design_tte("T1")
  s <- boundary_tte(des, 10, NULL, c(150, 100), NULL)
  expect_identical(s$cond, data.frame(exposure_t = c(100, 150)))
  expect_equal(s$values, list(0:10, 0:10))
  expect_true(s$discrete)
  expect_false(s$increasing)
  p <- s$prob(2L, 6)
  expect_equal(c(p$go, p$nogo), unname(ref_prob_tte$T1), tolerance = 1e-10)
  d <- s$derived(2L, 6)
  expect_equal(d$surv, exp(-3 * 6 / 150))
  expect_equal(d$effect, exp(-3 * 6 / 150) - 0.72)
  expect_match(s$titles[["surv"]], "exp\\(-3 d_t / E_t\\)")
})

test_that("boundary_tte() measures the hazard ratio", {
  s <- boundary_tte(ref_design_hr("H1"), 30, 0:30, 200, ref_data_hr$c)
  expect_true(s$increasing == FALSE)
  expect_null(s$derived(1L, 12)$surv)
  expect_equal(s$derived(1L, 12)$effect, (12 / 200) / (20 / 180))
  expect_match(s$titles[["effect"]], "observed in the control group")
  p <- s$prob(1L, 12)
  expect_equal(c(p$go, p$nogo), unname(ref_prob_hr$H1), tolerance = 1e-10)
  expect_match(s$axes$caption, "^Control data: ")
  # Fixed control: relative to the fixed hazard log(2) / 6
  s2 <- boundary_tte(ref_design_hr("H2"), 30, 12, 200, NULL)
  expect_equal(s2$derived(1L, 12)$effect, (12 / 200) / (log(2) / 6))
  expect_match(s2$titles[["effect"]], "fixed control value")
  # Milestone survival of a concurrent control: the observed value
  s3 <- boundary_tte(ref_design_tte("T2"), 20, 5, 40, data_tte(9, 35))
  expect_equal(s3$derived(1L, 5)$effect, exp(-3 * 5 / 40) - exp(-3 * 9 / 35))
  s5 <- boundary_tte(ref_design_tte("T5"), 30, 5, 150, NULL)
  expect_match(s5$titles[["effect"]], "mean of the control distribution")
  # A concurrent control without events: the observed hazard ratio is not
  # defined
  s6 <- boundary_tte(ref_design_hr("H1"), 30, 0:30, 200, data_tte(0, 180))
  expect_identical(s6$derived(1L, c(0, 12))$effect, c(NA_real_, NA_real_))
})

test_that("boundary_tte() validates its arguments", {
  des <- ref_design_tte("T2")
  expect_error(boundary_tte(des, 20, NULL, NULL, data_tte(9, 35)),
               "'exposure_t' must be given for a time-to-event endpoint.")
  expect_error(boundary_tte(des, 20, 1.5, 40, data_tte(9, 35)),
               "'events_t' must contain whole numbers only.")
  expect_error(boundary_tte(des, 20, 5, 0, data_tte(9, 35)),
               "'exposure_t' must be greater than 0.")
  expect_error(boundary_tte(des, 20, 5, 40, data_binary(9, 3)),
               "'data_c' must be created by data_tte\\(\\).")
})
