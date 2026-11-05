test_that("scenario_treatment() adds the parameter of the treatment group", {
  x <- data.frame(theta = c(0, 0.2), pi_c = 0.3)
  expect_equal(scenario_treatment(ref_design("D1"), x)$pi_t, c(0.3, 0.5))
  x <- data.frame(theta = 5, mu_c = -1)
  expect_equal(scenario_treatment(ref_design_cont("C1"), x)$mu_t, 4)
  x <- data.frame(theta1 = 0.1, theta2 = 0.2, pi_c1 = 0.2, pi_c2 = 0.3)
  out <- scenario_treatment(ref_design_bin2("B5"), x)
  expect_equal(c(out$pi_t1, out$pi_t2), c(0.3, 0.5))
})

test_that("scenario_treatment() rejects effects outside the range", {
  expect_error(qdm_scenarios(ref_design_hr("H1"), theta = c(0, 1),
                             median_c = 8),
               "'theta' must be positive for the hazard ratio.", fixed = TRUE)
  # 0.72 + 0.3 is above 1
  expect_error(qdm_scenarios(ref_design_tte("T1"), theta = 0.3),
               paste("'theta' gives milestone survival probabilities of the",
                     "treatment group outside (0, 1)."), fixed = TRUE)
})
