test_that("threshold_effects() gives no effect and the thresholds", {
  expect_identical(threshold_effects(crit_post()), list(c(0, 0.1, 0.2)))
  expect_identical(threshold_effects(ref_design("D5")$criterion),
                   list(c(0, 0.05)))
  # Hazard ratio: no effect is 1, equal to the minimum acceptable value
  expect_identical(threshold_effects(crit_hr(), "hazard_ratio"),
                   list(c(1, 0.7)))
  expect_identical(threshold_effects(crit_tte(), "milestone"),
                   list(c(0, 0.15)))
  expect_identical(threshold_effects(crit_cont2()),
                   list(c(0, 5, 10), c(0, 1, 2)))
  expect_identical(threshold_effects(ref_design_bin2("B1")$criterion),
                   list(c(0, 0.1), 0))
})
