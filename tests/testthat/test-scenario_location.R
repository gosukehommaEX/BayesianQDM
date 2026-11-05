test_that("scenario_location() names the columns that set the effect", {
  expect_identical(scenario_location(ref_design("D1")), "pi_t")
  expect_identical(scenario_location(ref_design_cont("C1")), "mu_t")
  expect_identical(scenario_location(ref_design_tte("T2")), "median_t")
  expect_identical(scenario_location(ref_design_bin2("B5")),
                   c("pi_t1", "pi_t2"))
  expect_identical(scenario_location(ref_design_cont2("K1")),
                   c("mu_t1", "mu_t2"))
})
