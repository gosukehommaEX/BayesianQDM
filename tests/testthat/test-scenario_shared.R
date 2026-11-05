test_that("scenario_shared() maps the names without the group to columns", {
  expect_identical(scenario_shared("continuous",
                                   c("mu_t", "sd_t", "mu_c", "sd_c")),
                   list(sd = c("sd_t", "sd_c")))
  # Without a concurrent control group, only the treatment group
  expect_identical(scenario_shared("continuous", c("mu_t", "sd_t")),
                   list(sd = "sd_t"))
  expect_identical(scenario_shared("binary2", c("pi_t1", "pi_t2", "rho_t")),
                   list(rho = "rho_t"))
  cols <- scenario_columns(ref_design_cont2("K1"))$names
  expect_identical(scenario_shared("continuous2", cols),
                   list(sd1 = c("sd_t1", "sd_c1"), sd2 = c("sd_t2", "sd_c2"),
                        rho = c("rho_t", "rho_c")))
  expect_identical(scenario_shared("binary", c("pi_t", "pi_c")), list())
  expect_identical(scenario_shared("tte", c("median_t", "median_c")), list())
})
