test_that("scenario_columns() lists the columns required by a design", {
  expect_identical(scenario_columns(ref_design("D1")),
                   list(names = c("pi_t", "pi_c"), lower = 0, upper = 1))
  expect_identical(scenario_columns(ref_design("D3"))$names, "pi_t")
  expect_identical(scenario_columns(ref_design("D4"))$names, "pi_t")
})
