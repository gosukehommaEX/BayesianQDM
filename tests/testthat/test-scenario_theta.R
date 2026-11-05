test_that("scenario_theta() gives one vector of effects per endpoint", {
  expect_identical(scenario_theta(ref_design("D1"), c(0, 0.1)),
                   list(c(0, 0.1)))
  expect_identical(scenario_theta(ref_design("D1"), 1L), list(1))
  expect_identical(scenario_theta(ref_design_cont2("K1"),
                                  list(a = 0, b = 1:2)),
                   list(0, c(1, 2)))
  expect_identical(scenario_theta(ref_design("D1"), "thresholds"),
                   list(c(0, 0.1, 0.2)))
})

test_that("scenario_theta() rejects other values", {
  expect_error(scenario_theta(ref_design("D1"), "tv"),
               paste("'theta' must be numeric (a list of two numeric vectors",
                     "for two endpoints) or \"thresholds\"."), fixed = TRUE)
  expect_error(scenario_theta(ref_design("D1"), list(0.1)),
               "'theta' must be numeric.", fixed = TRUE)
  expect_error(scenario_theta(ref_design_cont2("K1"), c(0, 1)),
               "For two endpoints, 'theta' must be a list of two numeric",
               fixed = TRUE)
  expect_error(scenario_theta(ref_design_cont2("K1"), list(0, "a")),
               "'theta[[2]]' must be numeric.", fixed = TRUE)
})
