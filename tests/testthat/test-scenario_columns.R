test_that("scenario_columns() lists the columns required by a binary design", {
  expect_identical(scenario_columns(ref_design("D1")),
                   list(names = c("pi_t", "pi_c"), lower = c(0, 0),
                        upper = c(1, 1), lower_open = c(FALSE, FALSE)))
  expect_identical(scenario_columns(ref_design("D3")),
                   list(names = "pi_t", lower = 0, upper = 1,
                        lower_open = FALSE))
  expect_identical(scenario_columns(ref_design("D4"))$names, "pi_t")
})

test_that("scenario_columns() lists the columns required by a continuous design", {
  expect_identical(scenario_columns(ref_design_cont("C1")),
                   list(names = c("mu_t", "sd_t", "mu_c", "sd_c"),
                        lower = c(-Inf, 0, -Inf, 0), upper = rep(Inf, 4),
                        lower_open = c(FALSE, TRUE, FALSE, TRUE)))
  expect_identical(scenario_columns(ref_design_cont("C3")),
                   list(names = c("mu_t", "sd_t"), lower = c(-Inf, 0),
                        upper = c(Inf, Inf), lower_open = c(FALSE, TRUE)))
  expect_identical(scenario_columns(ref_design_cont("C4"))$names,
                   c("mu_t", "sd_t"))
})
