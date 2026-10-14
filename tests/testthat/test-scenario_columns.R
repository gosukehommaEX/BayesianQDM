test_that("scenario_columns() lists the columns required by a binary design", {
  expect_identical(scenario_columns(ref_design("D1")),
                   list(names = c("pi_t", "pi_c"), lower = c(0, 0),
                        upper = c(1, 1), lower_open = c(FALSE, FALSE),
                        upper_open = c(FALSE, FALSE)))
  expect_identical(scenario_columns(ref_design("D3")),
                   list(names = "pi_t", lower = 0, upper = 1,
                        lower_open = FALSE, upper_open = FALSE))
  expect_identical(scenario_columns(ref_design("D4"))$names, "pi_t")
})

test_that("scenario_columns() lists the columns required by a continuous design", {
  expect_identical(scenario_columns(ref_design_cont("C1")),
                   list(names = c("mu_t", "sd_t", "mu_c", "sd_c"),
                        lower = c(-Inf, 0, -Inf, 0), upper = rep(Inf, 4),
                        lower_open = c(FALSE, TRUE, FALSE, TRUE),
                        upper_open = rep(FALSE, 4)))
  expect_identical(scenario_columns(ref_design_cont("C3")),
                   list(names = c("mu_t", "sd_t"), lower = c(-Inf, 0),
                        upper = c(Inf, Inf), lower_open = c(FALSE, TRUE),
                        upper_open = c(FALSE, FALSE)))
  expect_identical(scenario_columns(ref_design_cont("C4"))$names,
                   c("mu_t", "sd_t"))
})

test_that("scenario_columns() lists the columns for two binary endpoints", {
  cols <- scenario_columns(ref_design_bin2("B1"))
  expect_identical(cols$names, c("pi_t1", "pi_t2", "rho_t", "pi_c1",
                                 "pi_c2", "rho_c"))
  expect_identical(cols$lower, rep(c(0, 0, -1), 2))
  expect_identical(cols$upper_open, rep(c(TRUE, TRUE, FALSE), 2))
  expect_identical(scenario_columns(ref_design_bin2("B7"))$names,
                   c("pi_t1", "pi_t2", "rho_t"))
})
