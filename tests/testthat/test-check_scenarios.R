test_that("check_scenarios() returns a plain data frame", {
  des <- ref_design("D1")
  sc <- check_scenarios(des, scen_c)
  expect_identical(sc, scen_c)
  # A tibble-like data frame is converted
  tb <- structure(scen_c, class = c("tbl_df", "tbl", "data.frame"))
  expect_identical(class(check_scenarios(des, tb)), "data.frame")
})

test_that("check_scenarios() validates the scenarios", {
  des <- ref_design("D1")
  expect_error(check_scenarios(des, scen_c[0, ]),
               "'scenarios' must be a data frame with at least one row.")
  expect_error(check_scenarios(des, scen_c, "sc", one_row = TRUE),
               "'sc' must have exactly one row.")
  expect_error(check_scenarios(des, data.frame(pi_c = 0.3)),
               "must have the column\\(s\\) pi_t, pi_c.")
  expect_error(check_scenarios(ref_design("D3"), data.frame(x = 0.3)),
               "must have the column\\(s\\) pi_t.")
  expect_error(check_scenarios(des, data.frame(pi_t = 0.3,
                                               pi_c = NA_real_)),
               "'scenarios\\$pi_c' must not contain missing values.")
  expect_error(check_scenarios(des, data.frame(pi_t = -0.1, pi_c = 0.2)),
               "'scenarios\\$pi_t' must be at least 0 and at most 1.")
  expect_error(check_scenarios(des, data.frame(pi_t = "a", pi_c = 0.2)),
               "'scenarios\\$pi_t' must be numeric.")
  expect_error(check_scenarios(des, data.frame(pi_t = 0.3, pi_c = 0.2,
                                               Miss = 0)),
               "must not have columns named Go, NoGo, Gray or Miss.")
})

test_that("check_scenarios() validates the scenarios of continuous designs", {
  des <- ref_design_cont("C1")
  expect_identical(check_scenarios(des, scen_cont_c), scen_cont_c)
  expect_error(check_scenarios(des, data.frame(mu_t = 1, sd_t = 1, mu_c = 0)),
               "must have the column\\(s\\) mu_t, sd_t, mu_c, sd_c.")
  expect_error(check_scenarios(des, data.frame(mu_t = 1, sd_t = 1, mu_c = 0,
                                               sd_c = -1)),
               "'scenarios\\$sd_c' must be greater than 0.")
  expect_error(check_scenarios(des, data.frame(mu_t = Inf, sd_t = 1,
                                               mu_c = 0, sd_c = 1)),
               "'scenarios\\$mu_t' must be finite.")
  u <- data.frame(mu_t = c(-5, 5), sd_t = c(1, 2))
  expect_identical(check_scenarios(ref_design_cont("C3"), u), u)
})
