test_that("criterion_predictive() for one endpoint stores the settings", {
  cr <- criterion_predictive(theta_null = 0, m_t = 50, gamma_go = 0.9,
                             gamma_nogo = 0.7)
  expect_s3_class(cr, c("qdm_criterion_predictive", "qdm_criterion"),
                  exact = TRUE)
  expect_identical(cr$m_t, 50L)
  expect_identical(cr$m_c, 50L)
  expect_null(cr$go_regions)
})

test_that("criterion_predictive() for two endpoints stores the regions", {
  cr <- criterion_predictive(theta_null = c(0, 0), m_t = 40, m_c = 20,
                             gamma_go = 0.5, gamma_nogo = 0.5)
  expect_identical(cr$n_endpoints, 2L)
  expect_identical(cr$m_c, 20L)
  expect_identical(cr$go_regions, 1L)
  expect_identical(cr$nogo_regions, 4L)
})

test_that("criterion_predictive() rejects invalid input", {
  expect_error(criterion_predictive(0, m_t = 0, gamma_go = 0.5,
                                    gamma_nogo = 0.5),
               "'m_t' must be at least 1.")
  expect_error(criterion_predictive(0, m_t = 10, m_c = 2.5, gamma_go = 0.5,
                                    gamma_nogo = 0.5), "whole numbers")
  expect_error(criterion_predictive(c(0, 0), m_t = 10, gamma_go = 0.5,
                                    gamma_nogo = 0.5, go_regions = 5),
               "'go_regions' must be at least 1 and at most 4.")
})

test_that("criterion_predictive() validates the null values and thresholds", {
  expect_error(criterion_predictive(c(0, 0, 0), m_t = 10, gamma_go = 0.5,
                                    gamma_nogo = 0.5),
               "'theta_null' must have length 1")
  expect_error(criterion_predictive(NA_real_, m_t = 10, gamma_go = 0.5,
                                    gamma_nogo = 0.5),
               "'theta_null' must not contain missing")
  expect_error(criterion_predictive(0, m_t = 10, gamma_go = 0,
                                    gamma_nogo = 0.5),
               "'gamma_go' must be greater than 0 and less than 1.")
  expect_error(criterion_predictive(0, m_t = 10, gamma_go = 0.5,
                                    gamma_nogo = 1),
               "'gamma_nogo' must be greater than 0 and less than 1.")
  expect_error(criterion_predictive(c(0, 0), m_t = 10, gamma_go = 0.5,
                                    gamma_nogo = 0.5, go_regions = 1,
                                    nogo_regions = c(1, 4)),
               "must not share any region")
})
