test_that("decide_continuous() returns the posterior parameters", {
  # Python_reference/ref_continuous_qdm.py, design C2
  res <- decide_continuous(ref_design_cont("C2"), ref_data_cont$C2$t,
                           ref_data_cont$C2$c)
  expect_equal(res$posterior_t,
               c(mu = 19.3333333333, kappa = 15, nu = 17,
                 sigma = 23.6610045699), tolerance = 1e-10)
  expect_equal(res$posterior_c,
               c(mu = -10.1428571429, kappa = 21, nu = 23,
                 sigma = 19.6574705664), tolerance = 1e-10)
  expect_equal(res$prob_go, 0.472868606954, tolerance = 1e-9)
  expect_equal(res$prob_nogo, 0.0337620028253, tolerance = 1e-9)
})

test_that("decide_continuous() validates the data", {
  des <- ref_design_cont("C1")
  expect_error(decide_continuous(des, data_binary(12, 3), NULL),
               "'data_t' must be created by data_continuous\\(\\).")
  expect_error(decide_continuous(des, data_continuous(12, 1, 1), NULL),
               "'data_c' must be created by data_continuous\\(\\) for a")
  expect_error(decide_continuous(ref_design_cont("C3"),
                                 data_continuous(12, 1, 1),
                                 data_continuous(12, 1, 1)),
               "'data_c' must be NULL for a design without a concurrent")
})

test_that("decide_continuous() checks the degrees of freedom for MM", {
  # The observed sample sizes may be smaller than those of the design
  des <- ref_design_cont("C1", method = "MM")
  expect_error(decide_continuous(des, data_continuous(5, 18, 25),
                                 data_continuous(12, -12, 25)),
               "\\(here 4 and 11\\)")
})

test_that("decide_continuous() handles a fixed control", {
  # Reference: Python_reference/ref_fixed_tte_qdm.py, design C7
  res <- decide_continuous(ref_design_cont("C7"), ref_data_cont$C7$t, NULL)
  expect_equal(res$posterior_t, c(mu = 25, kappa = 12, nu = 11, sigma = 20))
  expect_identical(res$posterior_c, c(value = -10))
  expect_equal(res$prob_go, 0.797515237746, tolerance = 1e-10)
  expect_equal(res$prob_nogo, 0.00264736615439, tolerance = 1e-10)
  # No check of the degrees of freedom for MM, which is not used
  res_mm <- decide_continuous(ref_design_cont("C7", "MM"),
                              data_continuous(3, 25, 20), NULL)
  expect_equal(res_mm$prob_go, stats::pt(-5 / (20 / sqrt(3)), 2,
                                         lower.tail = FALSE))
})
