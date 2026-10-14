test_that("criterion_posterior() for one endpoint ignores the regions", {
  cr <- criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
                            gamma_go = 0.5, gamma_nogo = 0.8)
  expect_s3_class(cr, c("qdm_criterion_posterior", "qdm_criterion"),
                  exact = TRUE)
  expect_identical(cr$n_endpoints, 1L)
  expect_identical(cr$theta_tv, 0.2)
  expect_identical(cr$gamma_nogo, 0.8)
  expect_null(cr$go_regions)
  expect_null(cr$nogo_regions)
})

test_that("criterion_posterior() for two endpoints stores the regions", {
  cr <- criterion_posterior(theta_tv = c(0.2, 0.15), theta_mav = c(0.1, 0),
                            gamma_go = 0.3, gamma_nogo = 0.4,
                            go_regions = c(4, 1), nogo_regions = c(9, 8))
  expect_identical(cr$n_endpoints, 2L)
  expect_identical(cr$go_regions, c(1L, 4L))
  expect_identical(cr$nogo_regions, c(8L, 9L))
  default <- criterion_posterior(c(0.2, 0.15), c(0.1, 0), 0.3, 0.4)
  expect_identical(default$go_regions, 1L)
  expect_identical(default$nogo_regions, 9L)
})

test_that("criterion_posterior() rejects invalid input", {
  expect_error(criterion_posterior(c(1, 2, 3), c(0, 0, 0), 0.5, 0.5),
               "'theta_tv' must have length 1")
  expect_error(criterion_posterior(0.2, c(0.1, 0), 0.5, 0.5),
               "'theta_mav' must have length 1.")
  expect_error(criterion_posterior(0.1, 0.1, 0.5, 0.5),
               "must be greater than that of 'theta_mav'")
  expect_error(criterion_posterior(0.2, 0.1, 1, 0.5),
               "'gamma_go' must be greater than 0 and less than 1.")
  expect_error(criterion_posterior(0.2, 0.1, 0.5, 0),
               "'gamma_nogo' must be greater than 0 and less than 1.")
  expect_error(criterion_posterior(c(0.2, 0.2), c(0, 0), 0.5, 0.5,
                                   go_regions = 1, nogo_regions = 1),
               "must not share any region")
})

test_that("criterion_posterior() coerces region numbers to integers", {
  cr <- criterion_posterior(c(0.2, 0.15), c(0.1, 0), 0.3, 0.4,
                            go_regions = c(1, 2, 4), nogo_regions = 9)
  expect_type(cr$go_regions, "integer")
  expect_identical(cr$go_regions, c(1L, 2L, 4L))
  expect_error(criterion_posterior(c(0.2, 0.15), c(0.1, 0), 0.3, 0.4,
                                   go_regions = 0),
               "'go_regions' must be at least 1 and at most 9.")
})

test_that("criterion_posterior() supports marginal rules for two endpoints", {
  cr <- criterion_posterior(c(0.2, 0.15), c(0.1, 0), gamma_go = c(0.6, 0.5),
                            gamma_nogo = 0.7, rule = "and")
  expect_identical(cr$rule, "and")
  expect_identical(cr$gamma_go, c(0.6, 0.5))
  expect_identical(cr$gamma_nogo, c(0.7, 0.7))
  expect_null(cr$go_regions)
  expect_identical(criterion_posterior(0.2, 0.1, 0.5, 0.5, rule = "or")$rule,
                   "joint")
  expect_identical(criterion_posterior(c(0.2, 0.1), c(0, 0), 0.5,
                                       0.5)$rule, "joint")
  expect_error(criterion_posterior(c(0.2, 0.15), c(0.1, 0), 0.5, 0.5,
                                   go_regions = 1:2, rule = "or"),
               "used only with rule = \"joint\"")
  expect_error(criterion_posterior(c(0.2, 0.15), c(0.1, 0), c(0.5, 0.5),
                                   0.5),
               "'gamma_go' must have length 1.")
  expect_error(criterion_posterior(c(0.2, 0.15), c(0.1, 0), c(0.5, 0.5, 0.5),
                                   0.5, rule = "and"),
               "must have length 1 or 2 \\(one threshold per endpoint\\)")
  expect_error(criterion_posterior(c(0.2, 0.15), c(0.1, 0), c(0.5, 1),
                                   0.5, rule = "and"),
               "'gamma_go' must be greater than 0 and less than 1.")
  expect_error(criterion_posterior(0.2, 0.1, 0.5, 0.5, rule = "xor"),
               "'arg' should be one of")
})
