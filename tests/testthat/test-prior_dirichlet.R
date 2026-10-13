test_that("prior_dirichlet() stores a named parameter vector", {
  p <- prior_dirichlet(c(0.25, 0.5, 1, 2))
  expect_s3_class(p, c("qdm_prior_dirichlet", "qdm_prior"), exact = TRUE)
  expect_identical(p$alpha, c("00" = 0.25, "01" = 0.5, "10" = 1, "11" = 2))
})

test_that("prior_dirichlet() rejects invalid parameters", {
  expect_error(prior_dirichlet(c(1, 1, 1)), "'alpha' must have length 4.")
  expect_error(prior_dirichlet(c(1, 1, 1, 0)),
               "'alpha' must be greater than 0.")
})
