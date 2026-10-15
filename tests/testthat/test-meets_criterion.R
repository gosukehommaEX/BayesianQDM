test_that("meets_criterion() applies the joint and marginal rules", {
  expect_identical(meets_criterion(c(0.2, 0.6), 0.5), c(FALSE, TRUE))
  p <- list(c(0.7, 0.7, 0.2), c(0.6, 0.3, 0.6))
  expect_identical(meets_criterion(p, c(0.5, 0.5), "and"),
                   c(TRUE, FALSE, FALSE))
  expect_identical(meets_criterion(p, c(0.5, 0.5), "or"),
                   c(TRUE, TRUE, TRUE))
  expect_identical(meets_criterion(p, c(0.8, 0.5), "or"),
                   c(TRUE, FALSE, TRUE))
  # A single decision: one probability per endpoint
  expect_true(meets_criterion(c(0.7, 0.6), c(0.5, 0.6 + 1e-12), "and"))
  expect_false(meets_criterion(c(0.7, 0.4), c(0.5, 0.5), "and"))
})
