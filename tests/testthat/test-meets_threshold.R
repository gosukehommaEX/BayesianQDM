test_that("meets_threshold() compares with a tolerance for ties", {
  expect_true(meets_threshold(0.5, 0.5))
  expect_true(meets_threshold(0.5 - 1e-12, 0.5))
  expect_false(meets_threshold(0.5 - 1e-9, 0.5))
  expect_true(meets_threshold(0.7, 0.5))
  expect_identical(meets_threshold(c(0.1, 0.3, 0.29999999999999), 0.3),
                   c(FALSE, TRUE, TRUE))
  m <- matrix(c(0.2, 0.6, 0.4, 0.8), 2)
  expect_identical(meets_threshold(m, 0.5), m >= 0.5)
  expect_false(meets_threshold(0.49, 0.5, tol = 0.001))
  expect_true(meets_threshold(0.4995, 0.5, tol = 0.001))
})
