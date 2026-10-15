test_that("data_continuous2() stores the summary statistics", {
  S <- matrix(c(4, 1.2, 1.2, 2), 2)
  d <- data_continuous2(n = 12, mean = c(3.2, 1.5), cov = S)
  expect_s3_class(d, c("qdm_data_continuous2", "qdm_data"), exact = TRUE)
  expect_identical(d$n, 12L)
  expect_identical(d$mean, c(3.2, 1.5))
  expect_identical(d$cov, S)
})

test_that("data_continuous2() rejects invalid summaries", {
  expect_error(data_continuous2(2, c(0, 0), diag(2)), "'n' must be at least 3.")
  expect_error(data_continuous2(10, 0, diag(2)), "'mean' must have length 2.")
  expect_error(data_continuous2(10, c(0, 0), matrix(c(1, 1, 1, 1), 2)),
               "'cov' must be positive definite.")
})

test_that("data_continuous2() drops dimnames and checks symmetry", {
  S <- matrix(c(4, 1, 1, 2), 2, dimnames = list(c("a", "b"), c("a", "b")))
  expect_null(dimnames(data_continuous2(10, c(0, 0), S)$cov))
  expect_error(data_continuous2(10, c(0, 0), matrix(c(4, 1, 0, 2), 2)),
               "'cov' must be symmetric.")
})
