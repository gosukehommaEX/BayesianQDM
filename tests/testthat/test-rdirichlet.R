test_that("rdirichlet() returns a matrix of probability vectors", {
  set.seed(1)
  x <- rdirichlet(5, c(a = 1, b = 2, c = 3))
  expect_identical(dim(x), c(5L, 3L))
  expect_identical(colnames(x), c("a", "b", "c"))
  expect_equal(rowSums(x), rep(1, 5), tolerance = 1e-12)
  expect_true(all(x > 0))
  expect_identical(dim(rdirichlet(1, c(1, 1))), c(1L, 2L))
})

test_that("rdirichlet() has the Dirichlet means and variances", {
  set.seed(2)
  alpha <- c(1.5, 0.5, 2, 1)
  a0 <- sum(alpha)
  x <- rdirichlet(4e4, alpha)
  mean_true <- alpha / a0
  var_true <- alpha * (a0 - alpha) / (a0 ^ 2 * (a0 + 1))
  expect_true(all(abs(colMeans(x) - mean_true) < 4 * sqrt(var_true / 4e4)))
  expect_equal(apply(x, 2, stats::var), var_true, tolerance = 0.05)
})

test_that("rdirichlet() validates its arguments", {
  expect_error(rdirichlet(0, c(1, 1)), "'n' must be at least 1.")
  expect_error(rdirichlet(2, 1), "'alpha' must have length 2 or more.")
  expect_error(rdirichlet(2, c(1, 0)), "'alpha' must be greater than 0.")
})
