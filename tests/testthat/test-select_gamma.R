test_that("select_gamma() returns the smallest threshold meeting the target", {
  grid <- c(0.1, 0.2, 0.3, 0.4)
  expect_identical(select_gamma(grid, c(0.5, 0.2, 0.05, 0.01), 0.05), 0.3)
  expect_identical(select_gamma(grid, c(0.5, 0.2, 0.05 + 1e-13, 0.01), 0.05),
                   0.3)
  expect_identical(select_gamma(grid, c(0.5, 0.2, 0.1, 0.06), 0.05),
                   NA_real_)
  expect_identical(select_gamma(grid, c(0.01, 0, 0, 0), 0.05), 0.1)
})
