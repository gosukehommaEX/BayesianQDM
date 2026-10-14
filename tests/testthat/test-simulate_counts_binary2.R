test_that("simulate_counts_binary2() assigns subjects to the cells", {
  u <- rbind(c(0.05, 0.25, 0.55, 0.95), c(0.1, 0.1, 0.1, 0.1))
  x <- simulate_counts_binary2(u, c(0.2, 0.3, 0.4, 0.1))
  expect_identical(x, matrix(c(1L, 4L, 1L, 0L, 1L, 0L, 1L, 0L), 2L,
                             dimnames = list(NULL, c("00", "01", "10",
                                                     "11"))))
  set.seed(1)
  big <- simulate_counts_binary2(matrix(runif(4e4), 1L), c(0.1, 0.2, 0.3,
                                                         0.4))
  expect_equal(as.numeric(big) / 4e4, c(0.1, 0.2, 0.3, 0.4), tolerance = 0.03)
})
