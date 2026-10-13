test_that("sample_column() selects the column of a scenario", {
  expect_identical(sample_column(c(0.1, 0.2), 2L), c(0.1, 0.2))
  m <- cbind(c(1, 2), c(3, 4))
  expect_identical(sample_column(m, 2L), c(3, 4))
})
