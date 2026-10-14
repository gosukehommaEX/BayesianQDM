test_that("sample_column() selects the column of a scenario", {
  expect_identical(sample_column(c(0.1, 0.2), 2L), c(0.1, 0.2))
  m <- cbind(c(1, 2), c(3, 4))
  expect_identical(sample_column(m, 2L), c(3, 4))
})

test_that("sample_column() handles lists and shared matrices", {
  m <- cbind(c(1, 2), c(3, 4))
  expect_identical(sample_column(m, 2L, per_scenario = FALSE), m)
  expect_identical(sample_column(list(m, m * 2), 2L, per_scenario = TRUE),
                   list(c(3, 4), c(6, 8)))
  expect_identical(sample_column(list(m, m), 1L, per_scenario = FALSE),
                   list(m, m))
})
