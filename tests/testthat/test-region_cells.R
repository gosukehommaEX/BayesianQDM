test_that("region_cells() maps the regions to the nine-region grid", {
  expect_identical(region_cells(c(1, 5, 9), "posterior"), c(1L, 5L, 9L))
  expect_identical(region_cells(1:4, "predictive"), c(1L, 3L, 7L, 9L))
  expect_identical(region_cells(c(2, 4), "predictive"), c(3L, 9L))
})
