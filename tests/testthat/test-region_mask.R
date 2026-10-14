test_that("region_mask() sets one bit per region", {
  expect_identical(region_mask(1L), 1L)
  expect_identical(region_mask(c(1L, 2L)), 3L)
  expect_identical(region_mask(9L), 256L)
  expect_identical(region_mask(c(3L, 6L, 9L, 3L)), 4L + 32L + 256L)
})
