test_that("check_regions() returns sorted integer region sets", {
  out <- check_regions(c(4, 1), 9, n_regions = 9L)
  expect_identical(out, list(go = c(1L, 4L), nogo = 9L))
})

test_that("check_regions() rejects invalid region sets", {
  expect_error(check_regions(10, 9, 9L), "'go_regions' must be at least 1")
  expect_error(check_regions(1.5, 9, 9L), "whole numbers")
  expect_error(check_regions(c(1, 1), 9, 9L), "'go_regions' must not contain")
  expect_error(check_regions(1, c(9, 9), 9L), "'nogo_regions' must not")
  expect_error(check_regions(c(1, 9), 9, 9L), "must not share any region")
})
