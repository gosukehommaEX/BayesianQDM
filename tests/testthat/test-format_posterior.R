test_that("format_posterior() formats Beta parameters", {
  expect_identical(format_posterior(ref_design("D1"), c(a = 7.5, b = 5.5)),
                   "Beta(a = 7.5, b = 5.5)")
  expect_identical(format_posterior(ref_design("D1"), c(a = 1 / 3, b = 2),
                                    digits = 2),
                   "Beta(a = 0.33, b = 2)")
})
