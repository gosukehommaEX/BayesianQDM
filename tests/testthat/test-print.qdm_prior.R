test_that("print.qdm_prior() prints the description and returns invisibly", {
  p <- prior_beta(1, 2)
  expect_output(res <- withVisible(print(p)), "Prior: Beta\\(a = 1, b = 2\\)")
  expect_false(res$visible)
  expect_identical(res$value, p)
})
