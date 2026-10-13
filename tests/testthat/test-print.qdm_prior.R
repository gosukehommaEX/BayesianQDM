test_that("print.qdm_prior() prints the description and returns invisibly", {
  p <- prior_beta(1, 2)
  expect_output(print(p), "Prior: Beta\\(a = 1, b = 2\\)")
  expect_invisible(print(p))
})
