test_that("warn_miss() warns only for positive probabilities", {
  expect_silent(warn_miss(c(0, 0)))
  expect_warning(warn_miss(c(0, 0.1, 0.2)), "in 2 scenario\\(s\\)")
  expect_warning(warn_miss(1e-20), "gamma_go \\+ gamma_nogo > 1")
})
