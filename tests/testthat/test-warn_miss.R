test_that("warn_miss() warns only for positive probabilities", {
  expect_silent(warn_miss(c(0, 0)))
  expect_warning(warn_miss(c(0, 0.1, 0.2)), "in 2 scenario\\(s\\)")
  expect_warning(warn_miss(1e-20), "gamma_go \\+ gamma_nogo > 1\\.$")
  expect_warning(warn_miss(1e-20, "joint"), "gamma_go \\+ gamma_nogo > 1\\.$")
  expect_warning(warn_miss(0.1, "and"),
                 "gamma_go \\+ gamma_nogo > 1 for at least one endpoint\\.$")
  expect_warning(warn_miss(0.1, "or"),
                 "With rule = \"or\", a Miss decision can occur for any")
})

test_that("warn_miss() ignores scenarios that are not evaluable", {
  expect_silent(warn_miss(c(0, NA)))
  expect_warning(warn_miss(c(NA, 0.1)), "in 1 scenario\\(s\\)")
})
