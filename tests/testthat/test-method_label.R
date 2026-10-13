test_that("method_label() describes the methods", {
  expect_identical(method_label("NI", 100L), "numerical integration")
  expect_identical(method_label("MM", 100L), "moment-matching approximation")
  expect_identical(method_label("MC", 2000L),
                   "Monte Carlo simulation (nMC = 2000)")
})
