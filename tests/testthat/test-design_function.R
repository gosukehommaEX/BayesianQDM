test_that("design_function() names the design function of each endpoint", {
  expect_identical(design_function("binary"), "design_binary")
  expect_identical(design_function("continuous"), "design_continuous")
  expect_identical(design_function("binary2"), "design_binary2")
  expect_identical(design_function("continuous2"), "design_continuous2")
  expect_identical(design_function("tte"), "design_tte")
  expect_error(design_function("other"), "Unknown endpoint.", fixed = TRUE)
})
