test_that("effect_support() gives the range of the effect", {
  expect_identical(effect_support(ref_design("D1")), c(-1, 1))
  expect_identical(effect_support(ref_design("D8")), c(-0.2, 0.8))
  expect_identical(effect_support(ref_design_cont("C1")), c(-Inf, Inf))
  expect_identical(effect_support(ref_design_tte("T2")), c(-1, 1))
  expect_identical(effect_support(ref_design_tte("T1")), c(-0.72, 1 - 0.72))
  expect_identical(effect_support(ref_design_hr("H1")), c(0, Inf))
})
