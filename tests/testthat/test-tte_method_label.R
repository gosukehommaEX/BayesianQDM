test_that("tte_method_label() describes the computation", {
  expect_identical(tte_method_label(ref_design_tte("T1")),
                   "exact computation (fixed control)")
  expect_identical(tte_method_label(ref_design_tte("T2")),
                   "numerical integration")
  expect_identical(tte_method_label(ref_design_tte("T3", nMC = 400L)),
                   "simulation of future studies (nMC = 400)")
})
