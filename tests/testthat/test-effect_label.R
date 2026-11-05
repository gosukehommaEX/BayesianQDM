test_that("effect_label() names the effect of the endpoint", {
  expect_identical(effect_label(ref_design("D1")),
                   expression(theta == pi[t] - pi[c]))
  expect_identical(effect_label(ref_design_cont("C1")),
                   expression(theta == mu[t] - mu[c]))
  expect_identical(effect_label(ref_design_tte("T2")),
                   expression(theta == S[t](t[0]) - S[c](t[0])))
  expect_identical(effect_label(ref_design_hr("H1")),
                   expression(theta == lambda[t] / lambda[c]))
})
