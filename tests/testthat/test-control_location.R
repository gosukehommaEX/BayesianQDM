test_that("control_location() gives the center of the control distribution", {
  # Binary: means of Beta(5, 9) (D3) and Beta(4, 8) (D4)
  expect_equal(control_location(ref_design("D3")), 5 / 14)
  expect_equal(control_location(ref_design("D4")), 1 / 3)
  # Continuous: hypothetical control location (C3) and posterior mean given
  # external data with mean -12 (C4)
  expect_identical(control_location(ref_design_cont("C3")), -10)
  expect_equal(control_location(ref_design_cont("C4")), -12)
  des <- design_continuous(12, 0, prior_vague(), prior_nix(0, 6, 3, 10),
                           crit_cont(),
                           external_c = borrow(data_continuous(30, -12, 20),
                                               0.4))
  expect_equal(control_location(des), (6 * 0 + 12 * -12) / 18)
})
