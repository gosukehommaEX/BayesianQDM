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

test_that("control_location() gives the fixed value or the mean survival", {
  expect_identical(control_location(ref_design("D8")), 0.2)
  expect_identical(control_location(ref_design_cont("C7")), -10)
  expect_identical(control_location(ref_design_tte("T1")), 0.72)
  # Hypothetical control: Gamma(8.1, 61), E[exp(-3 lambda)] = (61 / 64)^8.1
  expect_equal(control_location(ref_design_tte("T5")), (61 / 64)^8.1)
  # External control only: Gamma(10.1, 51)
  expect_equal(control_location(ref_design_tte("T6")), (51 / 54)^10.1)
})

test_that("control_location() gives the control hazard for the hazard ratio", {
  expect_equal(control_location(ref_design_hr("H2")), log(2) / 6)
  # Mean of the Gamma(15.1, 101) and Gamma(15.1, 126) distributions
  expect_equal(control_location(ref_design_hr("H3")), 15.1 / 101)
  expect_equal(control_location(ref_design_hr("H4")), 15.1 / 126)
})

test_that("control_location() gives one center per endpoint for two", {
  # Dirichlet(0.25 + counts): the hypothetical counts (2, 1, 1, 0) give
  # (2.25, 1.25, 1.25, 0.25), and the external counts (14, 2, 3, 1) with
  # a0 = 0.5 give (7.25, 1.25, 1.75, 0.75); the response rates of endpoints
  # 1 and 2 are the sums of the cells "10" and "11", and "01" and "11"
  expect_equal(control_location(ref_design_bin2("B7")), c(1.5, 1.5) / 5)
  expect_equal(control_location(ref_design_bin2("B8")), c(2.5, 2) / 11)
  expect_identical(control_location(ref_design_bin2("BF1")), c(0.25, 0.5))
  # Hypothetical control location, and posterior mean given the external
  # control data (vague prior)
  expect_identical(control_location(ref_design_cont2("K4")), c(-10, -1))
  expect_equal(control_location(ref_design_cont2("K5")), c(-10, -1))
  expect_identical(control_location(ref_design_cont2("KF1")), c(-10, -1))
})
