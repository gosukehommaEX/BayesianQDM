test_that("psurv_fixed() uses the Gamma distribution of the hazard", {
  # Reference: Python_reference/ref_fixed_tte_qdm.py, design T1:
  # P(S_t > 0.87) and P(S_t <= 0.72) for Gamma(6.1, 151) at t0 = 3
  expect_equal(psurv_fixed(0.87, 3, 6.1, 151, lower.tail = FALSE),
               0.686259165773, tolerance = 1e-10)
  expect_equal(psurv_fixed(0.72, 3, 6.1, 151), 0.00106102273767,
               tolerance = 1e-10)
  expect_equal(psurv_fixed(0.5, 2, c(1, 2), c(3, 4), lower.tail = FALSE),
               stats::pgamma(-log(0.5) / 2, c(1, 2), c(3, 4)))
})

test_that("psurv_fixed() handles values outside (0, 1)", {
  expect_identical(psurv_fixed(c(-0.1, 0, 1, 1.2), 3, 2, 10,
                               lower.tail = FALSE), c(1, 1, 0, 0))
  expect_identical(psurv_fixed(c(-0.1, 0, 1, 1.2), 3, 2, 10), c(0, 0, 1, 1))
})
