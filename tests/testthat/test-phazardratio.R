test_that("phazardratio() matches numerical integration", {
  # Python_reference/ref_hazard_ratio.py: integration over the control
  # hazard and over the treatment hazard (relative difference at most
  # 1.1e-14), checked by 2e6 Monte Carlo draws
  q <- c(0.7, 1, 0.5, 2, 0.25)
  shape_t <- c(20.1, 5.5, 0.6, 12, 50)
  rate_t <- c(301, 40, 2, 100, 400)
  shape_c <- c(30.1, 9.5, 3, 4, 40)
  rate_c <- c(251, 35, 10, 60, 250)
  ref <- c(0.7923542854783, 0.9069641684317, 0.4607304126901,
           0.5320801962467, 4.761184427754e-08)
  p <- phazardratio(q, shape_t, rate_t, shape_c, rate_c)
  # Element by element, so that the small probability is checked with its
  # own relative error
  for (k in seq_along(ref)) {
    expect_equal(p[k], ref[k], tolerance = 1e-10, info = paste("case", k))
  }
  upper <- phazardratio(q, shape_t, rate_t, shape_c, rate_c,
                        lower.tail = FALSE)
  expect_equal(upper, 1 - ref, tolerance = 1e-12)
  # The upper tail is computed directly: 1 - upper keeps the small
  # probability up to the rounding of a number near 1
  expect_equal(1 - upper[5], ref[5], tolerance = 1e-7)
})

test_that("phazardratio() handles non-positive quantiles and recycling", {
  expect_identical(phazardratio(c(-1, 0), 5, 40, 9, 35), c(0, 0))
  expect_identical(phazardratio(-1, 5, 40, 9, 35, lower.tail = FALSE), 1)
  p <- phazardratio(c(0.5, 1, 2), 5, 40, 9, 35)
  expect_length(p, 3L)
  expect_true(all(diff(p) > 0))
  # The ratio of the treatment to the control hazard: exchanging the groups
  # gives the distribution of the reciprocal
  expect_equal(phazardratio(0.8, 5, 40, 9, 35),
               phazardratio(1 / 0.8, 9, 35, 5, 40, lower.tail = FALSE),
               tolerance = 1e-12)
})

test_that("phazardratio() validates its arguments", {
  expect_error(phazardratio("a", 1, 1, 1, 1), "'q' must be numeric.")
  expect_error(phazardratio(1, 0, 1, 1, 1), "'shape_t' must be greater than 0.")
  expect_error(phazardratio(1, 1, -1, 1, 1), "'rate_t' must be greater than 0.")
  expect_error(phazardratio(1, 1, 1, 0, 1), "'shape_c' must be greater than 0.")
  expect_error(phazardratio(1, 1, 1, 1, 0), "'rate_c' must be greater than 0.")
  expect_error(phazardratio(1, 1, 1, 1, 1, lower.tail = NA),
               "'lower.tail' must be TRUE or FALSE.")
  expect_error(phazardratio(c(1, 2), c(1, 2, 3), 1, 1, 1),
               "'q' must have length 1 or 3.")
})
