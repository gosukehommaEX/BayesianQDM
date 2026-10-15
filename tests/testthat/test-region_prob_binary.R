test_that("region_prob_binary() uses pbetadiff() for posterior criteria", {
  cr <- crit_post(tv = 0.15, mav = -0.05)
  p <- region_prob_binary(cr, c(7.5, 3), c(5.5, 9), 4.5, 8.5)
  expect_equal(p$go, pbetadiff(0.15, c(7.5, 3), c(5.5, 9), 4.5, 8.5,
                               lower.tail = FALSE))
  expect_equal(p$nogo, pbetadiff(-0.05, c(7.5, 3), c(5.5, 9), 4.5, 8.5))
})

test_that("region_prob_binary() uses pbetabinomdiff() for predictive criteria", {
  cr <- criterion_predictive(0.1, m_t = 30, m_c = 20, gamma_go = 0.8,
                             gamma_nogo = 0.6)
  p <- region_prob_binary(cr, 7.5, 5.5, c(4.5, 2), c(8.5, 3))
  expect_equal(p$go, pbetabinomdiff(0.1, 30, 7.5, 5.5, 20, c(4.5, 2),
                                    c(8.5, 3), lower.tail = FALSE))
  expect_equal(p$go + p$nogo, c(1, 1), tolerance = 1e-12)
})

test_that("region_prob_binary() uses the Beta distribution for a fixed control", {
  cr <- crit_post(tv = 0.2, mav = 0.05)
  # P(pi_t > 0.4) and P(pi_t <= 0.25) for Beta(12, 10) and Beta(3, 19)
  p <- region_prob_binary(cr, c(12, 3), c(10, 19), NULL, NULL, fixed = 0.2)
  expect_equal(p$go, stats::pbeta(0.4, c(12, 3), c(10, 19),
                                  lower.tail = FALSE))
  expect_equal(p$nogo, stats::pbeta(0.25, c(12, 3), c(10, 19)))
  # Reference: Python_reference/ref_fixed_tte_qdm.py (design D8, y = 11)
  expect_equal(p$go[1], 0.915075681321, tolerance = 1e-10)
  expect_equal(p$nogo[1], 0.00168707910052, tolerance = 1e-10)
  # Thresholds beyond the range of the response rate
  q <- region_prob_binary(crit_post(tv = 0.9, mav = -0.5), 2, 3, NULL, NULL,
                          fixed = 0.2)
  expect_identical(c(q$go, q$nogo), c(0, 0))
})

test_that("region_prob_binary() gives exact predictive probabilities for a fixed control", {
  cr <- criterion_predictive(0, m_t = 10, gamma_go = 0.6, gamma_nogo = 0.6)
  p <- region_prob_binary(cr, 4, 8, NULL, NULL, fixed = 0.3)
  expect_equal(p$go, 283 / 646, tolerance = 1e-12)
  expect_equal(p$go + p$nogo, 1, tolerance = 1e-12)
})
