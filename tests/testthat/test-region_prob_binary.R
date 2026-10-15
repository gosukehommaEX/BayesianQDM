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

test_that("region_prob_binary() reverses the regions for direction = \"less\"", {
  # Python_reference/ref_direction_qdm.py: posterior Beta(4, 10) and
  # Beta(9, 5), P(theta < -0.2) and P(theta >= 0)
  cr <- criterion_posterior(-0.2, 0, 0.6, 0.7, direction = "less")
  p <- region_prob_binary(cr, 4, 10, 9, 5)
  expect_equal(p$go, 0.820732075666555, tolerance = 1e-9)
  expect_equal(p$nogo, 0.023589985193162, tolerance = 1e-9)
  # The other tails of the distribution of the difference
  expect_equal(p$go, pbetadiff(-0.2, 4, 10, 9, 5), tolerance = 1e-10)
  expect_equal(p$nogo, pbetadiff(0, 4, 10, 9, 5, lower.tail = FALSE),
               tolerance = 1e-10)
  # Fixed control: the Beta distribution of the treatment rate
  pf <- region_prob_binary(cr, c(2, 5), c(20, 17), NULL, NULL, fixed = 0.3)
  expect_equal(pf$go, stats::pbeta(0.1, c(2, 5), c(20, 17)),
               tolerance = 1e-12)
  expect_equal(pf$nogo, stats::pbeta(0.3, c(2, 5), c(20, 17),
                                     lower.tail = FALSE), tolerance = 1e-12)
})

test_that("region_prob_binary() excludes ties from the Go region for direction = \"less\"", {
  # Python_reference/ref_direction_qdm.py (rational arithmetic)
  cl <- criterion_predictive(-0.1, m_t = 30, m_c = 20, gamma_go = 0.8,
                             gamma_nogo = 0.6, direction = "less")
  p <- region_prob_binary(cl, 8, 5, 5, 9)
  expect_equal(p$go, 0.0571131832789503, tolerance = 1e-12)
  expect_equal(p$go + p$nogo, 1, tolerance = 1e-12)
  # With NULL = 0, P(future diff < 0) = 0.1226 excludes the ties, which
  # P(future diff <= 0) = 0.1373 includes
  c0 <- criterion_predictive(0, m_t = 30, m_c = 20, gamma_go = 0.8,
                             gamma_nogo = 0.6, direction = "less")
  expect_equal(region_prob_binary(c0, 8, 5, 5, 9)$go, 0.122625537562555,
               tolerance = 1e-12)
  # Fixed control 0.3: P(Y / 10 - 0.3 < 0) = 1661 / 4522, whereas the Go
  # region of direction = "greater" has probability 283 / 646; Y = 3 is in
  # neither
  cf <- criterion_predictive(0, m_t = 10, gamma_go = 0.6, gamma_nogo = 0.6,
                             direction = "less")
  pf <- region_prob_binary(cf, 4, 8, NULL, NULL, fixed = 0.3)
  expect_equal(pf$go, 1661 / 4522, tolerance = 1e-12)
  expect_equal(pf$nogo, 1 - 1661 / 4522, tolerance = 1e-12)
})
