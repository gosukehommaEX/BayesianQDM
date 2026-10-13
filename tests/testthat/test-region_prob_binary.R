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
