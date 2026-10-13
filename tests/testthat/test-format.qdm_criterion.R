test_that("format.qdm_criterion() describes a one-endpoint criterion", {
  cr <- criterion_posterior(0.2, 0.1, gamma_go = 0.5, gamma_nogo = 0.8)
  expect_identical(format(cr),
                   c("Posterior probability criterion",
                     "  Go:   P(theta > 0.2 | data) >= 0.5",
                     "  NoGo: P(theta <= 0.1 | data) >= 0.8"))
  cp <- criterion_predictive(0, m_t = 50, m_c = 40, gamma_go = 0.9,
                             gamma_nogo = 0.7)
  expect_identical(format(cp)[1L], paste("Posterior predictive probability",
                                         "criterion (m_t = 50, m_c = 40)"))
  expect_identical(format(cp)[2L],
                   "  Go:   P(future theta > 0 | data) >= 0.9")
})

test_that("format.qdm_criterion() describes a two-endpoint criterion", {
  cr <- criterion_posterior(c(0.2, 0.15), c(0.1, 0.05), 0.3, 0.4,
                            go_regions = c(1, 4))
  expect_identical(format(cr)[2:4],
                   c("  Go:   P(regions {1, 4} | data) >= 0.3",
                     "  NoGo: P(regions {9} | data) >= 0.4",
                     "  TV = (0.20, 0.15), MAV = (0.10, 0.05)"))
})

test_that("format.qdm_criterion() describes a two-endpoint predictive criterion", {
  cp <- criterion_predictive(c(0, 0.1), m_t = 30, gamma_go = 0.6,
                             gamma_nogo = 0.8, nogo_regions = c(2, 4))
  expect_identical(format(cp),
                   c(paste("Posterior predictive probability criterion",
                           "(m_t = 30, m_c = 30)"),
                     "  Go:   P(regions {1} | data) >= 0.6",
                     "  NoGo: P(regions {2, 4} | data) >= 0.8",
                     "  NULL = (0.0, 0.1)"))
})
