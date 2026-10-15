test_that("event_labels() describes the Go and NoGo events", {
  expect_identical(event_labels(crit_post()),
                   list(go = "P(theta > 0.2 | data)",
                        nogo = "P(theta <= 0.1 | data)"))
  cp <- criterion_predictive(0.05, m_t = 30, gamma_go = 0.9,
                             gamma_nogo = 0.7)
  expect_identical(event_labels(cp),
                   list(go = "P(future theta > 0.05 | data)",
                        nogo = "P(future theta <= 0.05 | data)"))
  c2 <- criterion_posterior(c(0.2, 0.1), c(0, 0), 0.5, 0.5,
                            go_regions = c(1, 2), nogo_regions = 9)
  expect_identical(event_labels(c2),
                   list(go = "P(regions {1, 2} | data)",
                        nogo = "P(regions {9} | data)"))
  expect_identical(event_labels(crit_post(tv = 1 / 3), digits = 2)[["go"]],
                   "P(theta > 0.33 | data)")
})

test_that("event_labels() gives one label per endpoint for marginal rules", {
  cr <- criterion_posterior(c(0.2, 0.15), c(0.1, 0), 0.5, 0.5, rule = "and")
  expect_identical(event_labels(cr),
                   list(go = c("P(theta1 > 0.2 | data)",
                               "P(theta2 > 0.15 | data)"),
                        nogo = c("P(theta1 <= 0.1 | data)",
                                 "P(theta2 <= 0 | data)")))
  cp <- criterion_predictive(c(0, 0.1), m_t = 30, gamma_go = 0.5,
                             gamma_nogo = 0.5, rule = "or")
  expect_identical(event_labels(cp)$nogo,
                   c("P(future theta1 <= 0 | data)",
                     "P(future theta2 <= 0.1 | data)"))
})
