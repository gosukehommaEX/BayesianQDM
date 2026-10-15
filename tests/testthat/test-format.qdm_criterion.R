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

test_that("format.qdm_criterion() describes marginal rules", {
  cr <- criterion_posterior(c(0.2, 0.15), c(0.1, 0.05), c(0.6, 0.5), 0.7,
                            rule = "and")
  expect_identical(format(cr),
                   c("Posterior probability criterion, marginal rule \"and\"",
                     paste("  Go:   P(theta1 > 0.2 | data) >= 0.6 and",
                           "P(theta2 > 0.15 | data) >= 0.5"),
                     paste("  NoGo: P(theta1 <= 0.1 | data) >= 0.7 and",
                           "P(theta2 <= 0.05 | data) >= 0.7")))
  cp <- criterion_predictive(c(0, 0.1), m_t = 30, gamma_go = 0.6,
                             gamma_nogo = 0.8, rule = "or")
  expect_identical(format(cp)[1:2],
                   c(paste("Posterior predictive probability criterion",
                           "(m_t = 30, m_c = 30), marginal rule \"or\""),
                     paste("  Go:   P(future theta1 > 0 | data) >= 0.6 or",
                           "P(future theta2 > 0.1 | data) >= 0.6")))
})

test_that("format.qdm_criterion() shows the direction of two endpoints", {
  cr <- criterion_posterior(c(0.2, -0.1), c(0.1, 0), 0.5, 0.5,
                            direction = c("greater", "less"))
  expect_identical(format(cr)[4:5],
                   c("  TV = (0.2, -0.1), MAV = (0.1, 0.0)",
                     "  Favorable direction: theta1 greater, theta2 less"))
  cp <- criterion_predictive(c(0, 0), m_t = 30, gamma_go = 0.6,
                             gamma_nogo = 0.8, direction = "less")
  expect_identical(format(cp)[5L],
                   "  Favorable direction: theta1 less, theta2 less")
  # Not shown when larger effects are favorable on both endpoints, nor for
  # marginal rules, whose labels show the direction
  expect_length(format(criterion_posterior(c(0.2, 0.1), c(0.1, 0), 0.5,
                                           0.5)), 4L)
  expect_length(format(criterion_posterior(c(0.2, -0.1), c(0.1, 0), 0.5,
                                           0.5, rule = "and",
                                           direction = c("greater",
                                                         "less"))), 3L)
})

test_that("format.qdm_criterion() reverses the inequalities for one endpoint", {
  cr <- criterion_posterior(-0.1, 0, 0.6, 0.7, direction = "less")
  expect_identical(format(cr)[2:3],
                   c("  Go:   P(theta < -0.1 | data) >= 0.6",
                     "  NoGo: P(theta >= 0 | data) >= 0.7"))
})
