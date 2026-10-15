# References: Python_reference/ref_continuous_qdm.py (two numerical
# integrations of the whole real line, and the moment-matching formula).
post <- function(mu, kappa, nu, sigma2) {
  list(mu = mu, kappa = kappa, nu = nu, sigma2 = sigma2)
}
crit_c5 <- function() {
  criterion_predictive(10, m_t = 40, m_c = 30, gamma_go = 0.8,
                       gamma_nogo = 0.7)
}
crit_c6 <- function() {
  criterion_predictive(5, m_t = 30, m_c = 20, gamma_go = 0.7,
                       gamma_nogo = 0.6)
}

test_that("region_prob_continuous() matches the references (posterior)", {
  # Designs C1 and C3, as one vectorized call
  pt <- post(c(18, 25), 12, 11, c(625, 400))
  pc <- post(c(-12, -10), 12, 11, c(625, 600))
  ni <- region_prob_continuous(crit_cont(), pt, pc, "NI")
  expect_equal(ni$go, c(0.5, 0.697142834497), tolerance = 1e-9)
  expect_equal(ni$nogo, c(0.0878402872647, 0.0248900866027),
               tolerance = 1e-9)
  mm <- region_prob_continuous(crit_cont(), pt, pc, "MM")
  expect_equal(mm$go, c(0.5, 0.697418373575), tolerance = 1e-9)
  expect_equal(mm$nogo, c(0.0877885562533, 0.0249788420729),
               tolerance = 1e-9)
})

test_that("region_prob_continuous() matches the references (predictive)", {
  # The squared scale of the predictive distribution of the mean of m
  # future subjects is sigma^2 (1 / kappa + 1 / m)
  c5 <- region_prob_continuous(crit_c5(), post(18, 12, 11, 625),
                               post(-12, 12, 11, 625), "NI")
  expect_equal(c5$go, 0.938493080615, tolerance = 1e-9)
  expect_equal(c5$nogo, 0.0615069193852, tolerance = 1e-9)
  # Hypothetical control of design C6: sigma^2 = r sigma_nt^2 = 2 * 144
  c6 <- region_prob_continuous(crit_c6(), post(9, 15, 14, 144),
                               post(0, 15, 14, 288), "NI")
  expect_equal(c6$go, 0.7096513743, tolerance = 1e-9)
  expect_equal(c6$nogo, 0.2903486257, tolerance = 1e-9)
  mm <- region_prob_continuous(crit_c6(), post(9, 15, 14, 144),
                               post(0, 15, 14, 288), "MM")
  expect_equal(mm$go, 0.709872074443, tolerance = 1e-9)
  expect_equal(mm$nogo, 0.290127925557, tolerance = 1e-9)
})

test_that("region_prob_continuous() estimates both regions from the same draws", {
  n_mc <- 20000L
  se <- function(p) sqrt(p * (1 - p) / n_mc)
  set.seed(11)
  pred <- region_prob_continuous(crit_c5(), post(c(18, 18), 12, 11, 625),
                                 post(-12, 12, 11, 625), "MC", n_mc)
  # Predictive criterion: the two probabilities add up to 1 for every draw
  expect_equal(pred$go + pred$nogo, c(1, 1), tolerance = 1e-14)
  expect_lt(max(abs(pred$go - 0.938493080615)), 4.5 * se(0.9385))
  # Posterior criterion: the regions are disjoint
  set.seed(12)
  pst <- region_prob_continuous(crit_cont(), post(c(18, 25), 12, 11,
                                                  c(625, 400)),
                                post(c(-12, -10), 12, 11, c(625, 600)), "MC",
                                n_mc)
  expect_true(all(pst$go + pst$nogo <= 1))
  expect_lt(max(abs(pst$go - c(0.5, 0.697142834497))), 4.5 * se(0.5))
  expect_lt(max(abs(pst$nogo - c(0.0878402872647, 0.0248900866027))),
            4.5 * se(0.0878))
  # Reproducible with the seed
  set.seed(12)
  again <- region_prob_continuous(crit_cont(), post(c(18, 25), 12, 11,
                                                    c(625, 400)),
                                  post(c(-12, -10), 12, 11, c(625, 600)),
                                  "MC", n_mc)
  expect_identical(again, pst)
})

test_that("region_prob_continuous() uses the t distribution for a fixed control", {
  # Reference: Python_reference/ref_fixed_tte_qdm.py, designs C7 and C8
  p <- region_prob_continuous(crit_cont(), post(25, 12, 11, 400), NULL,
                              "NI", fixed = -10)
  expect_equal(p$go, 0.797515237746, tolerance = 1e-10)
  expect_equal(p$nogo, 0.00264736615439, tolerance = 1e-10)
  # The method is not used
  expect_identical(region_prob_continuous(crit_cont(), post(25, 12, 11, 400),
                                          NULL, "MC", nMC = 10L,
                                          fixed = -10), p)
  cr <- criterion_predictive(5, m_t = 30, gamma_go = 0.7, gamma_nogo = 0.6)
  q <- region_prob_continuous(cr, post(c(9, 0), 15, 14, 144), NULL, "MM",
                              fixed = 0)
  expect_equal(q$go[1], 0.845153264581, tolerance = 1e-10)
  expect_equal(q$go + q$nogo, c(1, 1), tolerance = 1e-12)
  expect_equal(q$go[2], stats::pt(5 / sqrt(144 * (1 / 15 + 1 / 30)), 14,
                                  lower.tail = FALSE))
})
