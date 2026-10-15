# References: Python_reference/ref_binary2_qdm.py (A: rational arithmetic
# over all future count vectors; C: Monte Carlo with 2e7 draws).
alpha_a1 <- rbind(c(2, 1, 3, 4) + 0.25, c(4, 2, 3, 1) + 0.25)

test_that("region_prob_binary2() computes predictive probabilities exactly", {
  cases <- list(
    list(t = c(2, 1, 3, 4), c = c(4, 2, 3, 1), m = c(4, 3), q = c(0.1, 0),
         ref = c(0.421763759173932, 0.221546787701068, 0.225275968676342,
                 0.131413484448658)),
    # Ties of the future effects with the null values 0 and 0.25
    list(t = c(1, 2, 2, 1), c = c(2, 1, 1, 2), m = c(4, 4), q = c(0, 0.25),
         ref = c(0.0734282326698303, 0.312585561275482, 0.116701772212982,
                 0.497284433841705)),
    list(t = c(0, 0, 0, 3), c = c(3, 0, 0, 0), m = c(2, 5),
         q = c(-0.2, 0.4),
         ref = c(0.838421426500593, 0.143884481702532, 0.00752583912440709,
                 0.0101682526724679))
  )
  for (k in seq_along(cases)) {
    cs <- cases[[k]]
    cr <- criterion_predictive(cs$q, m_t = cs$m[1], m_c = cs$m[2],
                               gamma_go = 0.5, gamma_nogo = 0.5)
    p <- region_prob_binary2(cr, matrix(cs$t + 0.25, 1L),
                             matrix(cs$c + 0.25, 1L), 1L, 1L, as.list(1:4),
                             10L)
    expect_equal(p[1L, ], cs$ref, tolerance = 1e-12, info = paste("case", k))
  }
  # Unions of regions and several pairs at once
  cr <- criterion_predictive(c(0.1, 0), m_t = 4, m_c = 3, gamma_go = 0.5,
                             gamma_nogo = 0.5)
  p2 <- region_prob_binary2(cr, alpha_a1, alpha_a1[2:1, ], c(1L, 2L),
                            c(1L, 1L), list(1:2, 1:4), 10L)
  expect_identical(dim(p2), c(2L, 2L))
  expect_equal(p2[1L, ], c(0.421763759173932 + 0.221546787701068, 1),
               tolerance = 1e-12)
  # Many pairs (the C++ code checks for user interrupts every 10000 pairs)
  p3 <- region_prob_binary2(cr, alpha_a1, alpha_a1, rep(1L, 10000L),
                            rep(2L, 10000L), list(1L), 10L)
  expect_equal(range(p3), rep(0.421763759173932, 2), tolerance = 1e-12)
  # Null values far outside [-1, 1]: every future effect on endpoint 1 is
  # larger and none on endpoint 2
  far <- criterion_predictive(c(-1e17, 1e17), m_t = 4, m_c = 3,
                              gamma_go = 0.5, gamma_nogo = 0.5)
  p4 <- region_prob_binary2(far, alpha_a1, alpha_a1, 1L, 2L, as.list(1:4),
                            10L)
  expect_equal(p4[1L, ], c(0, 1, 0, 0), tolerance = 1e-12)
})

test_that("region_prob_binary2() estimates posterior probabilities", {
  ref <- c(0.430209, 0.068027, 0.239321, 0.039777, 0.006862, 0.025173,
           0.101228, 0.018258, 0.071145)
  n_mc <- 200000L
  set.seed(2026)
  p <- region_prob_binary2(criterion_posterior(c(0.15, 0.15), c(0.1, 0.1),
                                               0.5, 0.5),
                           alpha_a1[1L, , drop = FALSE],
                           alpha_a1[2L, , drop = FALSE], 1L, 1L,
                           c(as.list(1:9), list(1:3, c(1L, 4L, 7L))), n_mc)
  se <- sqrt(ref * (1 - ref) * (1 / n_mc + 1 / 2e7))
  expect_true(all(abs(p[1L, 1:9] - ref) <= 4.5 * se))
  expect_equal(sum(p[1L, 1:9]), 1, tolerance = 1e-12)
  # The unions use the same draws as the regions
  expect_equal(p[1L, 10], sum(p[1L, 1:3]), tolerance = 1e-12)
  # Marginal probabilities: differences of two Beta variables
  exact1 <- pbetadiff(0.15, 7.5, 3.5, 4.5, 6.5, lower.tail = FALSE)
  exact2 <- pbetadiff(0.15, 5.5, 5.5, 3.5, 7.5, lower.tail = FALSE)
  expect_equal(exact1, 0.737498698318608, tolerance = 1e-9)
  expect_lt(abs(p[1L, 10] - exact1), 4.5 * sqrt(exact1 * (1 - exact1) / n_mc))
  expect_lt(abs(p[1L, 11] - exact2), 4.5 * sqrt(exact2 * (1 - exact2) / n_mc))
})

test_that("region_prob_binary2() handles a fixed control", {
  # Python_reference/ref_fixed_tte_qdm.py (design BF1): exact enumeration of
  # the future cell counts of the treatment group
  cr <- criterion_predictive(c(0.1, 0), m_t = 4, gamma_go = 0.4,
                             gamma_nogo = 0.3)
  p <- region_prob_binary2(cr, matrix(c(1, 0, 1, 1) + 0.25, 1L), NULL, 1L,
                           NULL, as.list(1:4), 10L, fixed = c(0.25, 0.5))
  expect_equal(p[1L, ], c(0.204310825893, 0.569126674107, 0.0222516741071,
                          0.204310825893), tolerance = 1e-10)
  # Posterior probabilities from draws of the treatment distribution only:
  # with fixed values 0, the regions are those of the treatment rates
  set.seed(2)
  cp <- criterion_posterior(c(0.5, 0.5), c(0.3, 0.3), 0.5, 0.5)
  a <- matrix(c(2, 1, 3, 4) + 0.25, 1L)
  p2 <- region_prob_binary2(cp, a, NULL, c(1L, 1L), NULL, as.list(1:9),
                            20000L, fixed = c(0, 0))
  # Both pairs use the same draws
  expect_identical(p2[1L, ], p2[2L, ])
  expect_equal(sum(p2[1L, ]), 1)
  set.seed(3)
  p3 <- region_prob_binary2(cp, a, NULL, 1L, NULL, list(1L), 20000L,
                            fixed = c(0, 0))
  set.seed(3)
  g <- matrix(stats::rgamma(4 * 20000, rep(a[1L, ], 20000)), ncol = 4,
              byrow = TRUE)
  pi1 <- (g[, 3] + g[, 4]) / rowSums(g)
  pi2 <- (g[, 2] + g[, 4]) / rowSums(g)
  expect_equal(p3[1L, 1L], mean(pi1 > 0.5 & pi2 > 0.5))
})
