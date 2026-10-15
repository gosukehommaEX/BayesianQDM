test_that("outcome_prob() uses full or factored weights", {
  s1 <- list(weights = cbind(c(0.2, 0.8), c(0.5, 0.5)))
  expect_equal(outcome_prob(c(TRUE, FALSE), s1), c(0.2, 0.5))
  # Pairs of outcomes: P(pair) = w_t * w_c
  s2 <- list(weights_t = cbind(c(0.3, 0.7), c(0.6, 0.4)),
             weights_c = cbind(c(0.1, 0.9), c(0.5, 0.5)))
  ind <- matrix(c(TRUE, FALSE, TRUE, TRUE), 2L)
  expect_equal(outcome_prob(ind, s2),
               c(0.3 * 0.1 + 0.3 * 0.9 + 0.7 * 0.9,
                 0.6 * 0.5 + 0.6 * 0.5 + 0.4 * 0.5))
})
