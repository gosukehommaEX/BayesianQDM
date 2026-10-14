test_that("oc_sample_binary2() enumerates pairs of outcomes", {
  des <- ref_design_bin2("B1")
  s <- oc_sample_binary2(des, scen_bin2, 1L)
  n_t <- choose(3 + 3, 3)
  n_c <- choose(2 + 3, 3)
  expect_identical(dim(s$p_go), as.integer(c(n_t, n_c)))
  expect_false(s$per_scenario)
  expect_null(s[["weights"]])
  expect_identical(dim(s$weights_t), as.integer(c(n_t, 3)))
  expect_equal(colSums(s$weights_t), rep(1, 3))
  expect_equal(colSums(s$weights_c), rep(1, 3))
  # Each element is the probability for one pair of outcomes
  x_t <- outcomes_binary2(3)
  x_c <- outcomes_binary2(2)
  res <- decide_binary2(des, data_binary2(x_t[7, ]), data_binary2(x_c[4, ]))
  expect_equal(s$p_go[7, 4], res$prob_go, tolerance = 1e-12)
  expect_equal(s$p_nogo[7, 4], res$prob_nogo, tolerance = 1e-12)
  expect_identical(s$method, paste("exact enumeration of the outcomes;",
                                   "predictive probabilities computed",
                                   "exactly"))
  s2 <- oc_sample_binary2(ref_design_bin2("B4", nMC = 100L), scen_bin2_d,
                          1L)
  expect_length(s2$p_go, 2L)
  expect_identical(dim(s2$p_go[[2]]), as.integer(rep(choose(7, 3), 2)))
  expect_match(s2$method, "Monte Carlo simulation \\(nMC = 100\\)")
  # Enumeration is limited to 1e7 pairs of outcomes (3276^2 for 25 subjects
  # in each group)
  big <- design_binary2(25, 25, dir_q(), dir_q(),
                        criterion_posterior(c(0.15, 0.15), c(0.1, 0.1), 0.3,
                                            0.4))
  expect_error(oc_sample_binary2(big, scen_bin2, 1L),
               "would need 10,732,176 pairs of treatment and control")
  expect_identical(dim(oc_sample_binary2(ref_design_bin2("B7"), scen_bin2,
                                         1L)$p_go), c(20L, 1L))
})

test_that("oc_sample_binary2() maps simulated studies to their probabilities", {
  des <- ref_design_bin2("B2", enumerate = FALSE)
  nsim <- 30L
  set.seed(9)
  s <- oc_sample_binary2(des, scen_bin2, nsim)
  expect_true(s$per_scenario)
  expect_identical(dim(s$p_go[[1]]), c(nsim, 3L))
  expect_equal(s$weights, matrix(1 / nsim, nsim, 3))
  # The uniform random numbers are drawn for the treatment group and then
  # for the control group
  set.seed(9)
  u_t <- matrix(runif(nsim * 3), nrow = nsim)
  u_c <- matrix(runif(nsim * 2), nrow = nsim)
  for (j in c(1, 3)) {
    x_t <- simulate_counts_binary2(u_t, cellprob_binary2(scen_bin2$pi_t1[j],
                                                         scen_bin2$pi_t2[j],
                                                         scen_bin2$rho_t[j]))
    x_c <- simulate_counts_binary2(u_c, cellprob_binary2(scen_bin2$pi_c1[j],
                                                         scen_bin2$pi_c2[j],
                                                         scen_bin2$rho_c[j]))
    for (i in c(1, 17, 30)) {
      res <- decide_binary2(des, data_binary2(x_t[i, ]),
                            data_binary2(x_c[i, ]))
      expect_equal(c(s$p_go[[1]][i, j], s$p_go[[2]][i, j]), res$prob_go,
                   tolerance = 1e-12)
      expect_equal(c(s$p_nogo[[1]][i, j], s$p_nogo[[2]][i, j]),
                   res$prob_nogo, tolerance = 1e-12)
    }
  }
  # Without a concurrent control group only the treatment group is simulated
  set.seed(4)
  s7 <- oc_sample_binary2(ref_design_bin2("B7", enumerate = FALSE),
                          scen_bin2[, 1:3], nsim)
  after <- runif(1)
  set.seed(4)
  invisible(runif(nsim * 3))
  expect_identical(runif(1), after)
  expect_identical(dim(s7$p_go), c(nsim, 3L))
})
