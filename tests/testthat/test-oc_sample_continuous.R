test_that("oc_sample_continuous() returns probabilities and weights", {
  set.seed(1)
  s <- oc_sample_continuous(ref_design_cont("C1"), scen_cont_c, 50L)
  expect_identical(dim(s$p_go), c(50L, 3L))
  expect_identical(dim(s$p_nogo), c(50L, 3L))
  expect_equal(s$weights, matrix(1 / 50, 50, 3))
  expect_true(all(s$p_go >= 0 & s$p_go <= 1))
  expect_true(all(s$p_go + s$p_nogo <= 1 + 1e-12))
  expect_identical(s$method, paste("simulation of 50 studies; probabilities",
                                   "by numerical integration"))
  s_mc <- oc_sample_continuous(ref_design_cont("C3", "MC", 200L),
                               scen_cont_u, 5L)
  expect_identical(s_mc$method, paste("simulation of 5 studies; probabilities",
                                      "by Monte Carlo simulation (nMC = 200)"))
})

test_that("oc_sample_continuous() uses common random numbers", {
  sc <- data.frame(mu_t = c(10, 30, 10), sd_t = c(20, 25, 20), mu_c = 0,
                   sd_c = 25)
  set.seed(2)
  s <- oc_sample_continuous(ref_design_cont("C1", "MM"), sc, 200L)
  expect_identical(s$p_go[, 1], s$p_go[, 3])
  expect_identical(s$p_nogo[, 1], s$p_nogo[, 3])
  # A larger true effect gives larger probabilities of the Go region for
  # every simulated study when only the location changes
  sc2 <- data.frame(mu_t = c(10, 30), sd_t = 20, mu_c = 0, sd_c = 25)
  set.seed(2)
  s2 <- oc_sample_continuous(ref_design_cont("C1", "MM"), sc2, 200L)
  expect_true(all(s2$p_go[, 2] > s2$p_go[, 1]))
})

test_that("oc_sample_continuous() simulates the sufficient statistics", {
  # The draws are rnorm(nsim) and rchisq(nsim, n_t - 1) for the treatment
  # group, then the same for the control group. Each simulated study is
  # recomputed here with qdm_decide().
  nsim <- 20L
  sc <- data.frame(mu_t = c(20, 35), sd_t = c(25, 30), mu_c = c(0, -5),
                   sd_c = c(25, 20))
  des <- ref_design_cont("C2")
  set.seed(3)
  s <- oc_sample_continuous(des, sc, nsim)
  set.seed(3)
  z_t <- rnorm(nsim)
  x_t <- rchisq(nsim, 9)
  z_c <- rnorm(nsim)
  x_c <- rchisq(nsim, 7)
  for (j in 1:2) {
    for (i in c(1, 7, 20)) {
      d_t <- data_continuous(10, sc$mu_t[j] + sc$sd_t[j] * z_t[i] / sqrt(10),
                             sc$sd_t[j] * sqrt(x_t[i] / 9))
      d_c <- data_continuous(8, sc$mu_c[j] + sc$sd_c[j] * z_c[i] / sqrt(8),
                             sc$sd_c[j] * sqrt(x_c[i] / 7))
      res <- qdm_decide(des, d_t, d_c)
      expect_equal(s$p_go[i, j], res$prob_go, tolerance = 1e-12)
      expect_equal(s$p_nogo[i, j], res$prob_nogo, tolerance = 1e-12)
    }
  }
  # Without a concurrent control group only the treatment group is simulated
  des3 <- ref_design_cont("C3")
  set.seed(4)
  s3 <- oc_sample_continuous(des3, scen_cont_u[2, ], nsim)
  after <- runif(1)
  set.seed(4)
  z_t <- rnorm(nsim)
  x_t <- rchisq(nsim, 11)
  expect_identical(runif(1), after)
  d_t <- data_continuous(12, 20 + 20 * z_t[5] / sqrt(12),
                         20 * sqrt(x_t[5] / 11))
  expect_equal(s3$p_go[5, 1], qdm_decide(des3, d_t)$prob_go,
               tolerance = 1e-12)
})
