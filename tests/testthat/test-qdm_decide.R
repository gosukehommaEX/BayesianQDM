# Exact references: Python_reference/ref_binary_qdm.py (rational arithmetic).
# In design D7 the probability of the Go region equals the threshold 0.5
# exactly (a tie), which meets the Go criterion.
ref_decide <- data.frame(
  design = c("D1", "D1", "D1", "D2", "D3",
             "D4", "D5", "D5", "D6", "D7"),
  n_t = c(12, 12, 12, 10, 12, 12, 10, 10, 8, 6),
  y_t = c(9, 6, 3, 7, 9, 7, 8, 4, 1, 3),
  n_c = c(12, 12, 12, 8, NA, NA, 10, 10, 8, 6),
  y_c = c(4, 5, 6, 2, NA, NA, 3, 4, 0, 3),
  a_t = c(10, 7, 4, 11, 10, 8, 9, 5, 2, 4),
  b_t = c(4, 7, 10, 6, 4, 6, 3, 7, 8, 4),
  a_c = c(5, 6, 7, 8, 5, 4, 4, 5, 1, 4),
  b_c = c(9, 8, 7, 15, 9, 8, 8, 7, 9, 4),
  prob_go = c(0.8207320756665553, 0.24509783589486162, 0.010407141009605608,
              0.7502162840687308, 0.8207320756665553, 0.5924430094360085,
              0.9402635338185402, 0.40597304686911645, 0.23621980701275141,
              0.5),
  prob_nogo = c(0.07292664584103187, 0.5578589165230211, 0.9606462168290966,
                0.09572468725900227, 0.07292664584103187, 0.22499768661993003,
                0.059736466181459745, 0.5940269531308836, 0.5205968501563143,
                0.34098075320453497),
  decision = c("Go", "Gray", "NoGo", "Go", "Go",
               "Go", "Go", "Gray", "Miss", "Go")
)

test_that("qdm_decide() matches the exact references", {
  r <- ref_decide
  for (i in seq_len(nrow(r))) {
    des <- ref_design(r$design[i])
    data_c <- if (is.na(r$n_c[i])) NULL else data_binary(r$n_c[i], r$y_c[i])
    res <- qdm_decide(des, data_binary(r$n_t[i], r$y_t[i]), data_c)
    info <- paste("case", i)
    expect_equal(res$posterior_t, c(a = r$a_t[i], b = r$b_t[i]), info = info)
    expect_equal(res$posterior_c, c(a = r$a_c[i], b = r$b_c[i]), info = info)
    expect_equal(res$prob_go, r$prob_go[i], tolerance = 1e-10, info = info)
    expect_equal(res$prob_nogo, r$prob_nogo[i], tolerance = 1e-10,
                 info = info)
    expect_identical(res$decision, r$decision[i], info = info)
  }
})

test_that("qdm_decide() matches high-precision references with Jeffreys priors", {
  # Python_reference/ref_distributions_mpmath.py: Beta(7.5, 5.5) versus
  # Beta(4.5, 8.5), the posteriors for 7/12 versus 4/12 responders
  jeff <- prior_beta(0.5, 0.5)
  d_post <- design_binary(12, 12, jeff, jeff, crit_post(tv = 0.2, mav = 0))
  res <- qdm_decide(d_post, data_binary(12, 7), data_binary(12, 4))
  expect_equal(res$prob_go, 0.57601763350543618, tolerance = 1e-10)
  expect_identical(res$decision, "Go")
  d_pred <- design_binary(12, 12, jeff, jeff,
                          criterion_predictive(0.05, m_t = 20,
                                               gamma_go = 0.7,
                                               gamma_nogo = 0.5))
  res <- qdm_decide(d_pred, data_binary(12, 7), data_binary(12, 4))
  expect_equal(res$prob_go, 0.74700508465629342, tolerance = 1e-12)
  expect_equal(res$prob_nogo, 1 - 0.74700508465629342, tolerance = 1e-12)
  expect_identical(res$decision, "Go")
})

test_that("qdm_decide() uses the sample sizes of the observed data", {
  des <- ref_design("D1")
  res <- qdm_decide(des, data_binary(20, 15), data_binary(5, 1))
  expect_equal(res$posterior_t, c(a = 16, b = 6))
  expect_equal(res$posterior_c, c(a = 2, b = 5))
  expect_equal(res$prob_go, pbetadiff(0.2, 16, 6, 2, 5, lower.tail = FALSE))
})

test_that("qdm_decide() returns the inputs and the criteria met", {
  des <- ref_design("D6")
  d_t <- data_binary(8, 1)
  d_c <- data_binary(8, 0)
  res <- qdm_decide(des, d_t, d_c)
  expect_s3_class(res, "qdm_decision", exact = TRUE)
  expect_identical(res$decision, "Miss")
  expect_true(res$met_go)
  expect_true(res$met_nogo)
  expect_identical(res$design, des)
  expect_identical(res$data_t, d_t)
  expect_identical(res$data_c, d_c)
  res_u <- qdm_decide(ref_design("D3"), data_binary(12, 2))
  expect_null(res_u$data_c)
  expect_identical(res_u$decision, "NoGo")
})

test_that("qdm_decide() validates the data", {
  des <- ref_design("D1")
  expect_error(qdm_decide(des, list(n = 12, y = 3), data_binary(12, 3)),
               "'data_t' must be created by data_binary().")
  expect_error(qdm_decide(des, data_binary(12, 3)),
               "'data_c' must be created by data_binary\\(\\) for a design")
  expect_error(qdm_decide(des, data_binary(12, 3),
                          data_continuous(12, 1, 1)),
               "'data_c' must be created by data_binary\\(\\) for a design")
  expect_error(qdm_decide(ref_design("D3"), data_binary(12, 3),
                          data_binary(12, 3)),
               "'data_c' must be NULL for a design without a concurrent")
  expect_error(qdm_decide("design", data_binary(12, 3)),
               "'design' must be created by a design_\\*\\(\\) function")
})

# Continuous endpoint. References: Python_reference/ref_continuous_qdm.py
# (closed-form posterior parameters, two numerical integrations and the
# moment-matching formula); designs and data in
# helper-designs-continuous.R.
ref_decide_cont <- data.frame(
  design = c("C1", "C2", "C3", "C4", "C5", "C6"),
  mu_t = c(18, 19.3333333333, 25, 18.4615384615, 18, 9),
  kappa_t = c(12, 15, 12, 13, 12, 15),
  nu_t = c(11, 17, 11, 14, 11, 14),
  sigma_t = c(25, 23.6610045699, 20, 22.2698238643, 25, 12),
  mu_c = c(-12, -10.1428571429, -10, -12, -12, 0),
  kappa_c = c(12, 21, 12, 12, 12, 15),
  nu_c = c(11, 23, 11, 11, 11, 14),
  sigma_c = c(25, 19.6574705664, 24.4948974278, 20.5382127221, 25,
              16.9705627485),
  go_ni = c(0.5, 0.472868606954, 0.697142834497, 0.520435789904,
            0.938493080615, 0.7096513743),
  nogo_ni = c(0.0878402872647, 0.0337620028253, 0.0248900866027,
              0.048372965683, 0.0615069193852, 0.2903486257),
  go_mm = c(0.5, 0.472839814489, 0.697418373575, 0.520460570391,
            0.938485309768, 0.709872074443),
  nogo_mm = c(0.0877885562533, 0.033819052643, 0.0249788420729,
              0.0484064726881, 0.0615146902318, 0.290127925557),
  decision = c("Gray", "Gray", "Go", "Gray", "Go", "Go")
)

test_that("qdm_decide() matches the references for continuous designs", {
  r <- ref_decide_cont
  for (i in seq_len(nrow(r))) {
    d <- ref_data_cont[[r$design[i]]]
    for (method in c("NI", "MM")) {
      res <- qdm_decide(ref_design_cont(r$design[i], method), d$t, d$c)
      info <- paste(r$design[i], method)
      expect_equal(res$posterior_t,
                   c(mu = r$mu_t[i], kappa = r$kappa_t[i], nu = r$nu_t[i],
                     sigma = r$sigma_t[i]), tolerance = 1e-10, info = info)
      expect_equal(res$posterior_c,
                   c(mu = r$mu_c[i], kappa = r$kappa_c[i], nu = r$nu_c[i],
                     sigma = r$sigma_c[i]), tolerance = 1e-10, info = info)
      go <- if (method == "NI") r$go_ni[i] else r$go_mm[i]
      nogo <- if (method == "NI") r$nogo_ni[i] else r$nogo_mm[i]
      expect_equal(res$prob_go, go, tolerance = 1e-9, info = info)
      expect_equal(res$prob_nogo, nogo, tolerance = 1e-9, info = info)
      expect_identical(res$decision, r$decision[i], info = info)
    }
  }
})

test_that("qdm_decide() reproduces Yamaguchi et al. (submitted), Section 5.2", {
  # Glaucoma example: P(Delta < 0 | data) = 0.9949 and
  # P(Delta <= -2 | data) = 0.7011 (four decimals) with vague priors
  des <- design_continuous(36, 31, prior_vague(), prior_vague(),
                           criterion_posterior(theta_tv = 0, theta_mav = -2,
                                               gamma_go = 0.5,
                                               gamma_nogo = 0.5))
  d_t <- data_continuous(36, -5.2, sqrt(15.62))
  d_c <- data_continuous(31, -2.7, sqrt(13.28))
  ni <- qdm_decide(des, d_t, d_c)
  expect_equal(ni$prob_go, 0.00510195902552982, tolerance = 1e-9)
  expect_equal(ni$prob_nogo, 0.701087911524607, tolerance = 1e-9)
  des_mm <- design_continuous(36, 31, prior_vague(), prior_vague(),
                              des$criterion, method = "MM")
  mm <- qdm_decide(des_mm, d_t, d_c)
  expect_equal(mm$prob_go, 0.00510318693692072, tolerance = 1e-9)
  expect_equal(mm$prob_nogo, 0.701095827640321, tolerance = 1e-9)
  for (res in list(ni, mm)) {
    expect_equal(round(1 - res$prob_go, 4), 0.9949)
    expect_equal(round(res$prob_nogo, 4), 0.7011)
  }
})

test_that("qdm_decide() reproduces Huang et al. (submitted), Table 3", {
  # P(delta < -2 | data) = 0.928, 0.954 and 0.970 for a0 = 0.25, 0.5, 0.75
  # (vague priors, power priors with equal weights for both groups)
  pub <- c(0.928, 0.954, 0.970)
  a0 <- c(0.25, 0.5, 0.75)
  for (i in seq_along(a0)) {
    des <- design_continuous(
      30, 30, prior_vague(), prior_vague(),
      criterion_posterior(theta_tv = 0, theta_mav = -2, gamma_go = 0.5,
                          gamma_nogo = 0.95),
      external_t = borrow(data_continuous(60, -0.7, 3.8), a0[i]),
      external_c = borrow(data_continuous(60, 2.5, 4), a0[i])
    )
    res <- qdm_decide(des, data_continuous(30, -0.5, 3.9),
                      data_continuous(30, 2.8, 4.3))
    expect_equal(round(res$prob_nogo, 3), pub[i], info = paste("a0", a0[i]))
  }
})

test_that("qdm_decide() uses the seed for Monte Carlo probabilities", {
  des <- ref_design_cont("C5", method = "MC", nMC = 20000L)
  d <- ref_data_cont$C5
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  r1 <- qdm_decide(des, d$t, d$c, seed = 1)
  expect_identical(runif(1), before)
  r2 <- qdm_decide(des, d$t, d$c, seed = 1)
  expect_identical(r1$prob_go, r2$prob_go)
  expect_equal(r1$prob_go + r1$prob_nogo, 1, tolerance = 1e-14)
  expect_lt(abs(r1$prob_go - 0.938493080615),
            4.5 * sqrt(0.9385 * 0.0615 / 20000))
  expect_identical(r1$decision, "Go")
  expect_error(qdm_decide(des, d$t, d$c, seed = "a"),
               "'seed' must be numeric.")
  expect_error(qdm_decide(des, d$t, d$c, seed = -3e9),
               "'seed' must be at least -2147483647 and at most 2147483647.")
})

test_that("qdm_decide() decides for two binary endpoints", {
  # Predictive probabilities (exact): Python_reference/ref_binary2_qdm.py,
  # case A1. With gamma_go = 0.4 and gamma_nogo = 0.3 (design B1) only the
  # Go criterion is met.
  d_t <- data_binary2(c(2, 1, 3, 4))
  d_c <- data_binary2(c(4, 2, 3, 1))
  res <- qdm_decide(ref_design_bin2("B1"), d_t, d_c)
  expect_identical(res$decision, "Go")
  expect_true(res$met_go)
  expect_false(res$met_nogo)
  expect_length(res$prob_regions, 4L)
  # Marginal rule "and" (design B2): Go needs 0.643 >= 0.6 and 0.647 >= 0.5
  res2 <- qdm_decide(ref_design_bin2("B2"), d_t, d_c)
  expect_identical(res2$decision, "Go")
  # Posterior probabilities by Monte Carlo simulation with a seed
  des <- ref_design_bin2("B5", nMC = 50000L)
  set.seed(5)
  before <- runif(1)
  set.seed(5)
  r1 <- qdm_decide(des, d_t, d_c, seed = 3)
  expect_identical(runif(1), before)
  r2 <- qdm_decide(des, d_t, d_c, seed = 3)
  expect_identical(r1$prob_regions, r2$prob_regions)
  expect_equal(sum(r1$prob_regions), 1, tolerance = 1e-12)
  expect_identical(r1$prob_go, r1$prob_regions[["R1"]])
  expect_identical(r1$prob_nogo, r1$prob_regions[["R9"]])
  expect_lt(abs(r1$prob_go - 0.430209), 4.5 * sqrt(0.43 * 0.57 / 50000))
  expect_identical(r1$decision, "Go")
})

test_that("qdm_decide() decides for two continuous endpoints", {
  # Moment matching: Python_reference/ref_continuous2_qdm.py, design K3
  # (predictive criterion), P(R1) = 0.626 and P(R4) = 0.0609
  d <- ref_data_cont2$K3
  res <- qdm_decide(ref_design_cont2("K3"), d$t, d$c)
  expect_identical(res$decision, "Go")
  expect_equal(res$prob_go, 0.626051840857435, tolerance = 1e-9)
  expect_equal(res$prob_nogo, 0.060905364287987, tolerance = 1e-9)
  expect_length(res$prob_regions, 4L)
  # Monte Carlo probabilities depend on the seed only
  des <- ref_design_cont2("K3", "MC", 4000L)
  a <- qdm_decide(des, d$t, d$c, seed = 21)
  b <- qdm_decide(des, d$t, d$c, seed = 21)
  expect_identical(a$prob_regions, b$prob_regions)
  expect_lt(abs(a$prob_go - 0.626051840857435), 0.04)
})

test_that("qdm_decide() decides with a fixed control for every endpoint", {
  # References: Python_reference/ref_fixed_tte_qdm.py
  res <- qdm_decide(ref_design("D8"), data_binary(20, 11))
  expect_identical(res$decision, "Go")
  expect_equal(res$prob_go, 0.915075681321, tolerance = 1e-10)
  expect_null(res$data_c)
  res <- qdm_decide(ref_design("D9"), data_binary(10, 3))
  expect_identical(res$decision, "Gray")
  expect_equal(res$prob_go, 283 / 646, tolerance = 1e-12)
  res <- qdm_decide(ref_design_cont("C7"), ref_data_cont$C7$t)
  expect_identical(res$decision, "Go")
  expect_equal(res$prob_go, 0.797515237746, tolerance = 1e-10)
  res <- qdm_decide(ref_design_cont("C8"), ref_data_cont$C8$t)
  expect_identical(res$decision, "Go")
  expect_equal(res$prob_go, 0.845153264581, tolerance = 1e-10)
  res <- qdm_decide(ref_design_bin2("BF1"), data_binary2(c(1, 0, 1, 1)))
  expect_identical(res$decision, "Gray")
  expect_equal(res$prob_go, 0.204310825893, tolerance = 1e-10)
  expect_equal(res$prob_nogo, 0.204310825893, tolerance = 1e-10)
  res <- qdm_decide(ref_design_cont2("KF1"), ref_data_cont2$KF1$t)
  expect_identical(res$decision, "Gray")
  expect_equal(res$prob_go, 0.448613893688, tolerance = 1e-9)
  res <- qdm_decide(ref_design_cont2("KF2"), ref_data_cont2$KF2$t)
  expect_identical(res$decision, "Go")
  expect_equal(res$prob_nogo, 0.13754578381, tolerance = 1e-9)
  expect_error(qdm_decide(ref_design("D8"), data_binary(20, 11),
                          data_binary(20, 4)),
               "'data_c' must be NULL for a design without a concurrent")
})

test_that("qdm_decide() decides for a time-to-event endpoint", {
  for (name in c("T1", "T2", "T5", "T6")) {
    d <- ref_data_tte[[name]]
    res <- qdm_decide(ref_design_tte(name), d$t, d$c)
    expect_identical(res$decision, "Go", info = name)
    expect_equal(c(res$prob_go, res$prob_nogo), unname(ref_prob_tte[[name]]),
                 tolerance = 1e-10, info = name)
  }
  # Predictive probabilities: Monte Carlo estimates that depend on the seed
  for (name in c("T3", "T4")) {
    d <- ref_data_tte[[name]]
    des <- ref_design_tte(name)
    a <- qdm_decide(des, d$t, d$c, seed = 8)
    b <- qdm_decide(des, d$t, d$c, seed = 8)
    expect_identical(a$prob_go, b$prob_go, info = name)
    ref <- ref_prob_tte[[name]][["go"]]
    expect_lt(abs(a$prob_go - ref), 4.5 * sqrt(ref * (1 - ref) / 1e4) + 0.002)
    expect_equal(a$prob_go + a$prob_nogo, 1)
  }
  expect_identical(qdm_decide(ref_design_tte("T3"), ref_data_tte$T3$t,
                              seed = 8)$decision, "Go")
})
