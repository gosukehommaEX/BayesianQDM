test_that("oc_sample_tte() matches the exact operating characteristics", {
  # Complete follow-up with a fixed control: each criterion is a threshold on
  # the total observed time, which is Gamma(30, log(2) / median)
  # (Python_reference/ref_fixed_tte_qdm.py, design T1)
  des <- ref_design_tte("T1")
  sc <- data.frame(median_t = c(6.3, 9, 12, 20))
  nsim <- 20000L
  set.seed(1)
  s <- oc_sample_tte(des, sc, nsim)
  oc <- oc_from_sample(s, des$criterion$gamma_go, des$criterion$gamma_nogo)
  ref_go <- c(2.48918966507e-09, 0.000437959295081, 0.0642164668685,
              0.899733183632)
  ref_nogo <- c(0.205961280815, 0.00511165154121, 6.01475591958e-05,
                2.15285582264e-09)
  tol <- function(p) 4.5 * sqrt(p * (1 - p) / nsim) + 0.002
  expect_true(all(abs(oc[, "Go"] - ref_go) <= tol(ref_go)))
  expect_true(all(abs(oc[, "NoGo"] - ref_nogo) <= tol(ref_nogo)))
  expect_identical(s$method, paste("simulation of 20000 studies;",
                                   "probabilities by exact computation",
                                   "(fixed control)"))
})

test_that("oc_sample_tte() computes the probabilities of each simulated study", {
  des <- ref_design_tte("T2")
  sc <- data.frame(median_t = c(6, 12), median_c = c(6, 8))
  nsim <- 15L
  set.seed(2)
  s <- oc_sample_tte(des, sc, nsim)
  expect_identical(dim(s$p_go), c(nsim, 2L))
  expect_true(s$per_scenario)
  expect_equal(s$weights, matrix(1 / nsim, nsim, 2))
  # The same random numbers, drawn by the C++ code
  set.seed(2)
  ca <- conduct_args(des$conduct)
  sim <- cpp_tte_oc_data(nsim, 20L, 20L, log(2) / c(6, 12), log(2) / c(6, 8),
                         ca$acc_time, ca$acc_cum, ca$drop_start,
                         ca$drop_hazard, ca$cut_time, ca$cut_events)
  k <- nsim + 4L
  res <- qdm_decide(des, data_tte(sim[k, 1], sim[k, 2]),
                    data_tte(sim[k, 3], sim[k, 4]))
  expect_equal(s$p_go[4, 2], res$prob_go, tolerance = 1e-12)
  expect_equal(s$p_nogo[4, 2], res$prob_nogo, tolerance = 1e-12)
  expect_identical(s$method, paste("simulation of 15 studies; probabilities",
                                   "by numerical integration"))
})

test_that("oc_sample_tte() handles the other control forms", {
  sc <- data.frame(median_t = c(6, 12))
  set.seed(3)
  s5 <- oc_sample_tte(ref_design_tte("T5"), sc, 10L)
  expect_identical(dim(s5$p_go), c(10L, 2L))
  # Only the treatment group is simulated: the C++ code uses the same
  # random numbers with no control group
  set.seed(3)
  ca <- conduct_args(conduct_example)
  sim <- cpp_tte_oc_data(10L, 30L, 0L, log(2) / c(6, 12), c(1, 1),
                         ca$acc_time, ca$acc_cum, ca$drop_start,
                         ca$drop_hazard, ca$cut_time, ca$cut_events)
  expect_true(all(sim[, 3:4] == 0))
  p <- region_prob_tte(ref_design_tte("T5"),
                       posterior_tte(prior_gamma(0.1, 1), sim[, 1], sim[, 2]),
                       list(shape = 8.1, rate = 61))
  expect_equal(as.vector(s5$p_go), p$go, tolerance = 1e-12)
  set.seed(4)
  s6 <- oc_sample_tte(ref_design_tte("T6"), sc, 5L)
  expect_true(all(s6$p_go >= 0 & s6$p_go <= 1))
  # Predictive criterion
  set.seed(5)
  s3 <- oc_sample_tte(ref_design_tte("T3", nMC = 50L), sc, 4L)
  expect_identical(s3$p_nogo, 1 - s3$p_go)
  expect_match(s3$method, "simulation of future studies \\(nMC = 50\\)")
})

test_that("oc_sample_tte() needs the conduct of the study", {
  des <- design_tte(10, 0, prior_gamma(1, 1), criterion = crit_tte(),
                    milestone = 3, fixed_c = 0.5)
  expect_error(oc_sample_tte(des, data.frame(median_t = 5), 10L),
               "The operating characteristics need the conduct of the study")
})

test_that("the simulated studies follow the model for other conducts", {
  # Closed forms for the simulation of cpp_tte_oc_data() with hazard 0.1 and
  # 30 subjects, within 4.5 simulation standard errors
  nsim <- 4000L
  lam <- 0.1
  close_to <- function(x, mu) {
    abs(mean(x) - mu) <= 4.5 * stats::sd(x) / sqrt(length(x))
  }
  sim_conduct <- function(conduct, seed) {
    ca <- conduct_args(conduct)
    set.seed(seed)
    cpp_tte_oc_data(nsim, 30L, 0L, lam, 1, ca$acc_time, ca$acc_cum,
                    ca$drop_start, ca$drop_hazard, ca$cut_time,
                    ca$cut_events)
  }
  # All subjects enrolled at once and an analysis at the 10th event: every
  # study has 10 events, and the total observed time is Gamma(10, lam)
  sim <- sim_conduct(conduct_tte(c(0, 1e-9), cut_events = 10), 21)
  expect_true(all(sim[, 1] == 10))
  expect_true(close_to(sim[, 2], 10 / lam))
  # Analysis at the 30th event with piecewise exponential dropout (hazard
  # 0.1 before time 2 and 0.5 after): a study with fewer than 30 events is
  # followed to the end, so the number of events and the total observed
  # time are the sums of 1(T <= C) and min(T, C) over the subjects
  sim <- sim_conduct(conduct_tte(c(0, 1), dropout_time = c(0, 2),
                                 dropout_hazard = c(0.1, 0.5),
                                 cut_events = 30), 22)
  e_min <- (1 - exp(-2 * (lam + 0.1))) / (lam + 0.1) +
    exp(-0.2 - 2 * lam) / (lam + 0.5)
  expect_true(any(sim[, 1] < 30))
  expect_true(close_to(sim[, 1], 30 * lam * e_min))
  expect_true(close_to(sim[, 2], 30 * e_min))
  # Uniform enrollment over 10 time units and an analysis at time 5:
  # subjects enrolled after the analysis contribute nothing
  sim <- sim_conduct(conduct_tte(c(0, 10), cut_time = 5), 23)
  e_r <- (5 - (1 - exp(-5 * lam)) / lam) / (10 * lam)
  expect_true(close_to(sim[, 1], 30 * lam * e_r))
  expect_true(close_to(sim[, 2], 30 * e_r))
})
