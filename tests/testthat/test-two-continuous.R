# Tests for two continuous endpoint functions
# Note: pbayespostpred2cont and pbayesdecisionprob2cont are computationally
#       intensive. nMC is required even for method = 'MM' (used internally
#       for predictive or validation); set to minimum values.
#       skip_on_cran() is applied to the most expensive tests.

# Shared test data
.S_small <- matrix(c(15, 3, 3, 8), 2, 2)
.Sigma   <- matrix(c(4.0, 1.2, 1.2, 2.0), 2, 2)

# ---------------------------------------------------------------------------
# pbayespostpred2cont
# nMC is always required (default = 10000L); method = 'MM' also uses nMC
# ---------------------------------------------------------------------------

test_that("pbayespostpred2cont posterior vague MC returns 9 named probs", {
  set.seed(1)
  result <- pbayespostpred2cont(
    prob = 'posterior', design = 'controlled', prior = 'vague',
    theta.TV1 = 1.5, theta.MAV1 = 0.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    n1 = 15L, n2 = 15L,
    ybar1 = c(3.0, 2.0), S1 = .S_small,
    ybar2 = c(1.0, 0.8), S2 = .S_small,
    m1 = NULL, m2 = NULL,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
    r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, Se1 = NULL, ybar_e2 = NULL, Se2 = NULL,
    nMC = 100L, method = 'MC'
  )
  expect_type(result, "double")
  expect_length(result, 9L)
  expect_named(result, paste0("R", 1:9))
  expect_true(all(result >= 0))
  expect_equal(sum(result), 1, tolerance = 1e-6)
})

test_that("pbayespostpred2cont posterior vague MM returns 9 named probs", {
  result <- pbayespostpred2cont(
    prob = 'posterior', design = 'controlled', prior = 'vague',
    theta.TV1 = 1.5, theta.MAV1 = 0.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    n1 = 15L, n2 = 15L,
    ybar1 = c(3.0, 2.0), S1 = .S_small,
    ybar2 = c(1.0, 0.8), S2 = .S_small,
    m1 = NULL, m2 = NULL,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
    r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, Se1 = NULL, ybar_e2 = NULL, Se2 = NULL,
    nMC = 100L, method = 'MM'
  )
  expect_type(result, "double")
  expect_length(result, 9L)
  expect_named(result, paste0("R", 1:9))
  expect_true(all(result >= 0))
  expect_equal(sum(result), 1, tolerance = 1e-6)
})

test_that("pbayespostpred2cont MC and MM posterior results are close", {
  set.seed(2)
  p_mc <- pbayespostpred2cont(
    prob = 'posterior', design = 'controlled', prior = 'vague',
    theta.TV1 = 1.5, theta.MAV1 = 0.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    n1 = 20L, n2 = 20L,
    ybar1 = c(3.0, 2.0), S1 = .S_small * 1.5,
    ybar2 = c(1.0, 0.8), S2 = .S_small * 1.5,
    m1 = NULL, m2 = NULL,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
    r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, Se1 = NULL, ybar_e2 = NULL, Se2 = NULL,
    nMC = 2000L, method = 'MC'
  )
  p_mm <- pbayespostpred2cont(
    prob = 'posterior', design = 'controlled', prior = 'vague',
    theta.TV1 = 1.5, theta.MAV1 = 0.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    n1 = 20L, n2 = 20L,
    ybar1 = c(3.0, 2.0), S1 = .S_small * 1.5,
    ybar2 = c(1.0, 0.8), S2 = .S_small * 1.5,
    m1 = NULL, m2 = NULL,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
    r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, Se1 = NULL, ybar_e2 = NULL, Se2 = NULL,
    nMC = 2000L, method = 'MM'
  )
  # R1 (best region) should agree within 10%
  expect_equal(p_mc["R1"], p_mm["R1"], tolerance = 0.10)
})

test_that("pbayespostpred2cont predictive vague MC returns 4 named probs", {
  set.seed(3)
  result <- pbayespostpred2cont(
    prob = 'predictive', design = 'controlled', prior = 'vague',
    theta.TV1 = NULL, theta.MAV1 = NULL,
    theta.TV2 = NULL, theta.MAV2 = NULL,
    theta.NULL1 = 1.0, theta.NULL2 = 0.5,
    n1 = 15L, n2 = 15L,
    ybar1 = c(3.0, 2.0), S1 = .S_small,
    ybar2 = c(1.0, 0.8), S2 = .S_small,
    m1 = 30L, m2 = 30L,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
    r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, Se1 = NULL, ybar_e2 = NULL, Se2 = NULL,
    nMC = 100L, method = 'MC'
  )
  expect_type(result, "double")
  expect_length(result, 4L)
  expect_named(result, paste0("R", 1:4))
  expect_true(all(result >= 0))
  expect_equal(sum(result), 1, tolerance = 1e-6)
})

test_that("pbayespostpred2cont input validation works", {
  # TV1 <= MAV1 for posterior
  expect_error(pbayespostpred2cont(
    prob = 'posterior', design = 'controlled', prior = 'vague',
    theta.TV1 = 0.5, theta.MAV1 = 1.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    n1 = 15L, n2 = 15L,
    ybar1 = c(3.0, 2.0), S1 = .S_small,
    ybar2 = c(1.0, 0.8), S2 = .S_small,
    m1 = NULL, m2 = NULL,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
    r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, Se1 = NULL, ybar_e2 = NULL, Se2 = NULL,
    nMC = 100L, method = 'MC'
  ))
})

test_that("pbayespostpred2cont posterior uncontrolled vague MM returns 9 named probs", {
  result <- pbayespostpred2cont(
    prob = 'posterior', design = 'uncontrolled', prior = 'vague',
    theta.TV1 = 1.5, theta.MAV1 = 0.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    n1 = 15L, n2 = NULL,
    ybar1 = c(3.0, 2.0), S1 = .S_small,
    ybar2 = NULL, S2 = NULL,
    m1 = NULL, m2 = NULL,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = c(1.0, 0.8), Lambda02 = NULL,
    r = 1.0,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, Se1 = NULL, ybar_e2 = NULL, Se2 = NULL,
    nMC = 100L, method = 'MM'
  )
  expect_type(result, "double")
  expect_length(result, 9L)
  expect_named(result, paste0("R", 1:9))
  expect_true(all(result >= 0))
  expect_equal(sum(result), 1, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# pbayesdecisionprob2cont
# ---------------------------------------------------------------------------

test_that("pbayesdecisionprob2cont posterior vague MM returns correct class", {
  skip_on_cran()
  set.seed(1)
  result <- pbayesdecisionprob2cont(
    nsim = 5L,
    prob = 'posterior', design = 'controlled', prior = 'vague',
    GoRegions = 1L, NoGoRegions = 9L,
    gamma1 = 0.60, gamma2 = 0.60,
    theta.TV1 = 1.5, theta.MAV1 = 0.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    n1 = 15L, n2 = 15L, m1 = NULL, m2 = NULL,
    mu1 = matrix(c(3.0, 2.0, 3.5, 2.5), nrow = 2, byrow = TRUE),
    Sigma1 = .Sigma,
    mu2 = matrix(c(1.0, 0.8, 1.0, 0.8), nrow = 2, byrow = TRUE),
    Sigma2 = .Sigma,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
    r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, Se1 = NULL, ybar_e2 = NULL, Se2 = NULL,
    method = 'MM', nMC = 100L,
    error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 1
  )
  expect_s3_class(result, "pbayesdecisionprob2cont")
  df <- as.data.frame(result)
  expect_true(all(c("mu1_ep1", "mu1_ep2", "Go", "Gray", "NoGo") %in% names(df)))
  expect_equal(nrow(df), 2L)
  expect_true(all(df$Go  >= 0 & df$Go  <= 1))
  expect_true(all(df$NoGo >= 0 & df$NoGo <= 1))
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

test_that("pbayesdecisionprob2cont predictive vague MM returns correct class", {
  skip_on_cran()
  set.seed(2)
  result <- pbayesdecisionprob2cont(
    nsim = 5L,
    prob = 'predictive', design = 'controlled', prior = 'vague',
    GoRegions = 1L, NoGoRegions = 4L,
    gamma1 = 0.60, gamma2 = 0.60,
    theta.TV1 = NULL, theta.MAV1 = NULL,
    theta.TV2 = NULL, theta.MAV2 = NULL,
    theta.NULL1 = 1.0, theta.NULL2 = 0.5,
    n1 = 15L, n2 = 15L, m1 = 30L, m2 = 30L,
    mu1 = matrix(c(3.0, 2.0, 3.5, 2.5), nrow = 2, byrow = TRUE),
    Sigma1 = .Sigma,
    mu2 = matrix(c(1.0, 0.8, 1.0, 0.8), nrow = 2, byrow = TRUE),
    Sigma2 = .Sigma,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
    r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, Se1 = NULL, ybar_e2 = NULL, Se2 = NULL,
    method = 'MM', nMC = 100L,
    error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 2
  )
  expect_s3_class(result, "pbayesdecisionprob2cont")
  df <- as.data.frame(result)
  expect_equal(nrow(df), 2L)
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

test_that("pbayesdecisionprob2cont input validation works", {
  # GoRegions and NoGoRegions overlap
  expect_error(pbayesdecisionprob2cont(
    nsim = 5L,
    prob = 'posterior', design = 'controlled', prior = 'vague',
    GoRegions = c(1L, 2L), NoGoRegions = c(2L, 9L),
    gamma1 = 0.60, gamma2 = 0.60,
    theta.TV1 = 1.5, theta.MAV1 = 0.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    n1 = 15L, n2 = 15L, m1 = NULL, m2 = NULL,
    mu1 = c(3.0, 2.0), Sigma1 = .Sigma,
    mu2 = c(1.0, 0.8), Sigma2 = .Sigma,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
    r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, Se1 = NULL, ybar_e2 = NULL, Se2 = NULL,
    method = 'MM', nMC = 100L,
    error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 1
  ))
})

test_that("pbayesdecisionprob2cont posterior uncontrolled vague MM returns correct class", {
  skip_on_cran()
  result <- pbayesdecisionprob2cont(
    nsim = 5L,
    prob = 'posterior', design = 'uncontrolled', prior = 'vague',
    GoRegions = 1L, NoGoRegions = 9L,
    gamma1 = 0.60, gamma2 = 0.60,
    theta.TV1 = 1.5, theta.MAV1 = 0.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    n1 = 15L, n2 = NULL, m1 = NULL, m2 = NULL,
    mu1 = matrix(c(3.0, 2.0, 3.5, 2.5), nrow = 2, byrow = TRUE),
    Sigma1 = .Sigma,
    mu2 = NULL, Sigma2 = NULL,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = c(1.0, 0.8), Lambda02 = NULL,
    r = 1.0,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, Se1 = NULL, ybar_e2 = NULL, Se2 = NULL,
    method = 'MM', nMC = 100L,
    error_if_Miss = FALSE, Gray_inc_Miss = FALSE, seed = 1
  )
  expect_s3_class(result, "pbayesdecisionprob2cont")
  df <- as.data.frame(result)
  expect_true(all(c("mu1_ep1", "mu1_ep2", "Go", "Gray", "NoGo") %in% names(df)))
  expect_false(any(c("mu2_ep1", "mu2_ep2") %in% names(df)))
  expect_equal(nrow(df), 2L)
  expect_true(all(df$Go  >= 0 & df$Go  <= 1))
  expect_true(all(df$NoGo >= 0 & df$NoGo <= 1))
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

# ---------------------------------------------------------------------------
# getgamma2cont
# ---------------------------------------------------------------------------

.Sigma_gc <- matrix(c(4.0, 0.8, 0.8, 1.0), 2, 2)

test_that("getgamma2cont posterior vague MM controlled returns correct class and structure", {
  result <- getgamma2cont(
    nsim = 50L, prob = 'posterior', design = 'controlled',
    prior = 'vague',
    GoRegions = 1L, NoGoRegions = 9L,
    mu1 = c(2.0, 1.0), Sigma1 = .Sigma_gc,
    mu2 = c(0.0, 0.0), Sigma2 = .Sigma_gc,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 10L, n2 = 10L,
    theta.TV1 = 1.5, theta.MAV1 = 0.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    m1 = NULL, m2 = NULL,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
    r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
    nMC = NULL, method = 'MM',
    gamma1_grid = seq(0.05, 0.95, by = 0.05),
    gamma2_grid = seq(0.05, 0.95, by = 0.05),
    seed = 1L
  )
  expect_s3_class(result, "getgamma2cont")
  expect_true(all(c("gamma1", "gamma2", "PrGo_at_gamma", "PrNoGo_at_gamma",
                    "gamma1_grid", "gamma2_grid",
                    "PrGo_grid", "PrNoGo_grid") %in% names(result)))
})

test_that("getgamma2cont posterior vague MM PrGo_grid and PrNoGo_grid in [0, 1]", {
  result <- getgamma2cont(
    nsim = 50L, prob = 'posterior', design = 'controlled',
    prior = 'vague',
    GoRegions = 1L, NoGoRegions = 9L,
    mu1 = c(2.0, 1.0), Sigma1 = .Sigma_gc,
    mu2 = c(0.0, 0.0), Sigma2 = .Sigma_gc,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 10L, n2 = 10L,
    theta.TV1 = 1.5, theta.MAV1 = 0.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    m1 = NULL, m2 = NULL,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
    r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
    nMC = NULL, method = 'MM',
    gamma1_grid = seq(0.05, 0.95, by = 0.05),
    gamma2_grid = seq(0.05, 0.95, by = 0.05),
    seed = 1L
  )
  expect_true(all(result$PrGo_grid >= 0 & result$PrGo_grid <= 1))
  expect_true(all(result$PrNoGo_grid >= 0 & result$PrNoGo_grid <= 1))
  expect_equal(dim(result$PrGo_grid),
               c(length(result$gamma1_grid), length(result$gamma2_grid)))
})

test_that("getgamma2cont posterior vague MM gamma values in (0, 1) or NA", {
  result <- getgamma2cont(
    nsim = 50L, prob = 'posterior', design = 'controlled',
    prior = 'vague',
    GoRegions = 1L, NoGoRegions = 9L,
    mu1 = c(2.0, 1.0), Sigma1 = .Sigma_gc,
    mu2 = c(0.0, 0.0), Sigma2 = .Sigma_gc,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 10L, n2 = 10L,
    theta.TV1 = 1.5, theta.MAV1 = 0.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    m1 = NULL, m2 = NULL,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
    r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
    nMC = NULL, method = 'MM',
    gamma1_grid = seq(0.05, 0.95, by = 0.05),
    gamma2_grid = seq(0.05, 0.95, by = 0.05),
    seed = 1L
  )
  if (!is.na(result$gamma1)) expect_true(result$gamma1 > 0 && result$gamma1 < 1)
  if (!is.na(result$gamma2)) expect_true(result$gamma2 > 0 && result$gamma2 < 1)
})

test_that("getgamma2cont posterior uncontrolled vague MM returns correct class", {
  result <- getgamma2cont(
    nsim = 50L, prob = 'posterior', design = 'uncontrolled',
    prior = 'vague',
    GoRegions = 1L, NoGoRegions = 9L,
    mu1 = c(2.0, 1.0), Sigma1 = .Sigma_gc,
    mu2 = NULL, Sigma2 = NULL,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 10L, n2 = NULL,
    theta.TV1 = 1.5, theta.MAV1 = 0.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    m1 = NULL, m2 = NULL,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = c(0.0, 0.0), Lambda02 = NULL,
    r = 1.0,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
    nMC = NULL, method = 'MM',
    gamma1_grid = seq(0.05, 0.95, by = 0.05),
    gamma2_grid = seq(0.05, 0.95, by = 0.05),
    seed = 2L
  )
  expect_s3_class(result, "getgamma2cont")
})

test_that("getgamma2cont predictive vague MC returns correct class", {
  result <- getgamma2cont(
    nsim = 30L, prob = 'predictive', design = 'controlled',
    prior = 'vague',
    GoRegions = 1L, NoGoRegions = 4L,
    mu1 = c(2.0, 1.0), Sigma1 = .Sigma_gc,
    mu2 = c(0.0, 0.0), Sigma2 = .Sigma_gc,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 10L, n2 = 10L,
    theta.TV1 = NULL, theta.MAV1 = NULL,
    theta.TV2 = NULL, theta.MAV2 = NULL,
    theta.NULL1 = 1.0, theta.NULL2 = 0.5,
    m1 = 20L, m2 = 20L,
    kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
    kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
    r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
    ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
    nMC = 50L, method = 'MC',
    gamma1_grid = seq(0.05, 0.95, by = 0.05),
    gamma2_grid = seq(0.05, 0.95, by = 0.05),
    seed = 3L
  )
  expect_s3_class(result, "getgamma2cont")
  if (!is.na(result$PrGo_at_gamma))
    expect_true(result$PrGo_at_gamma >= 0 && result$PrGo_at_gamma <= 1)
  if (!is.na(result$PrNoGo_at_gamma))
    expect_true(result$PrNoGo_at_gamma >= 0 && result$PrNoGo_at_gamma <= 1)
})

test_that("getgamma2cont input validation: invalid nsim", {
  expect_error(getgamma2cont(
    nsim = 0L, prob = 'posterior', design = 'controlled',
    prior = 'vague',
    GoRegions = 1L, NoGoRegions = 9L,
    mu1 = c(2.0, 1.0), Sigma1 = .Sigma_gc,
    mu2 = c(0.0, 0.0), Sigma2 = .Sigma_gc,
    target_go = 0.05, target_nogo = 0.20,
    n1 = 10L, n2 = 10L,
    theta.TV1 = 1.5, theta.MAV1 = 0.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3
  ))
})

test_that("getgamma2cont input validation: overlapping GoRegions and NoGoRegions", {
  expect_error(getgamma2cont(
    nsim = 50L, prob = 'posterior', design = 'controlled',
    prior = 'vague',
    GoRegions = c(1L, 2L), NoGoRegions = c(2L, 9L),
    mu1 = c(2.0, 1.0), Sigma1 = .Sigma_gc,
    mu2 = c(0.0, 0.0), Sigma2 = .Sigma_gc,
    target_go = 0.05, target_nogo = 0.20,
    n1 = 10L, n2 = 10L,
    theta.TV1 = 1.5, theta.MAV1 = 0.5,
    theta.TV2 = 1.0, theta.MAV2 = 0.3,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    seed = 1L
  ))
})
