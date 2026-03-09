# Tests for two continuous endpoint functions
# Note: pbayespostpred2cont and pbayesdecisionprob2cont are computationally
#       intensive. nMC is required even for CalcMethod = 'MM' (used internally
#       for predictive or validation); set to minimum values.
#       skip_on_cran() is applied to the most expensive tests.

# Shared test data
.S_small <- matrix(c(15, 3, 3, 8), 2, 2)
.Sigma   <- matrix(c(4.0, 1.2, 1.2, 2.0), 2, 2)

# ---------------------------------------------------------------------------
# pbayespostpred2cont
# nMC is always required (default = 10000L); CalcMethod = 'MM' also uses nMC
# ---------------------------------------------------------------------------

test_that("pbayespostpred2cont posterior vague MC returns 9 named probs", {
  set.seed(1)
  result <- pbayespostpred2cont(
    prob = 'posterior', design = 'controlled', prior = 'vague',
    theta_TV1 = 1.5, theta_MAV1 = 0.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    n_t = 15L, n_c = 15L,
    ybar_t = c(3.0, 2.0), S_t = .S_small,
    ybar_c = c(1.0, 0.8), S_c = .S_small,
    m_t = NULL, m_c = NULL,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
    r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL,
    nMC = 100L, CalcMethod = 'MC'
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
    theta_TV1 = 1.5, theta_MAV1 = 0.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    n_t = 15L, n_c = 15L,
    ybar_t = c(3.0, 2.0), S_t = .S_small,
    ybar_c = c(1.0, 0.8), S_c = .S_small,
    m_t = NULL, m_c = NULL,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
    r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL,
    nMC = 100L, CalcMethod = 'MM'
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
    theta_TV1 = 1.5, theta_MAV1 = 0.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    n_t = 20L, n_c = 20L,
    ybar_t = c(3.0, 2.0), S_t = .S_small * 1.5,
    ybar_c = c(1.0, 0.8), S_c = .S_small * 1.5,
    m_t = NULL, m_c = NULL,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
    r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL,
    nMC = 2000L, CalcMethod = 'MC'
  )
  p_mm <- pbayespostpred2cont(
    prob = 'posterior', design = 'controlled', prior = 'vague',
    theta_TV1 = 1.5, theta_MAV1 = 0.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    n_t = 20L, n_c = 20L,
    ybar_t = c(3.0, 2.0), S_t = .S_small * 1.5,
    ybar_c = c(1.0, 0.8), S_c = .S_small * 1.5,
    m_t = NULL, m_c = NULL,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
    r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL,
    nMC = 2000L, CalcMethod = 'MM'
  )
  # R1 (best region) should agree within 10%
  expect_equal(p_mc["R1"], p_mm["R1"], tolerance = 0.10)
})

test_that("pbayespostpred2cont predictive vague MC returns 4 named probs", {
  set.seed(3)
  result <- pbayespostpred2cont(
    prob = 'predictive', design = 'controlled', prior = 'vague',
    theta_TV1 = NULL, theta_MAV1 = NULL,
    theta_TV2 = NULL, theta_MAV2 = NULL,
    theta_NULL1 = 1.0, theta_NULL2 = 0.5,
    n_t = 15L, n_c = 15L,
    ybar_t = c(3.0, 2.0), S_t = .S_small,
    ybar_c = c(1.0, 0.8), S_c = .S_small,
    m_t = 30L, m_c = 30L,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
    r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL,
    nMC = 100L, CalcMethod = 'MC'
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
    theta_TV1 = 0.5, theta_MAV1 = 1.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    n_t = 15L, n_c = 15L,
    ybar_t = c(3.0, 2.0), S_t = .S_small,
    ybar_c = c(1.0, 0.8), S_c = .S_small,
    m_t = NULL, m_c = NULL,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
    r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL,
    nMC = 100L, CalcMethod = 'MC'
  ))
})

test_that("pbayespostpred2cont posterior uncontrolled vague MM returns 9 named probs", {
  result <- pbayespostpred2cont(
    prob = 'posterior', design = 'uncontrolled', prior = 'vague',
    theta_TV1 = 1.5, theta_MAV1 = 0.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    n_t = 15L, n_c = NULL,
    ybar_t = c(3.0, 2.0), S_t = .S_small,
    ybar_c = NULL, S_c = NULL,
    m_t = NULL, m_c = NULL,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = c(1.0, 0.8), Lambda0_c = NULL,
    r = 1.0,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL,
    nMC = 100L, CalcMethod = 'MM'
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
    gamma_go = 0.60, gamma_nogo = 0.60,
    theta_TV1 = 1.5, theta_MAV1 = 0.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    n_t = 15L, n_c = 15L, m_t = NULL, m_c = NULL,
    mu_t = matrix(c(3.0, 2.0, 3.5, 2.5), nrow = 2, byrow = TRUE),
    Sigma_t = .Sigma,
    mu_c = matrix(c(1.0, 0.8, 1.0, 0.8), nrow = 2, byrow = TRUE),
    Sigma_c = .Sigma,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
    r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL,
    CalcMethod = 'MM', nMC = 100L,
    error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 1
  )
  expect_s3_class(result, "pbayesdecisionprob2cont")
  df <- as.data.frame(result)
  expect_true(all(c("mu_t1", "mu_t2", "Go", "Gray", "NoGo") %in% names(df)))
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
    gamma_go = 0.60, gamma_nogo = 0.60,
    theta_TV1 = NULL, theta_MAV1 = NULL,
    theta_TV2 = NULL, theta_MAV2 = NULL,
    theta_NULL1 = 1.0, theta_NULL2 = 0.5,
    n_t = 15L, n_c = 15L, m_t = 30L, m_c = 30L,
    mu_t = matrix(c(3.0, 2.0, 3.5, 2.5), nrow = 2, byrow = TRUE),
    Sigma_t = .Sigma,
    mu_c = matrix(c(1.0, 0.8, 1.0, 0.8), nrow = 2, byrow = TRUE),
    Sigma_c = .Sigma,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
    r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL,
    CalcMethod = 'MM', nMC = 100L,
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
    gamma_go = 0.60, gamma_nogo = 0.60,
    theta_TV1 = 1.5, theta_MAV1 = 0.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    n_t = 15L, n_c = 15L, m_t = NULL, m_c = NULL,
    mu_t = c(3.0, 2.0), Sigma_t = .Sigma,
    mu_c = c(1.0, 0.8), Sigma_c = .Sigma,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
    r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL,
    CalcMethod = 'MM', nMC = 100L,
    error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 1
  ))
})

test_that("pbayesdecisionprob2cont posterior uncontrolled vague MM returns correct class", {
  skip_on_cran()
  result <- pbayesdecisionprob2cont(
    nsim = 5L,
    prob = 'posterior', design = 'uncontrolled', prior = 'vague',
    GoRegions = 1L, NoGoRegions = 9L,
    gamma_go = 0.60, gamma_nogo = 0.60,
    theta_TV1 = 1.5, theta_MAV1 = 0.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    n_t = 15L, n_c = NULL, m_t = NULL, m_c = NULL,
    mu_t = matrix(c(3.0, 2.0, 3.5, 2.5), nrow = 2, byrow = TRUE),
    Sigma_t = .Sigma,
    mu_c = NULL, Sigma_c = NULL,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = c(1.0, 0.8), Lambda0_c = NULL,
    r = 1.0,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL,
    CalcMethod = 'MM', nMC = 100L,
    error_if_Miss = FALSE, Gray_inc_Miss = FALSE, seed = 1
  )
  expect_s3_class(result, "pbayesdecisionprob2cont")
  df <- as.data.frame(result)
  expect_true(all(c("mu_t1", "mu_t2", "Go", "Gray", "NoGo") %in% names(df)))
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
    mu_t_go = c(2.0, 1.0), Sigma_t_go = .Sigma_gc,
    mu_c_go = c(0.0, 0.0), Sigma_c_go = .Sigma_gc,
    mu_t_nogo = c(2.0, 1.0), Sigma_t_nogo = .Sigma_gc,
    mu_c_nogo = c(0.0, 0.0), Sigma_c_nogo = .Sigma_gc,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 10L, n_c = 10L,
    theta_TV1 = 1.5, theta_MAV1 = 0.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    m_t = NULL, m_c = NULL,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
    r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
    nMC = NULL, CalcMethod = 'MM',
    gamma_go_grid   = seq(0.05, 0.95, by = 0.05),
    gamma_nogo_grid = seq(0.05, 0.95, by = 0.05),
    seed = 1L
  )
  expect_s3_class(result, "getgamma2cont")
  expect_true(all(c("gamma_go", "gamma_nogo", "PrGo_opt", "PrNoGo_opt",
                    "target_go", "target_nogo", "grid_results") %in% names(result)))
  expect_true(all(c("gamma_grid", "PrGo_grid", "PrNoGo_grid") %in%
                    names(result$grid_results)))
})

test_that("getgamma2cont posterior vague MM PrGo_grid and PrNoGo_grid in [0, 1]", {
  result <- getgamma2cont(
    nsim = 50L, prob = 'posterior', design = 'controlled',
    prior = 'vague',
    GoRegions = 1L, NoGoRegions = 9L,
    mu_t_go = c(2.0, 1.0), Sigma_t_go = .Sigma_gc,
    mu_c_go = c(0.0, 0.0), Sigma_c_go = .Sigma_gc,
    mu_t_nogo = c(2.0, 1.0), Sigma_t_nogo = .Sigma_gc,
    mu_c_nogo = c(0.0, 0.0), Sigma_c_nogo = .Sigma_gc,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 10L, n_c = 10L,
    theta_TV1 = 1.5, theta_MAV1 = 0.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    m_t = NULL, m_c = NULL,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
    r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
    nMC = NULL, CalcMethod = 'MM',
    gamma_go_grid   = seq(0.05, 0.95, by = 0.05),
    gamma_nogo_grid = seq(0.05, 0.95, by = 0.05),
    seed = 1L
  )
  expect_true(all(result$grid_results$PrGo_grid   >= 0 & result$grid_results$PrGo_grid   <= 1))
  expect_true(all(result$grid_results$PrNoGo_grid >= 0 & result$grid_results$PrNoGo_grid <= 1))
  expect_equal(length(result$grid_results$PrGo_grid), length(result$grid_results$gamma_grid))
})

test_that("getgamma2cont posterior vague MM gamma values in (0, 1) or NA", {
  result <- getgamma2cont(
    nsim = 50L, prob = 'posterior', design = 'controlled',
    prior = 'vague',
    GoRegions = 1L, NoGoRegions = 9L,
    mu_t_go = c(2.0, 1.0), Sigma_t_go = .Sigma_gc,
    mu_c_go = c(0.0, 0.0), Sigma_c_go = .Sigma_gc,
    mu_t_nogo = c(2.0, 1.0), Sigma_t_nogo = .Sigma_gc,
    mu_c_nogo = c(0.0, 0.0), Sigma_c_nogo = .Sigma_gc,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 10L, n_c = 10L,
    theta_TV1 = 1.5, theta_MAV1 = 0.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    m_t = NULL, m_c = NULL,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
    r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
    nMC = NULL, CalcMethod = 'MM',
    gamma_go_grid   = seq(0.05, 0.95, by = 0.05),
    gamma_nogo_grid = seq(0.05, 0.95, by = 0.05),
    seed = 1L
  )
  if (!is.na(result$gamma_go))   expect_true(result$gamma_go   > 0 && result$gamma_go   < 1)
  if (!is.na(result$gamma_nogo)) expect_true(result$gamma_nogo > 0 && result$gamma_nogo < 1)
})

test_that("getgamma2cont posterior uncontrolled vague MM returns correct class", {
  result <- getgamma2cont(
    nsim = 50L, prob = 'posterior', design = 'uncontrolled',
    prior = 'vague',
    GoRegions = 1L, NoGoRegions = 9L,
    mu_t_go = c(2.0, 1.0), Sigma_t_go = .Sigma_gc,
    mu_c_go = NULL, Sigma_c_go = NULL,
    mu_t_nogo = c(2.0, 1.0), Sigma_t_nogo = .Sigma_gc,
    mu_c_nogo = NULL, Sigma_c_nogo = NULL,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 10L, n_c = NULL,
    theta_TV1 = 1.5, theta_MAV1 = 0.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    m_t = NULL, m_c = NULL,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = c(0.0, 0.0), Lambda0_c = NULL,
    r = 1.0,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
    nMC = NULL, CalcMethod = 'MM',
    gamma_go_grid   = seq(0.05, 0.95, by = 0.05),
    gamma_nogo_grid = seq(0.05, 0.95, by = 0.05),
    seed = 2L
  )
  expect_s3_class(result, "getgamma2cont")
})

test_that("getgamma2cont predictive vague MC returns correct class", {
  result <- getgamma2cont(
    nsim = 30L, prob = 'predictive', design = 'controlled',
    prior = 'vague',
    GoRegions = 1L, NoGoRegions = 4L,
    mu_t_go = c(2.0, 1.0), Sigma_t_go = .Sigma_gc,
    mu_c_go = c(0.0, 0.0), Sigma_c_go = .Sigma_gc,
    mu_t_nogo = c(2.0, 1.0), Sigma_t_nogo = .Sigma_gc,
    mu_c_nogo = c(0.0, 0.0), Sigma_c_nogo = .Sigma_gc,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 10L, n_c = 10L,
    theta_TV1 = NULL, theta_MAV1 = NULL,
    theta_TV2 = NULL, theta_MAV2 = NULL,
    theta_NULL1 = 1.0, theta_NULL2 = 0.5,
    m_t = 20L, m_c = 20L,
    kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
    kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
    r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
    nMC = 50L, CalcMethod = 'MC',
    gamma_go_grid   = seq(0.05, 0.95, by = 0.05),
    gamma_nogo_grid = seq(0.05, 0.95, by = 0.05),
    seed = 3L
  )
  expect_s3_class(result, "getgamma2cont")
  if (!is.na(result$PrGo_opt))
    expect_true(result$PrGo_opt >= 0 && result$PrGo_opt <= 1)
  if (!is.na(result$PrNoGo_opt))
    expect_true(result$PrNoGo_opt >= 0 && result$PrNoGo_opt <= 1)
})

test_that("getgamma2cont input validation: invalid nsim", {
  expect_error(getgamma2cont(
    nsim = 0L, prob = 'posterior', design = 'controlled',
    prior = 'vague',
    GoRegions = 1L, NoGoRegions = 9L,
    mu_t_go = c(2.0, 1.0), Sigma_t_go = .Sigma_gc,
    mu_c_go = c(0.0, 0.0), Sigma_c_go = .Sigma_gc,
    mu_t_nogo = c(2.0, 1.0), Sigma_t_nogo = .Sigma_gc,
    mu_c_nogo = c(0.0, 0.0), Sigma_c_nogo = .Sigma_gc,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 10L, n_c = 10L,
    theta_TV1 = 1.5, theta_MAV1 = 0.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    seed = 1L
  ))
})

test_that("getgamma2cont input validation: overlapping GoRegions and NoGoRegions", {
  expect_error(getgamma2cont(
    nsim = 50L, prob = 'posterior', design = 'controlled',
    prior = 'vague',
    GoRegions = c(1L, 2L), NoGoRegions = c(2L, 9L),
    mu_t_go = c(2.0, 1.0), Sigma_t_go = .Sigma_gc,
    mu_c_go = c(0.0, 0.0), Sigma_c_go = .Sigma_gc,
    mu_t_nogo = c(2.0, 1.0), Sigma_t_nogo = .Sigma_gc,
    mu_c_nogo = c(0.0, 0.0), Sigma_c_nogo = .Sigma_gc,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 10L, n_c = 10L,
    theta_TV1 = 1.5, theta_MAV1 = 0.5,
    theta_TV2 = 1.0, theta_MAV2 = 0.3,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    seed = 1L
  ))
})
