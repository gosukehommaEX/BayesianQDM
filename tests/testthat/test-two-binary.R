# Tests for two binary endpoint functions
# Note: pbayesdecisionprob2bin is computationally intensive.
#       nsim is not applicable (uses exact enumeration); nMC is minimised.
#       skip_on_cran() is applied to the most expensive tests.

# ---------------------------------------------------------------------------
# pbayespostpred2bin
# ---------------------------------------------------------------------------

test_that("pbayespostpred2bin posterior controlled returns 9 named probs", {
  result <- pbayespostpred2bin(
    prob = 'posterior', design = 'controlled',
    theta_TV1 = 0.15, theta_MAV1 = 0.05,
    theta_TV2 = 0.15, theta_MAV2 = 0.05,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    x_t_00 = 3L, x_t_01 = 2L, x_t_10 = 3L, x_t_11 = 4L,
    x_c_00 = 5L, x_c_01 = 2L, x_c_10 = 2L, x_c_11 = 3L,
    a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
    a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
    m_t = NULL, m_c = NULL,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
    xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
    alpha0e_t = NULL, alpha0e_c = NULL, nMC = 50L
  )
  expect_type(result, "double")
  expect_length(result, 9L)
  expect_named(result, paste0("R", 1:9))
  expect_true(all(result >= 0))
  expect_equal(sum(result), 1, tolerance = 1e-6)
})

test_that("pbayespostpred2bin predictive controlled returns 4 named probs", {
  set.seed(1)
  result <- pbayespostpred2bin(
    prob = 'predictive', design = 'controlled',
    theta_TV1 = NULL, theta_MAV1 = NULL,
    theta_TV2 = NULL, theta_MAV2 = NULL,
    theta_NULL1 = 0.15, theta_NULL2 = 0.15,
    x_t_00 = 3L, x_t_01 = 2L, x_t_10 = 3L, x_t_11 = 4L,
    x_c_00 = 5L, x_c_01 = 2L, x_c_10 = 2L, x_c_11 = 3L,
    a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
    a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
    m_t = 30L, m_c = 30L,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
    xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
    alpha0e_t = NULL, alpha0e_c = NULL, nMC = 50L
  )
  expect_type(result, "double")
  expect_length(result, 4L)
  expect_named(result, paste0("R", 1:4))
  expect_true(all(result >= 0))
  expect_equal(sum(result), 1, tolerance = 1e-6)
})

test_that("pbayespostpred2bin posterior uncontrolled returns 9 named probs", {
  result <- pbayespostpred2bin(
    prob = 'posterior', design = 'uncontrolled',
    theta_TV1 = 0.15, theta_MAV1 = 0.05,
    theta_TV2 = 0.15, theta_MAV2 = 0.05,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    x_t_00 = 3L, x_t_01 = 2L, x_t_10 = 3L, x_t_11 = 4L,
    x_c_00 = NULL, x_c_01 = NULL, x_c_10 = NULL, x_c_11 = NULL,
    a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
    a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
    m_t = NULL, m_c = NULL,
    z00 = 2L, z01 = 1L, z10 = 2L, z11 = 1L,
    xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
    xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
    alpha0e_t = NULL, alpha0e_c = NULL, nMC = 50L
  )
  expect_type(result, "double")
  expect_length(result, 9L)
  expect_equal(sum(result), 1, tolerance = 1e-6)
})

test_that("pbayespostpred2bin input validation works", {
  # TV > MAV violated for posterior
  expect_error(pbayespostpred2bin(
    prob = 'posterior', design = 'controlled',
    theta_TV1 = 0.05, theta_MAV1 = 0.15,
    theta_TV2 = 0.15, theta_MAV2 = 0.05,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    x_t_00 = 3L, x_t_01 = 2L, x_t_10 = 3L, x_t_11 = 4L,
    x_c_00 = 5L, x_c_01 = 2L, x_c_10 = 2L, x_c_11 = 3L,
    a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
    a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
    m_t = NULL, m_c = NULL,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
    xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
    alpha0e_t = NULL, alpha0e_c = NULL, nMC = 50L
  ))
})

# ---------------------------------------------------------------------------
# pbayesdecisionprob2bin
# ---------------------------------------------------------------------------

test_that("pbayesdecisionprob2bin posterior controlled returns data.frame", {
  skip_on_cran()
  result <- pbayesdecisionprob2bin(
    prob = 'posterior', design = 'controlled',
    GoRegions = 1L, NoGoRegions = 9L,
    gamma_go = 0.60, gamma_nogo = 0.60,
    pi_t1 = c(0.30, 0.40), pi_t2 = c(0.35, 0.45),
    rho_t = rep(0.0, 2),
    pi_c1 = rep(0.15, 2), pi_c2 = rep(0.20, 2),
    rho_c = rep(0.0, 2),
    n_t = 10L, n_c = 10L,
    a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
    a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
    m_t = NULL, m_c = NULL,
    theta_TV1 = 0.10, theta_MAV1 = 0.05,
    theta_TV2 = 0.10, theta_MAV2 = 0.05,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
    xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
    alpha0e_t = NULL, alpha0e_c = NULL, nMC = 50L,
    error_if_Miss = TRUE, Gray_inc_Miss = FALSE
  )
  expect_s3_class(result, "pbayesdecisionprob2bin")
  df <- as.data.frame(result)
  expect_true(all(c("pi_t1", "pi_t2", "Go", "Gray", "NoGo") %in% names(df)))
  expect_equal(nrow(df), 2L)
  expect_true(all(df$Go  >= 0 & df$Go  <= 1))
  expect_true(all(df$NoGo >= 0 & df$NoGo <= 1))
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

test_that("pbayesdecisionprob2bin predictive controlled returns data.frame", {
  skip_on_cran()
  result <- pbayesdecisionprob2bin(
    prob = 'predictive', design = 'controlled',
    GoRegions = 1L, NoGoRegions = 4L,
    gamma_go = 0.60, gamma_nogo = 0.60,
    pi_t1 = c(0.30, 0.40), pi_t2 = c(0.35, 0.45),
    rho_t = rep(0.0, 2),
    pi_c1 = rep(0.15, 2), pi_c2 = rep(0.20, 2),
    rho_c = rep(0.0, 2),
    n_t = 10L, n_c = 10L,
    a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
    a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
    m_t = 30L, m_c = 30L,
    theta_TV1 = NULL, theta_MAV1 = NULL,
    theta_TV2 = NULL, theta_MAV2 = NULL,
    theta_NULL1 = 0.10, theta_NULL2 = 0.10,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
    xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
    alpha0e_t = NULL, alpha0e_c = NULL, nMC = 50L,
    error_if_Miss = TRUE, Gray_inc_Miss = FALSE
  )
  expect_s3_class(result, "pbayesdecisionprob2bin")
  df <- as.data.frame(result)
  expect_equal(nrow(df), 2L)
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

test_that("pbayesdecisionprob2bin input validation works", {
  # GoRegions and NoGoRegions overlap
  expect_error(pbayesdecisionprob2bin(
    prob = 'posterior', design = 'controlled',
    GoRegions = c(1L, 2L), NoGoRegions = c(2L, 9L),
    gamma_go = 0.60, gamma_nogo = 0.60,
    pi_t1 = 0.30, pi_t2 = 0.35, rho_t = 0.0,
    pi_c1 = 0.15, pi_c2 = 0.20, rho_c = 0.0,
    n_t = 10L, n_c = 10L,
    a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
    a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
    m_t = NULL, m_c = NULL,
    theta_TV1 = 0.10, theta_MAV1 = 0.05,
    theta_TV2 = 0.10, theta_MAV2 = 0.05,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
    xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
    alpha0e_t = NULL, alpha0e_c = NULL, nMC = 50L,
    error_if_Miss = TRUE, Gray_inc_Miss = FALSE
  ))
})

test_that("pbayesdecisionprob2bin posterior uncontrolled returns correct class", {
  skip_on_cran()
  result <- pbayesdecisionprob2bin(
    prob = 'posterior', design = 'uncontrolled',
    GoRegions = 1L, NoGoRegions = 9L,
    gamma_go = 0.60, gamma_nogo = 0.60,
    pi_t1 = c(0.30, 0.40), pi_t2 = c(0.35, 0.45),
    rho_t = rep(0.0, 2),
    pi_c1 = NULL, pi_c2 = NULL, rho_c = NULL,
    n_t = 10L, n_c = NULL,
    a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
    a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
    m_t = NULL, m_c = NULL,
    theta_TV1 = 0.10, theta_MAV1 = 0.05,
    theta_TV2 = 0.10, theta_MAV2 = 0.05,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    z00 = 2L, z01 = 1L, z10 = 2L, z11 = 1L,
    xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
    xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
    alpha0e_t = NULL, alpha0e_c = NULL, nMC = 50L,
    error_if_Miss = FALSE, Gray_inc_Miss = FALSE
  )
  expect_s3_class(result, "pbayesdecisionprob2bin")
  df <- as.data.frame(result)
  expect_true(all(c("pi_t1", "pi_t2", "Go", "Gray", "NoGo") %in% names(df)))
  expect_false(any(c("pi_c1", "pi_c2") %in% names(df)))
  expect_equal(nrow(df), 2L)
  expect_true(all(df$Go  >= 0 & df$Go  <= 1))
  expect_true(all(df$NoGo >= 0 & df$NoGo <= 1))
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

# ---------------------------------------------------------------------------
# getgamma2bin
# Note: n_t = n_c = 4L keeps allmultinom() at C(7,3) = 35 rows per arm,
#       giving 35 x 35 = 1225 combinations vs 165^2 = 27225 for n=8,
#       which is necessary to keep test runtime under 10s for CRAN.
# ---------------------------------------------------------------------------

test_that("getgamma2bin posterior controlled returns correct class and structure", {
  result <- getgamma2bin(
    prob = 'posterior', design = 'controlled',
    GoRegions = 1L, NoGoRegions = 9L,
    pi_t1_go = 0.30, pi_t2_go = 0.25, rho_t_go = 0.0,
    pi_c1_go = 0.15, pi_c2_go = 0.15, rho_c_go = 0.0,
    pi_t1_nogo = 0.30, pi_t2_nogo = 0.25, rho_t_nogo = 0.0,
    pi_c1_nogo = 0.15, pi_c2_nogo = 0.15, rho_c_nogo = 0.0,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 4L, n_c = 4L,
    a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
    a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
    theta_TV1   = 0.15, theta_MAV1  = 0.05,
    theta_TV2   = 0.10, theta_MAV2  = 0.03,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    m_t = NULL, m_c = NULL,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
    xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
    alpha0e_t = NULL, alpha0e_c = NULL,
    nMC = 50L,
    gamma_go_grid   = seq(0.05, 0.95, by = 0.05),
    gamma_nogo_grid = seq(0.05, 0.95, by = 0.05)
  )
  expect_s3_class(result, "getgamma2bin")
  expect_true(all(c("gamma_go", "gamma_nogo", "PrGo_opt", "PrNoGo_opt",
                    "target_go", "target_nogo", "grid_results") %in% names(result)))
  expect_true(all(c("gamma_grid", "PrGo_grid", "PrNoGo_grid") %in%
                    names(result$grid_results)))
})

test_that("getgamma2bin posterior controlled PrGo_grid and PrNoGo_grid in [0, 1]", {
  result <- getgamma2bin(
    prob = 'posterior', design = 'controlled',
    GoRegions = 1L, NoGoRegions = 9L,
    pi_t1_go = 0.30, pi_t2_go = 0.25, rho_t_go = 0.0,
    pi_c1_go = 0.15, pi_c2_go = 0.15, rho_c_go = 0.0,
    pi_t1_nogo = 0.30, pi_t2_nogo = 0.25, rho_t_nogo = 0.0,
    pi_c1_nogo = 0.15, pi_c2_nogo = 0.15, rho_c_nogo = 0.0,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 4L, n_c = 4L,
    a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
    a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
    theta_TV1   = 0.15, theta_MAV1  = 0.05,
    theta_TV2   = 0.10, theta_MAV2  = 0.03,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    m_t = NULL, m_c = NULL,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
    xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
    alpha0e_t = NULL, alpha0e_c = NULL,
    nMC = 50L,
    gamma_go_grid   = seq(0.05, 0.95, by = 0.05),
    gamma_nogo_grid = seq(0.05, 0.95, by = 0.05)
  )
  expect_true(all(result$grid_results$PrGo_grid   >= 0 & result$grid_results$PrGo_grid   <= 1))
  expect_true(all(result$grid_results$PrNoGo_grid >= 0 & result$grid_results$PrNoGo_grid <= 1))
  expect_equal(length(result$grid_results$PrGo_grid), length(result$grid_results$gamma_grid))
})

test_that("getgamma2bin posterior controlled gamma values in (0, 1) or NA", {
  result <- getgamma2bin(
    prob = 'posterior', design = 'controlled',
    GoRegions = 1L, NoGoRegions = 9L,
    pi_t1_go = 0.30, pi_t2_go = 0.25, rho_t_go = 0.0,
    pi_c1_go = 0.15, pi_c2_go = 0.15, rho_c_go = 0.0,
    pi_t1_nogo = 0.30, pi_t2_nogo = 0.25, rho_t_nogo = 0.0,
    pi_c1_nogo = 0.15, pi_c2_nogo = 0.15, rho_c_nogo = 0.0,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 4L, n_c = 4L,
    a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
    a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
    theta_TV1   = 0.15, theta_MAV1  = 0.05,
    theta_TV2   = 0.10, theta_MAV2  = 0.03,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    m_t = NULL, m_c = NULL,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
    xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
    alpha0e_t = NULL, alpha0e_c = NULL,
    nMC = 50L,
    gamma_go_grid   = seq(0.05, 0.95, by = 0.05),
    gamma_nogo_grid = seq(0.05, 0.95, by = 0.05)
  )
  if (!is.na(result$gamma_go))   expect_true(result$gamma_go   > 0 && result$gamma_go   < 1)
  if (!is.na(result$gamma_nogo)) expect_true(result$gamma_nogo > 0 && result$gamma_nogo < 1)
})

test_that("getgamma2bin predictive controlled returns correct class", {
  result <- getgamma2bin(
    prob = 'predictive', design = 'controlled',
    GoRegions = 1L, NoGoRegions = 4L,
    pi_t1_go = 0.30, pi_t2_go = 0.25, rho_t_go = 0.0,
    pi_c1_go = 0.15, pi_c2_go = 0.15, rho_c_go = 0.0,
    pi_t1_nogo = 0.30, pi_t2_nogo = 0.25, rho_t_nogo = 0.0,
    pi_c1_nogo = 0.15, pi_c2_nogo = 0.15, rho_c_nogo = 0.0,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 4L, n_c = 4L,
    a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
    a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
    theta_TV1   = NULL, theta_MAV1  = NULL,
    theta_TV2   = NULL, theta_MAV2  = NULL,
    theta_NULL1 = 0.10, theta_NULL2 = 0.10,
    m_t = 20L, m_c = 20L,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
    xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
    alpha0e_t = NULL, alpha0e_c = NULL,
    nMC = 50L,
    gamma_go_grid   = seq(0.05, 0.95, by = 0.05),
    gamma_nogo_grid = seq(0.05, 0.95, by = 0.05)
  )
  expect_s3_class(result, "getgamma2bin")
})

test_that("getgamma2bin posterior uncontrolled returns correct class", {
  result <- getgamma2bin(
    prob = 'posterior', design = 'uncontrolled',
    GoRegions = 1L, NoGoRegions = 9L,
    pi_t1_go = 0.30, pi_t2_go = 0.25, rho_t_go = 0.0,
    pi_c1_go = NULL, pi_c2_go = NULL, rho_c_go = NULL,
    pi_t1_nogo = 0.30, pi_t2_nogo = 0.25, rho_t_nogo = 0.0,
    pi_c1_nogo = NULL, pi_c2_nogo = NULL, rho_c_nogo = NULL,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 4L, n_c = 4L,
    a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
    a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
    theta_TV1   = 0.15, theta_MAV1  = 0.05,
    theta_TV2   = 0.10, theta_MAV2  = 0.03,
    theta_NULL1 = NULL, theta_NULL2 = NULL,
    m_t = NULL, m_c = NULL,
    z00 = 1L, z01 = 1L, z10 = 1L, z11 = 1L,
    xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
    xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
    alpha0e_t = NULL, alpha0e_c = NULL,
    nMC = 50L,
    gamma_go_grid   = seq(0.05, 0.95, by = 0.05),
    gamma_nogo_grid = seq(0.05, 0.95, by = 0.05)
  )
  expect_s3_class(result, "getgamma2bin")
  expect_true(all(c("gamma_go", "gamma_nogo", "PrGo_opt", "PrNoGo_opt",
                    "target_go", "target_nogo", "grid_results") %in% names(result)))
  expect_true(all(result$grid_results$PrGo_grid   >= 0 & result$grid_results$PrGo_grid   <= 1))
  expect_true(all(result$grid_results$PrNoGo_grid >= 0 & result$grid_results$PrNoGo_grid <= 1))
})

test_that("getgamma2bin input validation: invalid prob", {
  expect_error(getgamma2bin(
    prob = 'invalid', design = 'controlled',
    GoRegions = 1L, NoGoRegions = 9L,
    pi_t1_go = 0.30, pi_t2_go = 0.25, rho_t_go = 0.0,
    pi_c1_go = 0.15, pi_c2_go = 0.15, rho_c_go = 0.0,
    pi_t1_nogo = 0.30, pi_t2_nogo = 0.25, rho_t_nogo = 0.0,
    pi_c1_nogo = 0.15, pi_c2_nogo = 0.15, rho_c_nogo = 0.0,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 4L, n_c = 4L,
    a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
    a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
    theta_TV1 = 0.15, theta_MAV1 = 0.05,
    theta_TV2 = 0.10, theta_MAV2 = 0.03
  ))
})
