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
    theta.TV1 = 0.15, theta.MAV1 = 0.05,
    theta.TV2 = 0.15, theta.MAV2 = 0.05,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    x1_00 = 3L, x1_01 = 2L, x1_10 = 3L, x1_11 = 4L,
    x2_00 = 5L, x2_01 = 2L, x2_10 = 2L, x2_11 = 3L,
    a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
    a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
    m1 = NULL, m2 = NULL,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
    xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
    ae1 = NULL, ae2 = NULL, nMC = 50L
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
    theta.TV1 = NULL, theta.MAV1 = NULL,
    theta.TV2 = NULL, theta.MAV2 = NULL,
    theta.NULL1 = 0.15, theta.NULL2 = 0.15,
    x1_00 = 3L, x1_01 = 2L, x1_10 = 3L, x1_11 = 4L,
    x2_00 = 5L, x2_01 = 2L, x2_10 = 2L, x2_11 = 3L,
    a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
    a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
    m1 = 30L, m2 = 30L,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
    xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
    ae1 = NULL, ae2 = NULL, nMC = 50L
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
    theta.TV1 = 0.15, theta.MAV1 = 0.05,
    theta.TV2 = 0.15, theta.MAV2 = 0.05,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    x1_00 = 3L, x1_01 = 2L, x1_10 = 3L, x1_11 = 4L,
    x2_00 = NULL, x2_01 = NULL, x2_10 = NULL, x2_11 = NULL,
    a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
    a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
    m1 = NULL, m2 = NULL,
    z00 = 2L, z01 = 1L, z10 = 2L, z11 = 1L,
    xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
    xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
    ae1 = NULL, ae2 = NULL, nMC = 50L
  )
  expect_type(result, "double")
  expect_length(result, 9L)
  expect_equal(sum(result), 1, tolerance = 1e-6)
})

test_that("pbayespostpred2bin input validation works", {
  # TV > MAV violated for posterior
  expect_error(pbayespostpred2bin(
    prob = 'posterior', design = 'controlled',
    theta.TV1 = 0.05, theta.MAV1 = 0.15,
    theta.TV2 = 0.15, theta.MAV2 = 0.05,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    x1_00 = 3L, x1_01 = 2L, x1_10 = 3L, x1_11 = 4L,
    x2_00 = 5L, x2_01 = 2L, x2_10 = 2L, x2_11 = 3L,
    a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
    a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
    m1 = NULL, m2 = NULL,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
    xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
    ae1 = NULL, ae2 = NULL, nMC = 50L
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
    gamma1 = 0.60, gamma2 = 0.60,
    pi_t1 = c(0.30, 0.40), pi_t2 = c(0.35, 0.45),
    rho_t = rep(0.0, 2),
    pi_c1 = rep(0.15, 2), pi_c2 = rep(0.20, 2),
    rho_c = rep(0.0, 2),
    n1 = 10L, n2 = 10L,
    a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
    a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
    m1 = NULL, m2 = NULL,
    theta.TV1 = 0.10, theta.MAV1 = 0.05,
    theta.TV2 = 0.10, theta.MAV2 = 0.05,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
    xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
    ae1 = NULL, ae2 = NULL, nMC = 50L,
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
    gamma1 = 0.60, gamma2 = 0.60,
    pi_t1 = c(0.30, 0.40), pi_t2 = c(0.35, 0.45),
    rho_t = rep(0.0, 2),
    pi_c1 = rep(0.15, 2), pi_c2 = rep(0.20, 2),
    rho_c = rep(0.0, 2),
    n1 = 10L, n2 = 10L,
    a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
    a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
    m1 = 30L, m2 = 30L,
    theta.TV1 = NULL, theta.MAV1 = NULL,
    theta.TV2 = NULL, theta.MAV2 = NULL,
    theta.NULL1 = 0.10, theta.NULL2 = 0.10,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
    xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
    ae1 = NULL, ae2 = NULL, nMC = 50L,
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
    gamma1 = 0.60, gamma2 = 0.60,
    pi_t1 = 0.30, pi_t2 = 0.35, rho_t = 0.0,
    pi_c1 = 0.15, pi_c2 = 0.20, rho_c = 0.0,
    n1 = 10L, n2 = 10L,
    a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
    a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
    m1 = NULL, m2 = NULL,
    theta.TV1 = 0.10, theta.MAV1 = 0.05,
    theta.TV2 = 0.10, theta.MAV2 = 0.05,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
    xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
    ae1 = NULL, ae2 = NULL, nMC = 50L,
    error_if_Miss = TRUE, Gray_inc_Miss = FALSE
  ))
})

# ---------------------------------------------------------------------------
# getgamma2bin
# Note: n1 = n2 = 4L keeps allmultinom() at C(7,3) = 35 rows per arm,
#       giving 35 x 35 = 1225 combinations vs 165^2 = 27225 for n=8,
#       which is necessary to keep test runtime under 10s for CRAN.
# ---------------------------------------------------------------------------

test_that("getgamma2bin posterior controlled returns correct class and structure", {
  result <- getgamma2bin(
    prob = 'posterior', design = 'controlled',
    GoRegions = 1L, NoGoRegions = 9L,
    pi_t1 = 0.30, pi_t2 = 0.25, rho_t = 0.0,
    pi_c1 = 0.15, pi_c2 = 0.15, rho_c = 0.0,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 4L, n2 = 4L,
    a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
    a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
    theta.TV1   = 0.15, theta.MAV1  = 0.05,
    theta.TV2   = 0.10, theta.MAV2  = 0.03,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    m1 = NULL, m2 = NULL,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
    xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
    ae1 = NULL, ae2 = NULL,
    nMC = 50L, method = 'Exact',
    nsim = 100L,
    gamma1_grid = seq(0.05, 0.95, by = 0.05),
    gamma2_grid = seq(0.05, 0.95, by = 0.05)
  )
  expect_s3_class(result, "getgamma2bin")
  expect_true(all(c("gamma1", "gamma2", "PrGo_at_gamma", "PrNoGo_at_gamma",
                    "gamma1_grid", "gamma2_grid",
                    "PrGo_grid", "PrNoGo_grid") %in% names(result)))
})

test_that("getgamma2bin posterior controlled PrGo_grid and PrNoGo_grid in [0, 1]", {
  result <- getgamma2bin(
    prob = 'posterior', design = 'controlled',
    GoRegions = 1L, NoGoRegions = 9L,
    pi_t1 = 0.30, pi_t2 = 0.25, rho_t = 0.0,
    pi_c1 = 0.15, pi_c2 = 0.15, rho_c = 0.0,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 4L, n2 = 4L,
    a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
    a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
    theta.TV1   = 0.15, theta.MAV1  = 0.05,
    theta.TV2   = 0.10, theta.MAV2  = 0.03,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    m1 = NULL, m2 = NULL,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
    xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
    ae1 = NULL, ae2 = NULL,
    nMC = 50L, method = 'Exact',
    nsim = 100L,
    gamma1_grid = seq(0.05, 0.95, by = 0.05),
    gamma2_grid = seq(0.05, 0.95, by = 0.05)
  )
  expect_true(all(result$PrGo_grid >= 0 & result$PrGo_grid <= 1))
  expect_true(all(result$PrNoGo_grid >= 0 & result$PrNoGo_grid <= 1))
  expect_equal(dim(result$PrGo_grid),
               c(length(result$gamma1_grid), length(result$gamma2_grid)))
})

test_that("getgamma2bin posterior controlled gamma values in (0, 1) or NA", {
  result <- getgamma2bin(
    prob = 'posterior', design = 'controlled',
    GoRegions = 1L, NoGoRegions = 9L,
    pi_t1 = 0.30, pi_t2 = 0.25, rho_t = 0.0,
    pi_c1 = 0.15, pi_c2 = 0.15, rho_c = 0.0,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 4L, n2 = 4L,
    a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
    a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
    theta.TV1   = 0.15, theta.MAV1  = 0.05,
    theta.TV2   = 0.10, theta.MAV2  = 0.03,
    theta.NULL1 = NULL, theta.NULL2 = NULL,
    m1 = NULL, m2 = NULL,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
    xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
    ae1 = NULL, ae2 = NULL,
    nMC = 50L, method = 'Exact',
    nsim = 100L,
    gamma1_grid = seq(0.05, 0.95, by = 0.05),
    gamma2_grid = seq(0.05, 0.95, by = 0.05)
  )
  if (!is.na(result$gamma1)) expect_true(result$gamma1 > 0 && result$gamma1 < 1)
  if (!is.na(result$gamma2)) expect_true(result$gamma2 > 0 && result$gamma2 < 1)
})

test_that("getgamma2bin predictive controlled returns correct class", {
  result <- getgamma2bin(
    prob = 'predictive', design = 'controlled',
    GoRegions = 1L, NoGoRegions = 4L,
    pi_t1 = 0.30, pi_t2 = 0.25, rho_t = 0.0,
    pi_c1 = 0.15, pi_c2 = 0.15, rho_c = 0.0,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 4L, n2 = 4L,
    a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
    a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
    theta.TV1   = NULL, theta.MAV1  = NULL,
    theta.TV2   = NULL, theta.MAV2  = NULL,
    theta.NULL1 = 0.10, theta.NULL2 = 0.10,
    m1 = 20L, m2 = 20L,
    z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
    xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
    xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
    ae1 = NULL, ae2 = NULL,
    nMC = 50L, method = 'Exact',
    nsim = 100L,
    gamma1_grid = seq(0.05, 0.95, by = 0.05),
    gamma2_grid = seq(0.05, 0.95, by = 0.05)
  )
  expect_s3_class(result, "getgamma2bin")
})

test_that("getgamma2bin input validation: invalid prob", {
  expect_error(getgamma2bin(
    prob = 'invalid', design = 'controlled',
    GoRegions = 1L, NoGoRegions = 9L,
    pi_t1 = 0.30, pi_t2 = 0.25, rho_t = 0.0,
    pi_c1 = 0.15, pi_c2 = 0.15, rho_c = 0.0,
    target_go = 0.05, target_nogo = 0.20,
    n1 = 4L, n2 = 4L,
    a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
    a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
    theta.TV1 = 0.15, theta.MAV1 = 0.05,
    theta.TV2 = 0.10, theta.MAV2 = 0.03
  ))
})
