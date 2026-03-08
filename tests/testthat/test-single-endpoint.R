# Tests for single-endpoint functions (binary and continuous)

# ---------------------------------------------------------------------------
# pbayespostpred1bin
# Signature: pbayespostpred1bin(prob, design, theta0, n_t, n_c, y_t, y_c,
#            a_t, a_c, b_t, b_c, m_t, m_c, z, ne_t, ne_c, ye_t, ye_c, alpha0e_t, alpha0e_c,
#            lower.tail)
# Note: no nMC argument; y_t and y_c must have the same length
# ---------------------------------------------------------------------------

test_that("pbayespostpred1bin posterior controlled returns scalar in [0, 1]", {
  result <- pbayespostpred1bin(
    prob = 'posterior', design = 'controlled', theta0 = 0.15,
    n_t = 20, n_c = 20, y_t = 12, y_c = 8,
    a_t = 0.5, b_t = 0.5, a_c = 0.5, b_c = 0.5,
    m_t = NULL, m_c = NULL, z = NULL,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    lower.tail = FALSE
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1bin predictive controlled returns scalar in [0, 1]", {
  result <- pbayespostpred1bin(
    prob = 'predictive', design = 'controlled', theta0 = 0.15,
    n_t = 20, n_c = 20, y_t = 12, y_c = 8,
    a_t = 0.5, b_t = 0.5, a_c = 0.5, b_c = 0.5,
    m_t = 30, m_c = 30, z = NULL,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    lower.tail = FALSE
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1bin posterior uncontrolled returns scalar in [0, 1]", {
  result <- pbayespostpred1bin(
    prob = 'posterior', design = 'uncontrolled', theta0 = 0.15,
    n_t = 20, n_c = 20, y_t = 12, y_c = NULL,
    a_t = 0.5, b_t = 0.5, a_c = 0.5, b_c = 0.5,
    m_t = NULL, m_c = NULL, z = 5L,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    lower.tail = FALSE
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1bin posterior external returns scalar in [0, 1]", {
  result <- pbayespostpred1bin(
    prob = 'posterior', design = 'external', theta0 = 0.15,
    n_t = 20, n_c = 20, y_t = 12, y_c = 8,
    a_t = 0.5, b_t = 0.5, a_c = 0.5, b_c = 0.5,
    m_t = NULL, m_c = NULL, z = NULL,
    ne_t = 30, ne_c = 30, ye_t = 10, ye_c = 6, alpha0e_t = 0.5, alpha0e_c = 0.5,
    lower.tail = FALSE
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1bin vectorised y_t and y_c returns correct length", {
  result <- pbayespostpred1bin(
    prob = 'posterior', design = 'controlled', theta0 = 0.15,
    n_t = 20, n_c = 20, y_t = 8:12, y_c = rep(6, 5),
    a_t = 0.5, b_t = 0.5, a_c = 0.5, b_c = 0.5,
    m_t = NULL, m_c = NULL, z = NULL,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    lower.tail = FALSE
  )
  expect_length(result, 5L)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("pbayespostpred1bin input validation works", {
  # y_t > n_t
  expect_error(pbayespostpred1bin(
    prob = 'posterior', design = 'controlled', theta0 = 0.15,
    n_t = 20, n_c = 20, y_t = 25, y_c = 8,
    a_t = 0.5, b_t = 0.5, a_c = 0.5, b_c = 0.5,
    m_t = NULL, m_c = NULL, z = NULL,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL
  ))
  # m_t missing for predictive
  expect_error(pbayespostpred1bin(
    prob = 'predictive', design = 'controlled', theta0 = 0.15,
    n_t = 20, n_c = 20, y_t = 12, y_c = 8,
    a_t = 0.5, b_t = 0.5, a_c = 0.5, b_c = 0.5,
    m_t = NULL, m_c = NULL, z = NULL,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL
  ))
})

# ---------------------------------------------------------------------------
# pbayespostpred1cont
# ---------------------------------------------------------------------------

test_that("pbayespostpred1cont posterior vague NI returns scalar in [0, 1]", {
  result <- pbayespostpred1cont(
    prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
    theta0 = 1, nMC = NULL,
    n_t = 15, n_c = 15, m_t = NULL, m_c = NULL,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
    bar_y_t = 3, s_t = 1.5, bar_y_c = 1, s_c = 1.2, r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1cont posterior vague MM returns scalar in [0, 1]", {
  result <- pbayespostpred1cont(
    prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'MM',
    theta0 = 1, nMC = NULL,
    n_t = 15, n_c = 15, m_t = NULL, m_c = NULL,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
    bar_y_t = 3, s_t = 1.5, bar_y_c = 1, s_c = 1.2, r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1cont NI and MM agree closely", {
  p_ni <- pbayespostpred1cont(
    prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
    theta0 = 1, nMC = NULL,
    n_t = 15, n_c = 15, m_t = NULL, m_c = NULL,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
    bar_y_t = 3, s_t = 1.5, bar_y_c = 1, s_c = 1.2, r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL
  )
  p_mm <- pbayespostpred1cont(
    prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'MM',
    theta0 = 1, nMC = NULL,
    n_t = 15, n_c = 15, m_t = NULL, m_c = NULL,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
    bar_y_t = 3, s_t = 1.5, bar_y_c = 1, s_c = 1.2, r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL
  )
  expect_equal(p_ni, p_mm, tolerance = 0.02)
})

test_that("pbayespostpred1cont predictive NI returns scalar in [0, 1]", {
  result <- pbayespostpred1cont(
    prob = 'predictive', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
    theta0 = 1, nMC = NULL,
    n_t = 15, n_c = 15, m_t = 30, m_c = 30,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
    bar_y_t = 3, s_t = 1.5, bar_y_c = 1, s_c = 1.2, r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1cont N-Inv-Chisq prior returns scalar in [0, 1]", {
  result <- pbayespostpred1cont(
    prob = 'posterior', design = 'controlled', prior = 'N-Inv-Chisq',
    CalcMethod = 'NI', theta0 = 1, nMC = NULL,
    n_t = 15, n_c = 15, m_t = NULL, m_c = NULL,
    kappa0_t = 5, kappa0_c = 5, nu0_t = 5, nu0_c = 5,
    mu0_t = 3, mu0_c = 1, sigma0_t = 1.5, sigma0_c = 1.2,
    bar_y_t = 3, s_t = 1.5, bar_y_c = 1, s_c = 1.2, r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1cont input validation works", {
  # m_t missing for predictive
  expect_error(pbayespostpred1cont(
    prob = 'predictive', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
    theta0 = 1, nMC = NULL,
    n_t = 15, n_c = 15, m_t = NULL, m_c = NULL,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
    bar_y_t = 3, s_t = 1.5, bar_y_c = 1, s_c = 1.2, r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL
  ))
})

test_that("pbayespostpred1cont external N-Inv-Chisq posterior returns scalar in [0, 1]", {
  result <- pbayespostpred1cont(
    prob = 'posterior', design = 'external', prior = 'N-Inv-Chisq',
    CalcMethod = 'MM', theta0 = 1, nMC = NULL,
    n_t = 20, n_c = 10, m_t = NULL, m_c = NULL,
    kappa0_t = 5, kappa0_c = 5, nu0_t = 5, nu0_c = 5,
    mu0_t = 0, mu0_c = 0, sigma0_t = 30, sigma0_c = 30,
    bar_y_t = 5, s_t = 20, bar_y_c = 2, s_c = 18, r = NULL,
    ne_t = NULL, ne_c = 10L, alpha0e_t = NULL, alpha0e_c = 0.5,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = -2, se_c = 25
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1cont external N-Inv-Chisq predictive returns scalar in [0, 1]", {
  result <- pbayespostpred1cont(
    prob = 'predictive', design = 'external', prior = 'N-Inv-Chisq',
    CalcMethod = 'MM', theta0 = 0, nMC = NULL,
    n_t = 20, n_c = 10, m_t = 30, m_c = 30,
    kappa0_t = 5, kappa0_c = 5, nu0_t = 5, nu0_c = 5,
    mu0_t = 0, mu0_c = 0, sigma0_t = 30, sigma0_c = 30,
    bar_y_t = 5, s_t = 20, bar_y_c = 2, s_c = 18, r = NULL,
    ne_t = NULL, ne_c = 10L, alpha0e_t = NULL, alpha0e_c = 0.5,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = -2, se_c = 25
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1cont external N-Inv-Chisq with both arms returns scalar in [0, 1]", {
  result <- pbayespostpred1cont(
    prob = 'posterior', design = 'external', prior = 'N-Inv-Chisq',
    CalcMethod = 'NI', theta0 = 1, nMC = NULL,
    n_t = 15, n_c = 15, m_t = NULL, m_c = NULL,
    kappa0_t = 5, kappa0_c = 5, nu0_t = 5, nu0_c = 5,
    mu0_t = 3, mu0_c = 1, sigma0_t = 1.5, sigma0_c = 1.2,
    bar_y_t = 3, s_t = 1.5, bar_y_c = 1, s_c = 1.2, r = NULL,
    ne_t = 20L, ne_c = 20L, alpha0e_t = 0.5, alpha0e_c = 0.5,
    bar_ye_t = 3, se_t = 1.5, bar_ye_c = 1, se_c = 1.2
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1cont posterior uncontrolled vague MM returns scalar in [0, 1]", {
  result <- pbayespostpred1cont(
    prob = 'posterior', design = 'uncontrolled', prior = 'vague', CalcMethod = 'MM',
    theta0 = 1, nMC = NULL,
    n_t = 15, n_c = NULL, m_t = NULL, m_c = NULL,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = 1.0, sigma0_t = NULL, sigma0_c = NULL,
    bar_y_t = 3, s_t = 1.5, bar_y_c = NULL, s_c = NULL, r = 1.0,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1cont posterior uncontrolled N-Inv-Chisq NI returns scalar in [0, 1]", {
  result <- pbayespostpred1cont(
    prob = 'posterior', design = 'uncontrolled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
    theta0 = 1, nMC = NULL,
    n_t = 15, n_c = NULL, m_t = NULL, m_c = NULL,
    kappa0_t = 2, kappa0_c = NULL, nu0_t = 5, nu0_c = NULL,
    mu0_t = 3.0, mu0_c = 1.0, sigma0_t = 1.5, sigma0_c = NULL,
    bar_y_t = 3, s_t = 1.5, bar_y_c = NULL, s_c = NULL, r = 1.0,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1cont uncontrolled NI and MM agree closely", {
  p_ni <- pbayespostpred1cont(
    prob = 'posterior', design = 'uncontrolled', prior = 'vague', CalcMethod = 'NI',
    theta0 = 1, nMC = NULL,
    n_t = 15, n_c = NULL, m_t = NULL, m_c = NULL,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = 1.0, sigma0_t = NULL, sigma0_c = NULL,
    bar_y_t = 3, s_t = 1.5, bar_y_c = NULL, s_c = NULL, r = 1.0,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL
  )
  p_mm <- pbayespostpred1cont(
    prob = 'posterior', design = 'uncontrolled', prior = 'vague', CalcMethod = 'MM',
    theta0 = 1, nMC = NULL,
    n_t = 15, n_c = NULL, m_t = NULL, m_c = NULL,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = 1.0, sigma0_t = NULL, sigma0_c = NULL,
    bar_y_t = 3, s_t = 1.5, bar_y_c = NULL, s_c = NULL, r = 1.0,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL
  )
  expect_equal(p_ni, p_mm, tolerance = 0.02)
})

# ---------------------------------------------------------------------------
# pbayesdecisionprob1bin
# Note: no nMC argument
# ---------------------------------------------------------------------------

test_that("pbayesdecisionprob1bin posterior controlled returns correct class", {
  result <- pbayesdecisionprob1bin(
    prob = 'posterior', design = 'controlled',
    theta_TV = 0.2, theta_MAV = 0.05, theta_NULL = NULL,
    gamma_go = 0.8, gamma_nogo = 0.2,
    pi_t = c(0.3, 0.5), pi_c = rep(0.2, 2),
    n_t = 15, n_c = 15,
    a_t = 0.5, b_t = 0.5, a_c = 0.5, b_c = 0.5,
    z = NULL, m_t = NULL, m_c = NULL,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL
  )
  expect_s3_class(result, "pbayesdecisionprob1bin")
  df <- as.data.frame(result)
  expect_true(all(c("pi_t", "Go", "Gray", "NoGo") %in% names(df)))
  expect_equal(nrow(df), 2L)
  expect_true(all(df$Go  >= 0 & df$Go  <= 1))
  expect_true(all(df$NoGo >= 0 & df$NoGo <= 1))
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

test_that("pbayesdecisionprob1bin predictive controlled returns correct class", {
  result <- pbayesdecisionprob1bin(
    prob = 'predictive', design = 'controlled',
    theta_TV = NULL, theta_MAV = NULL, theta_NULL = 0.15,
    gamma_go = 0.8, gamma_nogo = 0.2,
    pi_t = c(0.3, 0.5), pi_c = rep(0.2, 2),
    n_t = 15, n_c = 15,
    a_t = 0.5, b_t = 0.5, a_c = 0.5, b_c = 0.5,
    z = NULL, m_t = 40, m_c = 40,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL
  )
  expect_s3_class(result, "pbayesdecisionprob1bin")
  df <- as.data.frame(result)
  expect_equal(nrow(df), 2L)
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

test_that("pbayesdecisionprob1bin input validation works", {
  # theta_TV missing for posterior
  expect_error(pbayesdecisionprob1bin(
    prob = 'posterior', design = 'controlled',
    theta_TV = NULL, theta_MAV = 0.05, theta_NULL = NULL,
    gamma_go = 0.8, gamma_nogo = 0.2,
    pi_t = 0.3, pi_c = 0.2, n_t = 15, n_c = 15,
    a_t = 0.5, b_t = 0.5, a_c = 0.5, b_c = 0.5,
    z = NULL, m_t = NULL, m_c = NULL,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL
  ))
})

test_that("pbayesdecisionprob1bin posterior uncontrolled returns correct class", {
  result <- pbayesdecisionprob1bin(
    prob = 'posterior', design = 'uncontrolled',
    theta_TV = 0.2, theta_MAV = 0.05, theta_NULL = NULL,
    gamma_go = 0.8, gamma_nogo = 0.2,
    pi_t = c(0.3, 0.5), pi_c = NULL,
    n_t = 15, n_c = 15,
    a_t = 0.5, b_t = 0.5, a_c = 0.5, b_c = 0.5,
    z = 3L, m_t = NULL, m_c = NULL,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL
  )
  expect_s3_class(result, "pbayesdecisionprob1bin")
  df <- as.data.frame(result)
  expect_true(all(c("pi_t", "Go", "Gray", "NoGo") %in% names(df)))
  expect_false("pi_c" %in% names(df))
  expect_equal(nrow(df), 2L)
  expect_true(all(df$Go  >= 0 & df$Go  <= 1))
  expect_true(all(df$NoGo >= 0 & df$NoGo <= 1))
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

# ---------------------------------------------------------------------------
# pbayesdecisionprob1cont
# ---------------------------------------------------------------------------

test_that("pbayesdecisionprob1cont posterior vague MM returns correct class", {
  set.seed(1)
  result <- pbayesdecisionprob1cont(
    nsim = 20L, prob = 'posterior', design = 'controlled',
    prior = 'vague', CalcMethod = 'MM',
    theta_TV = 1.5, theta_MAV = 0.5, theta_NULL = NULL, nMC = NULL,
    gamma_go = 0.8, gamma_nogo = 0.2,
    n_t = 15, n_c = 15, m_t = NULL, m_c = NULL,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
    mu_t = c(2, 3), mu_c = c(0, 0), sigma_t = 1.5, sigma_c = 1.2, r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL, seed = 1
  )
  expect_s3_class(result, "pbayesdecisionprob1cont")
  df <- as.data.frame(result)
  expect_true(all(c("mu_t", "mu_c", "Go", "Gray", "NoGo") %in% names(df)))
  expect_equal(nrow(df), 2L)
  expect_true(all(df$Go  >= 0 & df$Go  <= 1))
  expect_true(all(df$NoGo >= 0 & df$NoGo <= 1))
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

test_that("pbayesdecisionprob1cont predictive vague MM returns correct class", {
  set.seed(2)
  result <- pbayesdecisionprob1cont(
    nsim = 20L, prob = 'predictive', design = 'controlled',
    prior = 'vague', CalcMethod = 'MM',
    theta_TV = NULL, theta_MAV = NULL, theta_NULL = 1, nMC = NULL,
    gamma_go = 0.8, gamma_nogo = 0.2,
    n_t = 15, n_c = 15, m_t = 30, m_c = 30,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
    mu_t = c(2, 3), mu_c = c(0, 0), sigma_t = 1.5, sigma_c = 1.2, r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL, seed = 2
  )
  expect_s3_class(result, "pbayesdecisionprob1cont")
  df <- as.data.frame(result)
  expect_equal(nrow(df), 2L)
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

test_that("pbayesdecisionprob1cont input validation works", {
  # theta_TV missing for posterior
  expect_error(pbayesdecisionprob1cont(
    nsim = 10L, prob = 'posterior', design = 'controlled',
    prior = 'vague', CalcMethod = 'MM',
    theta_TV = NULL, theta_MAV = 0.5, theta_NULL = NULL, nMC = NULL,
    gamma_go = 0.8, gamma_nogo = 0.2,
    n_t = 15, n_c = 15, m_t = NULL, m_c = NULL,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
    mu_t = 2, mu_c = 0, sigma_t = 1.5, sigma_c = 1.2, r = NULL,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL, seed = 1
  ))
})


test_that("pbayesdecisionprob1cont posterior uncontrolled vague MM returns correct class", {
  result <- pbayesdecisionprob1cont(
    nsim = 20L, prob = 'posterior', design = 'uncontrolled',
    prior = 'vague', CalcMethod = 'MM',
    theta_TV = 1.5, theta_MAV = 0.5, theta_NULL = NULL, nMC = NULL,
    gamma_go = 0.8, gamma_nogo = 0.2,
    n_t = 15, n_c = NULL, m_t = NULL, m_c = NULL,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = 0.0, sigma0_t = NULL, sigma0_c = NULL,
    mu_t = c(2, 3), mu_c = NULL, sigma_t = 1.5, sigma_c = NULL, r = 1.0,
    ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL, seed = 1
  )
  expect_s3_class(result, "pbayesdecisionprob1cont")
  df <- as.data.frame(result)
  expect_true(all(c("mu_t", "Go", "Gray", "NoGo") %in% names(df)))
  expect_false("mu_c" %in% names(df))
  expect_equal(nrow(df), 2L)
  expect_true(all(df$Go  >= 0 & df$Go  <= 1))
  expect_true(all(df$NoGo >= 0 & df$NoGo <= 1))
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

# ---------------------------------------------------------------------------
# getgamma1bin
# ---------------------------------------------------------------------------

test_that("getgamma1bin posterior controlled returns correct class and structure", {
  result <- getgamma1bin(
    prob = 'posterior', design = 'controlled',
    theta_TV = 0.20, theta_MAV = 0.05, theta_NULL = NULL,
    pi_t = 0.15, pi_c = 0.15,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n_t = 12L, n_c = 12L,
    a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
    z = NULL, m_t = NULL, m_c = NULL,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL
  )
  expect_s3_class(result, "getgamma1bin")
  expect_true(all(c("gamma_go", "gamma_nogo", "PrGo_at_gamma_go", "PrNoGo_at_gamma_nogo",
                    "gamma_grid", "PrGo_grid", "PrNoGo_grid") %in% names(result)))
})

test_that("getgamma1bin posterior controlled gamma values in (0, 1) or NA", {
  result <- getgamma1bin(
    prob = 'posterior', design = 'controlled',
    theta_TV = 0.20, theta_MAV = 0.05, theta_NULL = NULL,
    pi_t = 0.15, pi_c = 0.15,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n_t = 12L, n_c = 12L,
    a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
    z = NULL, m_t = NULL, m_c = NULL,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL
  )
  if (!is.na(result$gamma_go)) expect_true(result$gamma_go > 0 && result$gamma_go < 1)
  if (!is.na(result$gamma_nogo)) expect_true(result$gamma_nogo > 0 && result$gamma_nogo < 1)
  expect_equal(length(result$PrGo_grid), length(result$gamma_grid))
  expect_equal(length(result$PrNoGo_grid), length(result$gamma_grid))
})

test_that("getgamma1bin posterior uncontrolled returns correct class", {
  result <- getgamma1bin(
    prob = 'posterior', design = 'uncontrolled',
    theta_TV = 0.20, theta_MAV = 0.05, theta_NULL = NULL,
    pi_t = 0.15, pi_c = NULL,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n_t = 12L, n_c = 12L,
    a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
    z = 3L, m_t = NULL, m_c = NULL,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL
  )
  expect_s3_class(result, "getgamma1bin")
})

test_that("getgamma1bin predictive controlled returns correct class", {
  result <- getgamma1bin(
    prob = 'predictive', design = 'controlled',
    theta_TV = NULL, theta_MAV = NULL, theta_NULL = 0.10,
    pi_t = 0.15, pi_c = 0.15,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n_t = 12L, n_c = 12L,
    a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
    z = NULL, m_t = 30L, m_c = 30L,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL
  )
  expect_s3_class(result, "getgamma1bin")
  if (!is.na(result$PrGo_at_gamma_go))
    expect_true(result$PrGo_at_gamma_go >= 0 && result$PrGo_at_gamma_go <= 1)
  if (!is.na(result$PrNoGo_at_gamma_nogo))
    expect_true(result$PrNoGo_at_gamma_nogo >= 0 && result$PrNoGo_at_gamma_nogo <= 1)
})

test_that("getgamma1bin external design returns correct class", {
  result <- getgamma1bin(
    prob = 'posterior', design = 'external',
    theta_TV = 0.20, theta_MAV = 0.05, theta_NULL = NULL,
    pi_t = 0.15, pi_c = 0.15,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n_t = 12L, n_c = 12L,
    a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
    z = NULL, m_t = NULL, m_c = NULL,
    ne_t = 15L, ne_c = 15L, ye_t = 6L, ye_c = 4L, alpha0e_t = 0.5, alpha0e_c = 0.5
  )
  expect_s3_class(result, "getgamma1bin")
})

test_that("getgamma1bin PrGo_grid and PrNoGo_grid values in [0, 1]", {
  result <- getgamma1bin(
    prob = 'posterior', design = 'controlled',
    theta_TV = 0.20, theta_MAV = 0.05, theta_NULL = NULL,
    pi_t = 0.15, pi_c = 0.15,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n_t = 12L, n_c = 12L,
    a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
    z = NULL, m_t = NULL, m_c = NULL,
    ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL
  )
  expect_true(all(result$PrGo_grid >= 0 & result$PrGo_grid <= 1))
  expect_true(all(result$PrNoGo_grid >= 0 & result$PrNoGo_grid <= 1))
})

test_that("getgamma1bin input validation: invalid prob", {
  expect_error(getgamma1bin(
    prob = 'bayes', design = 'controlled',
    theta_TV = 0.20, theta_MAV = 0.05, theta_NULL = NULL,
    pi_t = 0.15, pi_c = 0.15,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 12L, n_c = 12L,
    a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5
  ))
})

# ---------------------------------------------------------------------------
# getgamma1cont
# ---------------------------------------------------------------------------

test_that("getgamma1cont posterior vague MM controlled returns correct class and structure", {
  result <- getgamma1cont(
    nsim = 50L, prob = 'posterior', design = 'controlled',
    prior = 'vague', CalcMethod = 'MM',
    theta_TV = 1.0, theta_MAV = 0.0, theta_NULL = NULL, nMC = NULL,
    mu_t = 1.0, mu_c = 0.0, sigma_t = 1.5, sigma_c = 1.5,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n_t = 10L, n_c = 10L, m_t = NULL, m_c = NULL,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
    r = NULL, ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL, seed = 1L
  )
  expect_s3_class(result, "getgamma1cont")
  expect_true(all(c("gamma_go", "gamma_nogo", "PrGo_at_gamma_go", "PrNoGo_at_gamma_nogo",
                    "gamma_grid", "PrGo_grid", "PrNoGo_grid") %in% names(result)))
})

test_that("getgamma1cont posterior vague MM grid values in [0, 1]", {
  result <- getgamma1cont(
    nsim = 50L, prob = 'posterior', design = 'controlled',
    prior = 'vague', CalcMethod = 'MM',
    theta_TV = 1.0, theta_MAV = 0.0, theta_NULL = NULL, nMC = NULL,
    mu_t = 1.0, mu_c = 0.0, sigma_t = 1.5, sigma_c = 1.5,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n_t = 10L, n_c = 10L, m_t = NULL, m_c = NULL,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
    r = NULL, ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL, seed = 1L
  )
  expect_true(all(result$PrGo_grid >= 0 & result$PrGo_grid <= 1))
  expect_true(all(result$PrNoGo_grid >= 0 & result$PrNoGo_grid <= 1))
  expect_equal(length(result$PrGo_grid), length(result$gamma_grid))
})

test_that("getgamma1cont posterior uncontrolled vague MM returns correct class", {
  result <- getgamma1cont(
    nsim = 50L, prob = 'posterior', design = 'uncontrolled',
    prior = 'vague', CalcMethod = 'MM',
    theta_TV = 1.0, theta_MAV = 0.0, theta_NULL = NULL, nMC = NULL,
    mu_t = 1.0, mu_c = NULL, sigma_t = 1.5, sigma_c = NULL,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n_t = 10L, n_c = NULL, m_t = NULL, m_c = NULL,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = 0.0, sigma0_t = NULL, sigma0_c = NULL,
    r = 1.0, ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL, seed = 2L
  )
  expect_s3_class(result, "getgamma1cont")
})

test_that("getgamma1cont predictive vague MM returns correct class", {
  result <- getgamma1cont(
    nsim = 50L, prob = 'predictive', design = 'controlled',
    prior = 'vague', CalcMethod = 'MM',
    theta_TV = NULL, theta_MAV = NULL, theta_NULL = 0.0, nMC = NULL,
    mu_t = 1.0, mu_c = 0.0, sigma_t = 1.5, sigma_c = 1.5,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n_t = 10L, n_c = 10L, m_t = 30L, m_c = 30L,
    kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
    mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
    r = NULL, ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
    bar_ye_t = NULL, se_t = NULL, bar_ye_c = NULL, se_c = NULL, seed = 3L
  )
  expect_s3_class(result, "getgamma1cont")
  if (!is.na(result$PrGo_at_gamma_go))
    expect_true(result$PrGo_at_gamma_go >= 0 && result$PrGo_at_gamma_go <= 1)
  if (!is.na(result$PrNoGo_at_gamma_nogo))
    expect_true(result$PrNoGo_at_gamma_nogo >= 0 && result$PrNoGo_at_gamma_nogo <= 1)
})

test_that("getgamma1cont input validation: invalid nsim", {
  expect_error(getgamma1cont(
    nsim = -1L, prob = 'posterior', design = 'controlled',
    prior = 'vague', CalcMethod = 'MM',
    theta_TV = 1.0, theta_MAV = 0.0, theta_NULL = NULL, nMC = NULL,
    mu_t = 1.0, mu_c = 0.0, sigma_t = 1.5, sigma_c = 1.5,
    target_go = 0.05, target_nogo = 0.20,
    n_t = 10L, n_c = 10L, seed = 1L
  ))
})
