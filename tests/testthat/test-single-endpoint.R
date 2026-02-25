# Tests for single-endpoint functions (binary and continuous)

# ---------------------------------------------------------------------------
# pbayespostpred1bin
# Signature: pbayespostpred1bin(prob, design, theta0, n1, n2, y1, y2,
#            a1, a2, b1, b2, m1, m2, z, ne1, ne2, ye1, ye2, ae1, ae2,
#            lower.tail)
# Note: no nMC argument; y1 and y2 must have the same length
# ---------------------------------------------------------------------------

test_that("pbayespostpred1bin posterior controlled returns scalar in [0, 1]", {
  result <- pbayespostpred1bin(
    prob = 'posterior', design = 'controlled', theta0 = 0.15,
    n1 = 20, n2 = 20, y1 = 12, y2 = 8,
    a1 = 0.5, b1 = 0.5, a2 = 0.5, b2 = 0.5,
    m1 = NULL, m2 = NULL, z = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
    lower.tail = FALSE
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1bin predictive controlled returns scalar in [0, 1]", {
  result <- pbayespostpred1bin(
    prob = 'predictive', design = 'controlled', theta0 = 0.15,
    n1 = 20, n2 = 20, y1 = 12, y2 = 8,
    a1 = 0.5, b1 = 0.5, a2 = 0.5, b2 = 0.5,
    m1 = 30, m2 = 30, z = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
    lower.tail = FALSE
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1bin posterior uncontrolled returns scalar in [0, 1]", {
  result <- pbayespostpred1bin(
    prob = 'posterior', design = 'uncontrolled', theta0 = 0.15,
    n1 = 20, n2 = 20, y1 = 12, y2 = NULL,
    a1 = 0.5, b1 = 0.5, a2 = 0.5, b2 = 0.5,
    m1 = NULL, m2 = NULL, z = 5L,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
    lower.tail = FALSE
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1bin posterior external returns scalar in [0, 1]", {
  result <- pbayespostpred1bin(
    prob = 'posterior', design = 'external', theta0 = 0.15,
    n1 = 20, n2 = 20, y1 = 12, y2 = 8,
    a1 = 0.5, b1 = 0.5, a2 = 0.5, b2 = 0.5,
    m1 = NULL, m2 = NULL, z = NULL,
    ne1 = 30, ne2 = 30, ye1 = 10, ye2 = 6, ae1 = 0.5, ae2 = 0.5,
    lower.tail = FALSE
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1bin vectorised y1 and y2 returns correct length", {
  result <- pbayespostpred1bin(
    prob = 'posterior', design = 'controlled', theta0 = 0.15,
    n1 = 20, n2 = 20, y1 = 8:12, y2 = rep(6, 5),
    a1 = 0.5, b1 = 0.5, a2 = 0.5, b2 = 0.5,
    m1 = NULL, m2 = NULL, z = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
    lower.tail = FALSE
  )
  expect_length(result, 5L)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("pbayespostpred1bin input validation works", {
  # y1 > n1
  expect_error(pbayespostpred1bin(
    prob = 'posterior', design = 'controlled', theta0 = 0.15,
    n1 = 20, n2 = 20, y1 = 25, y2 = 8,
    a1 = 0.5, b1 = 0.5, a2 = 0.5, b2 = 0.5,
    m1 = NULL, m2 = NULL, z = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  ))
  # m1 missing for predictive
  expect_error(pbayespostpred1bin(
    prob = 'predictive', design = 'controlled', theta0 = 0.15,
    n1 = 20, n2 = 20, y1 = 12, y2 = 8,
    a1 = 0.5, b1 = 0.5, a2 = 0.5, b2 = 0.5,
    m1 = NULL, m2 = NULL, z = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  ))
})

# ---------------------------------------------------------------------------
# pbayespostpred1cont
# ---------------------------------------------------------------------------

test_that("pbayespostpred1cont posterior vague NI returns scalar in [0, 1]", {
  result <- pbayespostpred1cont(
    prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
    theta0 = 1, nMC = NULL,
    n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
    bar.y1 = 3, s1 = 1.5, bar.y2 = 1, s2 = 1.2, r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, se1 = NULL, bar.ye2 = NULL, se2 = NULL
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1cont posterior vague MM returns scalar in [0, 1]", {
  result <- pbayespostpred1cont(
    prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'MM',
    theta0 = 1, nMC = NULL,
    n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
    bar.y1 = 3, s1 = 1.5, bar.y2 = 1, s2 = 1.2, r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, se1 = NULL, bar.ye2 = NULL, se2 = NULL
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1cont NI and MM agree closely", {
  p_ni <- pbayespostpred1cont(
    prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
    theta0 = 1, nMC = NULL,
    n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
    bar.y1 = 3, s1 = 1.5, bar.y2 = 1, s2 = 1.2, r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, se1 = NULL, bar.ye2 = NULL, se2 = NULL
  )
  p_mm <- pbayespostpred1cont(
    prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'MM',
    theta0 = 1, nMC = NULL,
    n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
    bar.y1 = 3, s1 = 1.5, bar.y2 = 1, s2 = 1.2, r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, se1 = NULL, bar.ye2 = NULL, se2 = NULL
  )
  expect_equal(p_ni, p_mm, tolerance = 0.02)
})

test_that("pbayespostpred1cont predictive NI returns scalar in [0, 1]", {
  result <- pbayespostpred1cont(
    prob = 'predictive', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
    theta0 = 1, nMC = NULL,
    n1 = 15, n2 = 15, m1 = 30, m2 = 30,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
    bar.y1 = 3, s1 = 1.5, bar.y2 = 1, s2 = 1.2, r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, se1 = NULL, bar.ye2 = NULL, se2 = NULL
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1cont N-Inv-Chisq prior returns scalar in [0, 1]", {
  result <- pbayespostpred1cont(
    prob = 'posterior', design = 'controlled', prior = 'N-Inv-Chisq',
    CalcMethod = 'NI', theta0 = 1, nMC = NULL,
    n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
    kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5,
    mu01 = 3, mu02 = 1, sigma01 = 1.5, sigma02 = 1.2,
    bar.y1 = 3, s1 = 1.5, bar.y2 = 1, s2 = 1.2, r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, se1 = NULL, bar.ye2 = NULL, se2 = NULL
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbayespostpred1cont input validation works", {
  # m1 missing for predictive
  expect_error(pbayespostpred1cont(
    prob = 'predictive', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
    theta0 = 1, nMC = NULL,
    n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
    bar.y1 = 3, s1 = 1.5, bar.y2 = 1, s2 = 1.2, r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, se1 = NULL, bar.ye2 = NULL, se2 = NULL
  ))
})

# ---------------------------------------------------------------------------
# pbayesdecisionprob1bin
# Note: no nMC argument
# ---------------------------------------------------------------------------

test_that("pbayesdecisionprob1bin posterior controlled returns correct class", {
  result <- pbayesdecisionprob1bin(
    prob = 'posterior', design = 'controlled',
    theta.TV = 0.2, theta.MAV = 0.05, theta.NULL = NULL,
    gamma1 = 0.8, gamma2 = 0.2,
    pi1 = c(0.3, 0.5), pi2 = rep(0.2, 2),
    n1 = 15, n2 = 15,
    a1 = 0.5, b1 = 0.5, a2 = 0.5, b2 = 0.5,
    z = NULL, m1 = NULL, m2 = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  expect_s3_class(result, "pbayesdecisionprob1bin")
  df <- as.data.frame(result)
  expect_true(all(c("pi1", "Go", "Gray", "NoGo") %in% names(df)))
  expect_equal(nrow(df), 2L)
  expect_true(all(df$Go  >= 0 & df$Go  <= 1))
  expect_true(all(df$NoGo >= 0 & df$NoGo <= 1))
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

test_that("pbayesdecisionprob1bin predictive controlled returns correct class", {
  result <- pbayesdecisionprob1bin(
    prob = 'predictive', design = 'controlled',
    theta.TV = NULL, theta.MAV = NULL, theta.NULL = 0.15,
    gamma1 = 0.8, gamma2 = 0.2,
    pi1 = c(0.3, 0.5), pi2 = rep(0.2, 2),
    n1 = 15, n2 = 15,
    a1 = 0.5, b1 = 0.5, a2 = 0.5, b2 = 0.5,
    z = NULL, m1 = 40, m2 = 40,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  expect_s3_class(result, "pbayesdecisionprob1bin")
  df <- as.data.frame(result)
  expect_equal(nrow(df), 2L)
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

test_that("pbayesdecisionprob1bin input validation works", {
  # theta.TV missing for posterior
  expect_error(pbayesdecisionprob1bin(
    prob = 'posterior', design = 'controlled',
    theta.TV = NULL, theta.MAV = 0.05, theta.NULL = NULL,
    gamma1 = 0.8, gamma2 = 0.2,
    pi1 = 0.3, pi2 = 0.2, n1 = 15, n2 = 15,
    a1 = 0.5, b1 = 0.5, a2 = 0.5, b2 = 0.5,
    z = NULL, m1 = NULL, m2 = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  ))
  # gamma2 >= gamma1
  expect_error(pbayesdecisionprob1bin(
    prob = 'posterior', design = 'controlled',
    theta.TV = 0.2, theta.MAV = 0.05, theta.NULL = NULL,
    gamma1 = 0.5, gamma2 = 0.8,
    pi1 = 0.3, pi2 = 0.2, n1 = 15, n2 = 15,
    a1 = 0.5, b1 = 0.5, a2 = 0.5, b2 = 0.5,
    z = NULL, m1 = NULL, m2 = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  ))
})

# ---------------------------------------------------------------------------
# pbayesdecisionprob1cont
# ---------------------------------------------------------------------------

test_that("pbayesdecisionprob1cont posterior vague MM returns correct class", {
  set.seed(1)
  result <- pbayesdecisionprob1cont(
    nsim = 20L, prob = 'posterior', design = 'controlled',
    prior = 'vague', CalcMethod = 'MM',
    theta.TV = 1.5, theta.MAV = 0.5, theta.NULL = NULL, nMC = NULL,
    gamma1 = 0.8, gamma2 = 0.2,
    n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
    mu1 = c(2, 3), mu2 = c(0, 0), sigma1 = 1.5, sigma2 = 1.2, r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, se1 = NULL, bar.ye2 = NULL, se2 = NULL, seed = 1
  )
  expect_s3_class(result, "pbayesdecisionprob1cont")
  df <- as.data.frame(result)
  expect_true(all(c("mu1", "mu2", "Go", "Gray", "NoGo") %in% names(df)))
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
    theta.TV = NULL, theta.MAV = NULL, theta.NULL = 1, nMC = NULL,
    gamma1 = 0.8, gamma2 = 0.2,
    n1 = 15, n2 = 15, m1 = 30, m2 = 30,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
    mu1 = c(2, 3), mu2 = c(0, 0), sigma1 = 1.5, sigma2 = 1.2, r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, se1 = NULL, bar.ye2 = NULL, se2 = NULL, seed = 2
  )
  expect_s3_class(result, "pbayesdecisionprob1cont")
  df <- as.data.frame(result)
  expect_equal(nrow(df), 2L)
  expect_true(all(abs(df$Go + df$Gray + df$NoGo - 1) < 1e-6))
})

test_that("pbayesdecisionprob1cont input validation works", {
  # theta.TV missing for posterior
  expect_error(pbayesdecisionprob1cont(
    nsim = 10L, prob = 'posterior', design = 'controlled',
    prior = 'vague', CalcMethod = 'MM',
    theta.TV = NULL, theta.MAV = 0.5, theta.NULL = NULL, nMC = NULL,
    gamma1 = 0.8, gamma2 = 0.2,
    n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
    mu1 = 2, mu2 = 0, sigma1 = 1.5, sigma2 = 1.2, r = NULL,
    ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, se1 = NULL, bar.ye2 = NULL, se2 = NULL, seed = 1
  ))
})


# ---------------------------------------------------------------------------
# getgamma1bin
# ---------------------------------------------------------------------------

test_that("getgamma1bin posterior controlled returns correct class and structure", {
  result <- getgamma1bin(
    prob = 'posterior', design = 'controlled',
    theta.TV = 0.20, theta.MAV = 0.05, theta.NULL = NULL,
    pi1 = 0.15, pi2 = 0.15,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 12L, n2 = 12L,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    z = NULL, m1 = NULL, m2 = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  expect_s3_class(result, "getgamma1bin")
  expect_true(all(c("gamma1", "gamma2", "PrGo_at_gamma1", "PrNoGo_at_gamma2",
                    "gamma_grid", "PrGo_grid", "PrNoGo_grid") %in% names(result)))
})

test_that("getgamma1bin posterior controlled gamma values in (0, 1) or NA", {
  result <- getgamma1bin(
    prob = 'posterior', design = 'controlled',
    theta.TV = 0.20, theta.MAV = 0.05, theta.NULL = NULL,
    pi1 = 0.15, pi2 = 0.15,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 12L, n2 = 12L,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    z = NULL, m1 = NULL, m2 = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  if (!is.na(result$gamma1)) expect_true(result$gamma1 > 0 && result$gamma1 < 1)
  if (!is.na(result$gamma2)) expect_true(result$gamma2 > 0 && result$gamma2 < 1)
  expect_equal(length(result$PrGo_grid), length(result$gamma_grid))
  expect_equal(length(result$PrNoGo_grid), length(result$gamma_grid))
})

test_that("getgamma1bin posterior uncontrolled returns correct class", {
  result <- getgamma1bin(
    prob = 'posterior', design = 'uncontrolled',
    theta.TV = 0.20, theta.MAV = 0.05, theta.NULL = NULL,
    pi1 = 0.15, pi2 = NULL,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 12L, n2 = 12L,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    z = 3L, m1 = NULL, m2 = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  expect_s3_class(result, "getgamma1bin")
})

test_that("getgamma1bin predictive controlled returns correct class", {
  result <- getgamma1bin(
    prob = 'predictive', design = 'controlled',
    theta.TV = NULL, theta.MAV = NULL, theta.NULL = 0.10,
    pi1 = 0.15, pi2 = 0.15,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 12L, n2 = 12L,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    z = NULL, m1 = 30L, m2 = 30L,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  expect_s3_class(result, "getgamma1bin")
  if (!is.na(result$PrGo_at_gamma1))
    expect_true(result$PrGo_at_gamma1 >= 0 && result$PrGo_at_gamma1 <= 1)
  if (!is.na(result$PrNoGo_at_gamma2))
    expect_true(result$PrNoGo_at_gamma2 >= 0 && result$PrNoGo_at_gamma2 <= 1)
})

test_that("getgamma1bin external design returns correct class", {
  result <- getgamma1bin(
    prob = 'posterior', design = 'external',
    theta.TV = 0.20, theta.MAV = 0.05, theta.NULL = NULL,
    pi1 = 0.15, pi2 = 0.15,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 12L, n2 = 12L,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    z = NULL, m1 = NULL, m2 = NULL,
    ne1 = 15L, ne2 = 15L, ye1 = 6L, ye2 = 4L, ae1 = 0.5, ae2 = 0.5
  )
  expect_s3_class(result, "getgamma1bin")
})

test_that("getgamma1bin PrGo_grid and PrNoGo_grid values in [0, 1]", {
  result <- getgamma1bin(
    prob = 'posterior', design = 'controlled',
    theta.TV = 0.20, theta.MAV = 0.05, theta.NULL = NULL,
    pi1 = 0.15, pi2 = 0.15,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 12L, n2 = 12L,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    z = NULL, m1 = NULL, m2 = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  expect_true(all(result$PrGo_grid >= 0 & result$PrGo_grid <= 1))
  expect_true(all(result$PrNoGo_grid >= 0 & result$PrNoGo_grid <= 1))
})

test_that("getgamma1bin input validation: invalid prob", {
  expect_error(getgamma1bin(
    prob = 'bayes', design = 'controlled',
    theta.TV = 0.20, theta.MAV = 0.05, theta.NULL = NULL,
    pi1 = 0.15, pi2 = 0.15,
    target_go = 0.05, target_nogo = 0.20,
    n1 = 12L, n2 = 12L,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5
  ))
})

# ---------------------------------------------------------------------------
# getgamma1cont
# ---------------------------------------------------------------------------

test_that("getgamma1cont posterior vague MM controlled returns correct class and structure", {
  result <- getgamma1cont(
    nsim = 50L, prob = 'posterior', design = 'controlled',
    prior = 'vague', CalcMethod = 'MM',
    theta.TV = 1.0, theta.MAV = 0.0, theta.NULL = NULL, nMC = NULL,
    mu1 = 1.0, mu2 = 0.0, sigma1 = 1.5, sigma2 = 1.5,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 10L, n2 = 10L, m1 = NULL, m2 = NULL,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
    r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, se1 = NULL, bar.ye2 = NULL, se2 = NULL, seed = 1L
  )
  expect_s3_class(result, "getgamma1cont")
  expect_true(all(c("gamma1", "gamma2", "PrGo_at_gamma1", "PrNoGo_at_gamma2",
                    "gamma_grid", "PrGo_grid", "PrNoGo_grid") %in% names(result)))
})

test_that("getgamma1cont posterior vague MM grid values in [0, 1]", {
  result <- getgamma1cont(
    nsim = 50L, prob = 'posterior', design = 'controlled',
    prior = 'vague', CalcMethod = 'MM',
    theta.TV = 1.0, theta.MAV = 0.0, theta.NULL = NULL, nMC = NULL,
    mu1 = 1.0, mu2 = 0.0, sigma1 = 1.5, sigma2 = 1.5,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 10L, n2 = 10L, m1 = NULL, m2 = NULL,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
    r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, se1 = NULL, bar.ye2 = NULL, se2 = NULL, seed = 1L
  )
  expect_true(all(result$PrGo_grid >= 0 & result$PrGo_grid <= 1))
  expect_true(all(result$PrNoGo_grid >= 0 & result$PrNoGo_grid <= 1))
  expect_equal(length(result$PrGo_grid), length(result$gamma_grid))
})

test_that("getgamma1cont posterior uncontrolled vague MM returns correct class", {
  result <- getgamma1cont(
    nsim = 50L, prob = 'posterior', design = 'uncontrolled',
    prior = 'vague', CalcMethod = 'MM',
    theta.TV = 1.0, theta.MAV = 0.0, theta.NULL = NULL, nMC = NULL,
    mu1 = 1.0, mu2 = NULL, sigma1 = 1.5, sigma2 = NULL,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 10L, n2 = 10L, m1 = NULL, m2 = NULL,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = 0.0, sigma01 = NULL, sigma02 = NULL,
    r = 1.0, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, se1 = NULL, bar.ye2 = NULL, se2 = NULL, seed = 2L
  )
  expect_s3_class(result, "getgamma1cont")
})

test_that("getgamma1cont predictive vague MM returns correct class", {
  result <- getgamma1cont(
    nsim = 50L, prob = 'predictive', design = 'controlled',
    prior = 'vague', CalcMethod = 'MM',
    theta.TV = NULL, theta.MAV = NULL, theta.NULL = 0.0, nMC = NULL,
    mu1 = 1.0, mu2 = 0.0, sigma1 = 1.5, sigma2 = 1.5,
    target_go = 0.05, target_nogo = 0.20,
    crit_go = '<', crit_nogo = '<',
    sel_go = 'smallest', sel_nogo = 'largest',
    n1 = 10L, n2 = 10L, m1 = 30L, m2 = 30L,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
    r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, se1 = NULL, bar.ye2 = NULL, se2 = NULL, seed = 3L
  )
  expect_s3_class(result, "getgamma1cont")
  if (!is.na(result$PrGo_at_gamma1))
    expect_true(result$PrGo_at_gamma1 >= 0 && result$PrGo_at_gamma1 <= 1)
  if (!is.na(result$PrNoGo_at_gamma2))
    expect_true(result$PrNoGo_at_gamma2 >= 0 && result$PrNoGo_at_gamma2 <= 1)
})

test_that("getgamma1cont input validation: invalid nsim", {
  expect_error(getgamma1cont(
    nsim = -1L, prob = 'posterior', design = 'controlled',
    prior = 'vague', CalcMethod = 'MM',
    theta.TV = 1.0, theta.MAV = 0.0, theta.NULL = NULL, nMC = NULL,
    mu1 = 1.0, mu2 = 0.0, sigma1 = 1.5, sigma2 = 1.5,
    target_go = 0.05, target_nogo = 0.20,
    n1 = 10L, n2 = 10L, seed = 1L
  ))
})
