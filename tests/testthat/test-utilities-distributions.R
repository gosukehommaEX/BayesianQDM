# Tests for utility, sampling, and distribution functions

# ---------------------------------------------------------------------------
# rdirichlet
# ---------------------------------------------------------------------------

test_that("rdirichlet returns correct structure", {
  set.seed(1)
  result <- rdirichlet(n = 10L, alpha = c(1, 2, 3))
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(10L, 3L))
  expect_true(all(result >= 0))
  expect_true(all(abs(rowSums(result) - 1) < 1e-10))
})

test_that("rdirichlet row sums are 1 for various alpha", {
  set.seed(2)
  result <- rdirichlet(n = 100L, alpha = c(0.5, 0.5, 0.5, 0.5))
  expect_true(all(abs(rowSums(result) - 1) < 1e-10))
})

test_that("rdirichlet input validation works", {
  expect_error(rdirichlet(n = 0L,  alpha = c(1, 1, 1)))
  expect_error(rdirichlet(n = 10L, alpha = c(-1, 1, 1)))
  expect_error(rdirichlet(n = 10L, alpha = c(0, 1, 1)))
})

# ---------------------------------------------------------------------------
# rnsbt
# ---------------------------------------------------------------------------

test_that("rnsbt returns correct structure", {
  set.seed(1)
  mu <- c(2.0, -1.0)
  V  <- matrix(c(4.0, 1.2, 1.2, 1.0), 2, 2)
  result <- rnsbt(n = 500L, df = 10, mu = mu, V = V)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(500L, 2L))
})

test_that("rnsbt sample mean is close to mu for large n and large df", {
  set.seed(2)
  mu <- c(1.5, -0.5)
  V  <- matrix(c(1.0, 0.3, 0.3, 0.5), 2, 2)
  result <- rnsbt(n = 5000L, df = 100, mu = mu, V = V)
  expect_true(abs(mean(result[, 1L]) - mu[1L]) < 0.1)
  expect_true(abs(mean(result[, 2L]) - mu[2L]) < 0.1)
})

test_that("rnsbt input validation works", {
  V <- matrix(c(1, 0, 0, 1), 2, 2)
  expect_error(rnsbt(n = 0L,  df = 5,  mu = c(0, 0), V = V))
  expect_error(rnsbt(n = 10L, df = -1, mu = c(0, 0), V = V))
  expect_error(rnsbt(n = 10L, df = 5,  mu = c(0),    V = V))
  expect_error(rnsbt(n = 10L, df = 5,  mu = c(0, 0),
                     V = matrix(c(1, 0, 0, 1, 0, 0, 0, 0, 1), 3, 3)))
})

# ---------------------------------------------------------------------------
# getjointbin
# Returns a named numeric vector of length 4: (p00, p01, p10, p11)
# ---------------------------------------------------------------------------

test_that("getjointbin returns a named vector of length 4", {
  result <- getjointbin(pi1 = 0.3, pi2 = 0.4, rho = 0.0)
  expect_type(result, "double")
  expect_length(result, 4L)
  expect_named(result, c("p00", "p01", "p10", "p11"))
})

test_that("getjointbin elements are non-negative and sum to 1", {
  result <- getjointbin(pi1 = 0.3, pi2 = 0.4, rho = 0.0)
  expect_true(all(result >= 0))
  expect_equal(sum(result), 1, tolerance = 1e-10)
})

test_that("getjointbin with rho = 0: p11 equals product of marginals", {
  pi1 <- 0.3; pi2 <- 0.4
  result <- getjointbin(pi1 = pi1, pi2 = pi2, rho = 0.0)
  expect_equal(unname(result["p11"]), pi1 * pi2, tolerance = 1e-10)
})

test_that("getjointbin marginals are consistent with pi1 and pi2", {
  pi1 <- 0.3; pi2 <- 0.4
  result <- getjointbin(pi1 = pi1, pi2 = pi2, rho = 0.0)
  # P(Y1 = 1) = p10 + p11
  expect_equal(unname(result["p10"] + result["p11"]), pi1, tolerance = 1e-10)
  # P(Y2 = 1) = p01 + p11
  expect_equal(unname(result["p01"] + result["p11"]), pi2, tolerance = 1e-10)
})

test_that("getjointbin input validation works", {
  expect_error(getjointbin(pi1 = -0.1, pi2 = 0.3, rho = 0))
  expect_error(getjointbin(pi1 = 0.3,  pi2 = 1.1, rho = 0))
  expect_error(getjointbin(pi1 = 0.5,  pi2 = 0.5, rho = 2))
})

# ---------------------------------------------------------------------------
# allmultinom
# Signature: allmultinom(n) — always returns a 4-column matrix
# ---------------------------------------------------------------------------

test_that("allmultinom returns correct structure", {
  result <- allmultinom(n = 5L)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 4L)
  expect_named(result[1L, ], c("x00", "x01", "x10", "x11"))
  expect_true(all(rowSums(result) == 5L))
  expect_true(all(result >= 0L))
})

test_that("allmultinom row count matches stars-and-bars formula", {
  # C(n+3, 3) = (n+1)(n+2)(n+3)/6
  n <- 5L
  result <- allmultinom(n)
  expect_equal(nrow(result), choose(n + 3L, 3L))
})

test_that("allmultinom n = 0 returns 1 row of zeros", {
  result <- allmultinom(0L)
  expect_equal(nrow(result), 1L)
  expect_true(all(result == 0L))
})

test_that("allmultinom input validation works", {
  expect_error(allmultinom(-1L))
})

# ---------------------------------------------------------------------------
# pbetadiff
# ---------------------------------------------------------------------------

test_that("pbetadiff returns a value in [0, 1]", {
  result <- pbetadiff(q = 0.1, alpha1 = 6, beta1 = 5, alpha2 = 3, beta2 = 8)
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbetadiff lower.tail = TRUE and FALSE sum to 1", {
  q <- 0.1; a1 <- 6; b1 <- 5; a2 <- 3; b2 <- 8
  p_u <- pbetadiff(q, a1, b1, a2, b2, lower.tail = FALSE)
  p_l <- pbetadiff(q, a1, b1, a2, b2, lower.tail = TRUE)
  expect_equal(p_u + p_l, 1, tolerance = 1e-6)
})

test_that("pbetadiff at q = 0 with equal beta params gives ~0.5", {
  result <- pbetadiff(q = 0, alpha1 = 5, beta1 = 5, alpha2 = 5, beta2 = 5,
                      lower.tail = FALSE)
  expect_equal(result, 0.5, tolerance = 0.01)
})

test_that("pbetadiff input validation works", {
  expect_error(pbetadiff(q = 0.1, alpha1 = -1, beta1 = 5,  alpha2 = 3, beta2 = 8))
  expect_error(pbetadiff(q = 0.1, alpha1 = 6,  beta1 = 0,  alpha2 = 3, beta2 = 8))
  expect_error(pbetadiff(q = 0.1, alpha1 = 6,  beta1 = 5,  alpha2 = 3, beta2 = -1))
})

# ---------------------------------------------------------------------------
# pbetabinomdiff
# Signature: pbetabinomdiff(q, m1, m2, alpha1, alpha2, beta1, beta2, lower.tail)
# Posterior predictive: Y_k ~ BetaBinomial(m_k, alpha_k, beta_k)
# ---------------------------------------------------------------------------

test_that("pbetabinomdiff returns a value in [0, 1]", {
  result <- pbetabinomdiff(q = 0.1, m1 = 30, m2 = 30,
                           alpha1 = 10.5, alpha2 = 6.5,
                           beta1 = 10.5, beta2 = 14.5)
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("pbetabinomdiff lower.tail = TRUE and FALSE sum to 1", {
  p_u <- pbetabinomdiff(q = 0.1, m1 = 30, m2 = 30,
                        alpha1 = 10.5, alpha2 = 6.5,
                        beta1 = 10.5, beta2 = 14.5, lower.tail = FALSE)
  p_l <- pbetabinomdiff(q = 0.1, m1 = 30, m2 = 30,
                        alpha1 = 10.5, alpha2 = 6.5,
                        beta1 = 10.5, beta2 = 14.5, lower.tail = TRUE)
  expect_equal(p_u + p_l, 1, tolerance = 1e-10)
})

test_that("pbetabinomdiff input validation works", {
  expect_error(pbetabinomdiff(q = 0.1, m1 = 0,  m2 = 30,
                              alpha1 = 1, alpha2 = 1, beta1 = 1, beta2 = 1))
  expect_error(pbetabinomdiff(q = 0.1, m1 = 30, m2 = 30,
                              alpha1 = -1, alpha2 = 1, beta1 = 1, beta2 = 1))
})

# ---------------------------------------------------------------------------
# ptdiff_NI
# ---------------------------------------------------------------------------

test_that("ptdiff_NI returns a value in [0, 1]", {
  result <- ptdiff_NI(q = 1, mu.t1 = 3, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1,
                      nu.t1 = 15, nu.t2 = 15, lower.tail = FALSE)
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("ptdiff_NI symmetric case returns ~0.5", {
  result <- ptdiff_NI(q = 0, mu.t1 = 1, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1,
                      nu.t1 = 20, nu.t2 = 20, lower.tail = FALSE)
  expect_equal(result, 0.5, tolerance = 0.01)
})

test_that("ptdiff_NI lower.tail = TRUE and FALSE sum to 1", {
  p_u <- ptdiff_NI(q = 2, mu.t1 = 3, mu.t2 = 1, sd.t1 = 1.5, sd.t2 = 1.2,
                   nu.t1 = 15, nu.t2 = 18, lower.tail = FALSE)
  p_l <- ptdiff_NI(q = 2, mu.t1 = 3, mu.t2 = 1, sd.t1 = 1.5, sd.t2 = 1.2,
                   nu.t1 = 15, nu.t2 = 18, lower.tail = TRUE)
  expect_equal(p_u + p_l, 1, tolerance = 1e-6)
})

test_that("ptdiff_NI input validation works", {
  expect_error(ptdiff_NI(q = 1, mu.t1 = 3, mu.t2 = 1, sd.t1 = -1, sd.t2 = 1,
                         nu.t1 = 15, nu.t2 = 15))
  expect_error(ptdiff_NI(q = 1, mu.t1 = 3, mu.t2 = 1, sd.t1 = 1,  sd.t2 = 1,
                         nu.t1 = 1,  nu.t2 = 15))
})

# ---------------------------------------------------------------------------
# ptdiff_MC
# ---------------------------------------------------------------------------

test_that("ptdiff_MC returns a value in [0, 1]", {
  set.seed(1)
  result <- ptdiff_MC(nMC = 1000L, q = 1, mu.t1 = 3, mu.t2 = 1,
                      sd.t1 = 1, sd.t2 = 1, nu.t1 = 15, nu.t2 = 15,
                      lower.tail = FALSE)
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("ptdiff_MC symmetric case returns ~0.5", {
  set.seed(2)
  result <- ptdiff_MC(nMC = 5000L, q = 0, mu.t1 = 1, mu.t2 = 1,
                      sd.t1 = 1, sd.t2 = 1, nu.t1 = 20, nu.t2 = 20,
                      lower.tail = FALSE)
  expect_equal(result, 0.5, tolerance = 0.05)
})

test_that("ptdiff_MC and ptdiff_NI agree closely", {
  set.seed(3)
  p_mc <- ptdiff_MC(nMC = 10000L, q = 2, mu.t1 = 3, mu.t2 = 1,
                    sd.t1 = 1.5, sd.t2 = 1.2, nu.t1 = 15, nu.t2 = 18,
                    lower.tail = FALSE)
  p_ni <- ptdiff_NI(q = 2, mu.t1 = 3, mu.t2 = 1,
                    sd.t1 = 1.5, sd.t2 = 1.2, nu.t1 = 15, nu.t2 = 18,
                    lower.tail = FALSE)
  expect_equal(p_mc, p_ni, tolerance = 0.03)
})

test_that("ptdiff_MC input validation works", {
  expect_error(ptdiff_MC(nMC = 0L,    q = 1, mu.t1 = 1, mu.t2 = 0,
                         sd.t1 = 1, sd.t2 = 1, nu.t1 = 10, nu.t2 = 10))
  expect_error(ptdiff_MC(nMC = 1000L, q = 1, mu.t1 = 1, mu.t2 = 0,
                         sd.t1 = 1, sd.t2 = 1, nu.t1 = 1,  nu.t2 = 10))
})

# ---------------------------------------------------------------------------
# ptdiff_MM
# ---------------------------------------------------------------------------

test_that("ptdiff_MM returns a value in [0, 1]", {
  result <- ptdiff_MM(q = 1, mu.t1 = 3, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1,
                      nu.t1 = 15, nu.t2 = 15, lower.tail = FALSE)
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(result >= 0 && result <= 1)
})

test_that("ptdiff_MM symmetric case returns ~0.5", {
  result <- ptdiff_MM(q = 0, mu.t1 = 1, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1,
                      nu.t1 = 20, nu.t2 = 20, lower.tail = FALSE)
  expect_equal(result, 0.5, tolerance = 0.01)
})

test_that("ptdiff_MM and ptdiff_NI agree closely", {
  p_mm <- ptdiff_MM(q = 2, mu.t1 = 3, mu.t2 = 1,
                    sd.t1 = 1.5, sd.t2 = 1.2, nu.t1 = 15, nu.t2 = 18,
                    lower.tail = FALSE)
  p_ni <- ptdiff_NI(q = 2, mu.t1 = 3, mu.t2 = 1,
                    sd.t1 = 1.5, sd.t2 = 1.2, nu.t1 = 15, nu.t2 = 18,
                    lower.tail = FALSE)
  expect_equal(p_mm, p_ni, tolerance = 0.02)
})

test_that("ptdiff_MM vectorised input returns correct length", {
  result <- ptdiff_MM(q = 1, mu.t1 = c(2, 3, 4), mu.t2 = c(0, 1, 2),
                      sd.t1 = c(1, 1.2, 1.5), sd.t2 = c(1, 1.1, 1.3),
                      nu.t1 = 10, nu.t2 = 10, lower.tail = FALSE)
  expect_length(result, 3L)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("ptdiff_MM input validation works", {
  expect_error(ptdiff_MM(q = 1, mu.t1 = 1, mu.t2 = 0, sd.t1 = 1,  sd.t2 = 1,
                         nu.t1 = 4, nu.t2 = 10))
  expect_error(ptdiff_MM(q = 1, mu.t1 = 1, mu.t2 = 0, sd.t1 = -1, sd.t2 = 1,
                         nu.t1 = 10, nu.t2 = 10))
})
