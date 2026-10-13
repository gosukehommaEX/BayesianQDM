test_that("posterior_continuous() gives t(n - 1, ybar, s^2 / n) with the vague prior", {
  post <- posterior_continuous(prior_vague(), 12, c(18, -3), c(25, 4))
  expect_equal(post$mu, c(18, -3))
  expect_equal(post$kappa, 12)
  expect_equal(post$nu, 11)
  expect_equal(post$sigma2, c(25, 4) ^ 2)
})

test_that("posterior_continuous() matches the closed forms of the NIX posterior", {
  # Python_reference/ref_continuous_qdm.py, design C2: closed form with the
  # prior mean, the external data and the current data as three weighted
  # groups
  pt <- posterior_continuous(prior_nix(5, 2, 4, 20), 10, 22, 24,
                             external = borrow(data_continuous(6, 20, 30),
                                               0.5))
  expect_equal(pt$mu, 19.3333333333, tolerance = 1e-10)
  expect_equal(pt$kappa, 15)
  expect_equal(pt$nu, 17)
  expect_equal(sqrt(pt$sigma2), 23.6610045699, tolerance = 1e-10)
  pc <- posterior_continuous(prior_nix(-5, 1, 3, 25), 8, -8, 18,
                             external = borrow(data_continuous(30, -12, 20),
                                               0.4))
  expect_equal(pc$mu, -10.1428571429, tolerance = 1e-10)
  expect_equal(pc$kappa, 21)
  expect_equal(pc$nu, 23)
  expect_equal(sqrt(pc$sigma2), 19.6574705664, tolerance = 1e-10)
  # Without data, the prior parameters are returned
  p0 <- posterior_continuous(prior_nix(1, 2, 3, 4), 0)
  expect_identical(p0, list(mu = 1, kappa = 2, nu = 3, sigma2 = 16))
})

test_that("posterior_continuous() reproduces Huang et al. (submitted), Table 3", {
  # Vague prior, power prior with weight a0 for 60 adult subjects and 30
  # pediatric subjects. Closed form of Huang et al., Section 6.2:
  # kappa = 60 a0 + 30, nu = 60 a0 + 29,
  # nu sigma^2 = a0 S_x + S_y + 60 a0 30 / (60 a0 + 30) (ybar - ybar_e)^2.
  a0 <- c(0.25, 0.5, 0.75)
  pub_mean_c <- c(2.70, 2.65, 2.62)
  pub_sd_c <- c(0.640, 0.543, 0.480)
  pub_mean_t <- c(-0.57, -0.60, -0.62)
  pub_sd_t <- c(0.589, 0.504, 0.448)
  for (i in seq_along(a0)) {
    pc <- posterior_continuous(prior_vague(), 30, 2.8, 4.3,
                               external = borrow(data_continuous(60, 2.5, 4),
                                                 a0[i]))
    pt <- posterior_continuous(prior_vague(), 30, -0.5, 3.9,
                               external = borrow(data_continuous(60, -0.7,
                                                                 3.8),
                                                 a0[i]))
    k <- 60 * a0[i] + 30
    expect_equal(pc$kappa, k)
    expect_equal(pc$nu, k - 1)
    expect_equal(pc$mu, (60 * a0[i] * 2.5 + 30 * 2.8) / k)
    expect_equal(pc$sigma2, (a0[i] * 59 * 4 ^ 2 + 29 * 4.3 ^ 2 +
                               60 * a0[i] * 30 / k * (2.8 - 2.5) ^ 2) /
                   (k - 1))
    expect_equal(pt$sigma2, (a0[i] * 59 * 3.8 ^ 2 + 29 * 3.9 ^ 2 +
                               60 * a0[i] * 30 / k * (-0.5 + 0.7) ^ 2) /
                   (k - 1))
    # Posterior means (2 decimals) and standard deviations (3 decimals) of
    # the published table; the standard deviation of t(nu, mu, s^2) is
    # s sqrt(nu / (nu - 2)). Two of the treatment standard deviations differ
    # from the exact values (0.58846 and 0.44738) by up to 0.0006, slightly
    # more than rounding; the article summarizes Monte Carlo samples of the
    # posterior distributions (its Section 6.3).
    post_sd <- function(p) sqrt(p$sigma2 / p$kappa * p$nu / (p$nu - 2))
    expect_equal(round(pc$mu, 2), pub_mean_c[i])
    expect_equal(round(pt$mu, 2), pub_mean_t[i])
    expect_equal(round(post_sd(pc), 3), pub_sd_c[i])
    expect_lt(abs(post_sd(pt) - pub_sd_t[i]), 0.001)
  }
})

test_that("posterior_continuous() handles external data only and a0 = 0", {
  ext <- borrow(data_continuous(30, -12, 20), 0.4)
  # Huang et al., Section 6.1: NIX(ybar_e, a0 n_e, a0 n_e - 1,
  # a0 (n_e - 1) s_e^2 / (a0 n_e - 1)) with the vague prior
  p <- posterior_continuous(prior_vague(), 0, external = ext)
  expect_equal(p$mu, -12)
  expect_equal(p$kappa, 12)
  expect_equal(p$nu, 11)
  expect_equal(p$sigma2, 0.4 * 29 * 400 / 11)
  # a0 = 0 ignores the external data, and a small a0 is close to it
  none <- posterior_continuous(prior_vague(), 12, 3, 2)
  expect_identical(posterior_continuous(prior_vague(), 12, 3, 2,
                                        borrow(data_continuous(30, 9, 5), 0)),
                   none)
  tiny <- posterior_continuous(prior_vague(), 12, 3, 2,
                               borrow(data_continuous(30, 9, 5), 1e-8))
  expect_equal(unlist(tiny), unlist(none), tolerance = 1e-6)
})
