test_that("stats_continuous2() builds the means and covariances", {
  # One study: compare with L z / sqrt(n) and L A A^T L^T / (n - 1)
  d <- list(z1 = 0.3, z2 = -1.2, x1 = 7.5, x2 = 4.1, z12 = 0.8)
  sd <- c(3, 2)
  rho <- 0.4
  s <- stats_continuous2(d, 9, c(1, -1), sd, rho)
  L <- t(chol(cov2(sd[1], sd[2], rho)))
  A <- matrix(c(sqrt(7.5), 0.8, 0, sqrt(4.1)), 2L)
  expect_equal(s$mean[1L, ], c(1, -1) + drop(L %*% c(0.3, -1.2)) / 3,
               tolerance = 1e-14)
  S <- L %*% A %*% t(A) %*% t(L) / 8
  expect_equal(unname(s$cov[1L, ]), c(S[1, 1], S[1, 2], S[2, 2]),
               tolerance = 1e-14)
  expect_identical(colnames(s$cov), c("c11", "c12", "c22"))
})

test_that("stats_continuous2() has the sampling distribution of the data", {
  set.seed(8)
  n <- 6
  d <- draws_continuous2(40000L, n)
  sig <- cov2(3, 2, -0.6)
  s <- stats_continuous2(d, n, c(5, 1), c(3, 2), -0.6)
  # Means: N(mu, Sigma / n); covariances: unbiased for Sigma, with
  # Var(c11) = 2 sigma1^4 / (n - 1) (Wishart)
  expect_equal(colMeans(s$mean), c(5, 1), tolerance = 0.03)
  expect_equal(c(stats::cov(s$mean)), c(sig / n), tolerance = 0.05)
  expect_equal(unname(colMeans(s$cov)), c(9, -3.6, 4), tolerance = 0.03)
  expect_equal(stats::var(s$cov[, "c11"]), 2 * 81 / (n - 1),
               tolerance = 0.05)
  # The sample means and covariances are independent
  expect_lt(abs(stats::cor(s$mean[, 1L], s$cov[, "c11"])), 0.03)
})
