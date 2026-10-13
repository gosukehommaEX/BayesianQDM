test_that("control_posterior_continuous() uses the concurrent control data", {
  des <- ref_design_cont("C2")
  post_t <- posterior_continuous(des$prior_t, 10, 22, 24, des$external_t)
  expect_identical(
    control_posterior_continuous(des, post_t, 8, c(-8, -6), c(18, 20)),
    posterior_continuous(des$prior_c, 8, c(-8, -6), c(18, 20),
                         external = des$external_c)
  )
})

test_that("control_posterior_continuous() uses the external data only", {
  # Python_reference/ref_continuous_qdm.py, design C4
  des <- ref_design_cont("C4")
  post_c <- control_posterior_continuous(des, list())
  expect_equal(post_c$mu, -12)
  expect_equal(post_c$kappa, 12)
  expect_equal(post_c$nu, 11)
  expect_equal(sqrt(post_c$sigma2), 20.5382127221, tolerance = 1e-10)
})

test_that("control_posterior_continuous() builds the hypothetical control", {
  # Kang et al. (2026), Section 4.2.4:
  # mu_c ~ t(nu_nt, mu0c, r sigma_nt^2 / kappa_nt)
  des <- ref_design_cont("C3")
  post_t <- posterior_continuous(des$prior_t, 12, c(25, 30), c(20, 10))
  post_c <- control_posterior_continuous(des, post_t)
  expect_identical(post_c, list(mu = -10, kappa = 12, nu = 11,
                                sigma2 = 1.5 * c(400, 100)))
})
