# Reference: Python_reference/ref_continuous2_qdm.py, section A (control
# groups: concurrent with external data, hypothetical, external only).
test_that("control_posterior_continuous2() covers the three control forms", {
  ref <- list(
    K1 = list(
      mu = c(-10, -1), kappa = 12, nu = 11,
      Lambda = c(53900, 924, 396)
    ),
    K2 = list(
      mu = c(-10, -1), kappa = 30.1, nu = 34,
      Lambda = c(169525, 4256.25, 1261.25)
    ),
    K3 = list(
      mu = c(-0.67910447761194, -1.32089552238806), kappa = 20.1, nu = 24,
      Lambda = c(430.585223880597, 31.604776119403, 198.075223880597)
    ),
    K4 = list(
      mu = c(-10, -1), kappa = 15, nu = 14,
      Lambda = c(18900, 945, 189)
    ),
    K5 = list(
      mu = c(-10, -1), kappa = 10, nu = 9,
      Lambda = c(46550, 997.5, 237.5)
    ),
    K6 = list(
      mu = c(0, 0), kappa = 15, nu = 14,
      Lambda = c(2800, -168, 112)
    )
  )
  for (name in names(ref)) {
    des <- ref_design_cont2(name)
    d <- ref_data_cont2[[name]]
    pt <- posterior_continuous2(des$prior_t, d$t$n, d$t$mean, d$t$cov,
                                external = des$external_t)
    pc <- if (is.null(d$c)) {
      control_posterior_continuous2(des, pt)
    } else {
      control_posterior_continuous2(des, pt, d$c$n, d$c$mean, d$c$cov)
    }
    e <- ref[[name]]
    expect_equal(pc$mu, matrix(e$mu, 1L), tolerance = 1e-12, info = name)
    expect_equal(pc$kappa, e$kappa, tolerance = 1e-12, info = name)
    expect_equal(pc$nu, e$nu, tolerance = 1e-12, info = name)
    expect_equal(unname(pc$Lambda), matrix(e$Lambda, 1L),
                 tolerance = 1e-12, info = name)
  }
  # The hypothetical control follows the rows of the treatment parameters
  des <- ref_design_cont2("K4")
  pt <- posterior_continuous2(des$prior_t, 15, rbind(c(1, 2), c(3, 4)),
                              rbind(c(4, 1, 2), c(5, 1, 3)))
  pc <- control_posterior_continuous2(des, pt)
  expect_equal(pc$mu, rbind(c(-10, -1), c(-10, -1)))
  expect_equal(pc$Lambda, 1.5 * pt$Lambda)
})

test_that("control_posterior_continuous2() returns NULL for a fixed control", {
  expect_null(control_posterior_continuous2(ref_design_cont2("KF1"), list()))
})
