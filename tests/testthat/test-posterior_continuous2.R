# Reference: Python_reference/ref_continuous2_qdm.py, section A (closed form
# with three weighted groups: prior mean, external data and current data;
# equal to the formulas of Yamaguchi et al., submitted, A.10 and A.2.5).
test_that("posterior_continuous2() matches the closed form", {
  ref <- list(
    K1 = list(
      mu = c(18, 2.5), kappa = 12, nu = 11,
      Lambda = c(39600, 1320, 275)
    ),
    K2 = list(
      mu = c(4.98338870431894, 1.99667774086379), kappa = 30.1, nu = 34,
      Lambda = c(192002.491694352, 5583.49833887043, 1080.09966777409)
    ),
    K3 = list(
      mu = c(2.37257142857143, 1.00118633540373), kappa = 16.1, nu = 20,
      Lambda = c(388.265193142857, 60.9925782857143, 139.503625440994)
    ),
    K4 = list(
      mu = c(4, 1.5), kappa = 15, nu = 14,
      Lambda = c(12600, 630, 126)
    ),
    K5 = list(
      mu = c(2.97520661157025, 1), kappa = 12.1, nu = 16,
      Lambda = c(68275.8925619835, 1830, 368.75)
    ),
    K6 = list(
      mu = c(7, 2), kappa = 15, nu = 14,
      Lambda = c(1400, -84, 56)
    )
  )
  for (name in names(ref)) {
    des <- ref_design_cont2(name)
    d <- ref_data_cont2[[name]]$t
    p <- posterior_continuous2(des$prior_t, d$n, d$mean, d$cov,
                               external = des$external_t)
    e <- ref[[name]]
    expect_equal(p$mu, matrix(e$mu, 1L), tolerance = 1e-12, info = name)
    expect_equal(p$kappa, e$kappa, tolerance = 1e-12, info = name)
    expect_equal(p$nu, e$nu, tolerance = 1e-12, info = name)
    expect_equal(unname(p$Lambda), matrix(e$Lambda, 1L), tolerance = 1e-12,
                 info = name)
  }
})

test_that("posterior_continuous2() updates many studies at once", {
  des <- ref_design_cont2("K3")
  d <- ref_data_cont2$K3$t
  one <- posterior_continuous2(des$prior_t, d$n, d$mean, d$cov,
                               external = des$external_t)
  cov3 <- c(d$cov[1L, 1L], d$cov[1L, 2L], d$cov[2L, 2L])
  many <- posterior_continuous2(des$prior_t, d$n,
                                rbind(d$mean, d$mean + c(1, -1)),
                                rbind(cov3, cov3 * 2),
                                external = des$external_t)
  expect_identical(dim(many$mu), c(2L, 2L))
  expect_identical(colnames(many$Lambda), c("L11", "L12", "L22"))
  expect_equal(many$mu[1L, ], one$mu[1L, ], tolerance = 1e-14)
  expect_equal(many$Lambda[1L, ], one$Lambda[1L, ], tolerance = 1e-14)
  other <- posterior_continuous2(des$prior_t, d$n, d$mean + c(1, -1),
                                 2 * d$cov, external = des$external_t)
  expect_equal(many$Lambda[2L, ], other$Lambda[1L, ], tolerance = 1e-14)
})

test_that("posterior_continuous2() handles the vague prior and a0 = 0", {
  # Vague prior: mu = ybar, kappa = n, nu = n - 1, Lambda = (n - 1) cov
  cv <- cov2(3, 2, 0.5)
  p <- posterior_continuous2(prior_vague(), 10, c(1, 2), cv)
  expect_equal(p$mu, matrix(c(1, 2), 1L))
  expect_identical(c(p$kappa, p$nu), c(10, 9))
  expect_equal(unname(p$Lambda[1L, ]), 9 * c(9, 3, 4))
  # External data with a0 = 0 are ignored
  q <- posterior_continuous2(prior_vague(), 10, c(1, 2), cv,
                             external = borrow(data_continuous2(5, c(0, 0),
                                                                cv), 0))
  expect_identical(q, p)
  # No current data: the prior
  pr <- prior_niw(c(1, -1), 2, 5, lambda_cs4)
  r0 <- posterior_continuous2(pr, 0)
  expect_equal(r0$mu, matrix(c(1, -1), 1L))
  expect_equal(unname(r0$Lambda[1L, ]), c(25, 1, 9))
  expect_identical(c(r0$kappa, r0$nu), c(2, 5))
})
