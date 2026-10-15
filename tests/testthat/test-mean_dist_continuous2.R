# Reference: Python_reference/ref_continuous2_qdm.py, section C (posterior
# scale Lambda / (kappa (nu - 1)); predictive scale (1 / kappa + 1 / m)
# Lambda / (nu - 1), checked in section B against direct simulation of the
# hierarchical model).
test_that("mean_dist_continuous2() gives the posterior and predictive scales", {
  ref_t <- list(
    K1 = list(
      loc = c(18, 2.5), df = 10,
      scale = c(330, 11, 2.29166666666667)
    ),
    K2 = list(
      loc = c(4.98338870431894, 1.99667774086379), df = 33,
      scale = c(193.297585517318, 5.62116011161827, 1.08738514826748)
    ),
    K3 = list(
      loc = c(2.37257142857143, 1.00118633540373), df = 19,
      scale = c(1.67795549308551, 0.263590024497268, 0.60288915604117)
    ),
    K4 = list(
      loc = c(4, 1.5), df = 13,
      scale = c(64.6153846153846, 3.23076923076923, 0.646153846153846)
    ),
    K5 = list(
      loc = c(2.97520661157025, 1), df = 15,
      scale = c(376.17571659495, 10.0826446280992, 2.03168044077135)
    ),
    K6 = list(
      loc = c(7, 2), df = 13,
      scale = c(9.87179487179487, -0.592307692307692, 0.394871794871795)
    )
  )
  ref_c <- list(
    K1 = list(
      loc = c(-10, -1), df = 10,
      scale = c(449.166666666667, 7.7, 3.3)
    ),
    K2 = list(
      loc = c(-10, -1), df = 33,
      scale = c(170.668478808014, 4.28495922681969, 1.26975737440854)
    ),
    K3 = list(
      loc = c(-0.67910447761194, -1.32089552238806), df = 23,
      scale = c(1.30581977910577, 0.095846627989191, 0.600695357734365)
    ),
    K4 = list(
      loc = c(-10, -1), df = 13,
      scale = c(96.9230769230769, 4.84615384615385, 0.969230769230769)
    ),
    K5 = list(
      loc = c(-10, -1), df = 8,
      scale = c(581.875, 12.46875, 2.96875)
    ),
    K6 = list(
      loc = c(0, 0), df = 13,
      scale = c(21.5384615384615, -1.29230769230769, 0.861538461538462)
    )
  )
  for (name in names(ref_t)) {
    des <- ref_design_cont2(name)
    crit <- des$criterion
    d <- ref_data_cont2[[name]]
    pt <- posterior_continuous2(des$prior_t, d$t$n, d$t$mean, d$t$cov,
                                external = des$external_t)
    pc <- if (is.null(d$c)) {
      control_posterior_continuous2(des, pt)
    } else {
      control_posterior_continuous2(des, pt, d$c$n, d$c$mean, d$c$cov)
    }
    for (g in c("t", "c")) {
      dd <- mean_dist_continuous2(if (g == "t") pt else pc, crit$type,
                                  crit[[paste0("m_", g)]])
      e <- (if (g == "t") ref_t else ref_c)[[name]]
      expect_equal(dd$loc, matrix(e$loc, 1L), tolerance = 1e-12,
                   info = paste(name, g))
      expect_equal(unname(dd$scale), matrix(e$scale, 1L), tolerance = 1e-12,
                   info = paste(name, g))
      expect_equal(dd$df, e$df, info = paste(name, g))
    }
  }
  dd <- mean_dist_continuous2(pt, "posterior")
  expect_identical(colnames(dd$scale), c("V11", "V12", "V22"))
})
