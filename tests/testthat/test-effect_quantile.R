test_that("effect_quantile() finds the quantiles of the effect", {
  expect_equal(effect_quantile(stats::pnorm, 0.975, c(-Inf, Inf)),
               stats::qnorm(0.975), tolerance = 1e-6)
  # The interval is searched far from 0 when the range is infinite
  expect_equal(effect_quantile(function(q) stats::pnorm(q, 1000), 0.5,
                               c(-Inf, Inf)), 1000, tolerance = 1e-6)
  expect_equal(effect_quantile(function(q) stats::pnorm(q, -1000), 0.5,
                               c(-Inf, Inf)), -1000, tolerance = 1e-6)
  expect_equal(effect_quantile(function(q) stats::punif(q, -1, 1), 0.25,
                               c(-1, 1)), -0.5, tolerance = 1e-6)
  expect_equal(effect_quantile(stats::pexp, 0.5, c(0, Inf)), log(2),
               tolerance = 1e-6)
})
