test_that("boundary_binary() specifies the grid of responders by sample size", {
  s <- boundary_binary(ref_design("D8"), c(20, 10, 10))
  expect_identical(s$cond, data.frame(n_t = c(10L, 20L)))
  expect_identical(s$values, list(0:10, 0:20))
  expect_true(s$discrete)
  expect_true(s$increasing)
  p <- s$prob(2L, 11)
  dec <- qdm_decide(ref_design("D8"), data_binary(20, 11))
  expect_equal(c(p$go, p$nogo), c(dec$prob_go, dec$prob_nogo))
  d <- s$derived(c(1L, 2L), c(3, 11))
  expect_equal(d$rate, c(0.3, 0.55))
  expect_equal(d$effect, c(0.1, 0.35))
  expect_identical(s$axes$x, "n_t")
  expect_null(s$axes$sec_scale)
  expect_match(s$titles[["effect"]], "fixed control rate \\(y_t / n_t - 0.2\\)")
})

test_that("boundary_binary() conditions on the control responders", {
  s <- boundary_binary(ref_design("D1"), 12)
  expect_identical(s$cond$y_c, 0:12)
  expect_identical(s$cond$n_c, rep(12L, 13))
  p <- s$prob(c(4L, 4L), c(9, 2))
  dec <- qdm_decide(ref_design("D1"), data_binary(12, 9), data_binary(12, 3))
  expect_equal(p$go[1L], dec$prob_go)
  expect_equal(p$nogo[1L], dec$prob_nogo)
  expect_equal(s$derived(4L, 9)$effect, 9 / 12 - 3 / 12)
  expect_identical(s$axes$x, "y_c")
  expect_identical(s$axes$sec_scale, 12L)
  expect_error(boundary_binary(ref_design("D1"), c(10, 12)),
               "'n_t' must be a single number for a design with a concurrent")
})

test_that("boundary_binary() measures the effect from the control distribution", {
  # Hypothetical control Beta(1 + 4, 1 + 8): mean 5 / 14
  s <- boundary_binary(ref_design("D3"), 12)
  expect_equal(s$derived(1L, 6)$effect, 6 / 12 - 5 / 14)
  expect_match(s$titles[["effect"]], "mean of the control distribution")
  # Smaller effects are favorable
  des <- design_binary(20, 0, prior_beta(1, 1),
                       criterion = criterion_posterior(0.1, 0.25, 0.8, 0.8,
                                                       direction = "less"),
                       fixed_c = 0)
  expect_false(boundary_binary(des, 20)$increasing)
})
