test_that("boundary_continuous() searches a grid of observed means", {
  des <- ref_design_cont("C1")
  s <- boundary_continuous(des, 12, NULL, 25, ref_data_cont$C1$c)
  expect_identical(s$cond, data.frame(n_t = 12L, sd_t = 25))
  expect_false(s$discrete)
  expect_length(s$values[[1L]], 205L)
  # The default grid covers the boundary values (references of
  # Python_reference/ref_boundary_qdm.py)
  expect_true(s$range[1L] < -2.69 && s$range[2L] > 20.74)
  # Two values far beyond each end of the grid
  span <- diff(s$range)
  expect_equal(s$values[[1L]][c(1, 2, 204, 205)],
               c(s$range[1L] - c(1e4, 1e2) * span,
                 s$range[2L] + c(1e2, 1e4) * span))
  p <- s$prob(1L, 18)
  dec <- qdm_decide(des, data_continuous(12, 18, 25), ref_data_cont$C1$c)
  expect_equal(c(p$go, p$nogo), c(dec$prob_go, dec$prob_nogo))
  expect_equal(s$derived(1L, 18)$effect, 30)
  expect_identical(s$axes$y, "effect")
  expect_identical(s$axes$offset, -12)
  expect_match(s$titles[["effect"]], "control mean -12")
})

test_that("boundary_continuous() takes a grid and several sample sizes", {
  des <- ref_design_cont("C7")
  s <- boundary_continuous(des, c(20, 10), seq(-10, 40, by = 5), 20, NULL)
  expect_identical(s$cond$n_t, c(10L, 20L))
  expect_identical(s$range, c(-10, 40))
  expect_identical(s$axes$x, "n_t")
  expect_true(s$axes$discrete_x)
  expect_identical(s$axes$offset, 0)
  p <- s$prob(c(1L, 2L), c(25, 25))
  expect_equal(p$go[1L], qdm_decide(des, data_continuous(10, 25, 20))$prob_go)
  expect_equal(p$go[2L], qdm_decide(des, data_continuous(20, 25, 20))$prob_go)
  # Effect from the fixed control mean -10 and the hypothetical mean -10
  expect_equal(s$derived(1L, 25)$effect, 35)
  s3 <- boundary_continuous(ref_design_cont("C3"), 12, 0, c(20, 30), NULL)
  expect_identical(s3$axes$x, "sd_t")
  expect_equal(s3$derived(1L, 25)$effect, 35)
  expect_match(s3$titles[["effect"]], "hypothetical control")
  # A single value of the grid
  expect_identical(s3$range, c(0, 0))
  s4 <- boundary_continuous(ref_design_cont("C4"), 12, NULL, 20, NULL)
  expect_match(s4$titles[["effect"]], "posterior mean of the control mean")
  # Predictive criterion: the default grid covers the null value and the
  # closed-form boundary values of test-qdm_boundary.R (6.697 and 4.183)
  s8 <- boundary_continuous(ref_design_cont("C8"), 15, NULL, 10, NULL)
  expect_true(s8$range[1L] < 4.18 && s8$range[2L] > 6.70)
})

test_that("boundary_continuous() validates its arguments", {
  des <- ref_design_cont("C1")
  dc <- ref_data_cont$C1$c
  expect_error(boundary_continuous(des, 12, NULL, NULL, dc),
               "'sd_t' must be given for a continuous endpoint.")
  expect_error(boundary_continuous(des, c(10, 12), NULL, c(20, 25), dc),
               "Only one of 'n_t' and 'sd_t' may have more than one value.")
  expect_error(boundary_continuous(des, 12, NULL, 25, data_binary(12, 3)),
               "'data_c' must be created by data_continuous\\(\\).")
  expect_error(boundary_continuous(des, 12, "a", 25, dc),
               "'mean_t' must be numeric.")
  expect_error(boundary_continuous(ref_design_cont("C1", method = "MM"), 4,
                                   NULL, 25, dc),
               "method = \"MM\" requires the distributions of both group means")
})
