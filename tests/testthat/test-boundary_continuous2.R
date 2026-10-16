kf1_cov <- matrix(c(6400, 190, 190, 36), 2L)

test_that("boundary_continuous2() conditions on endpoint 1", {
  des <- ref_design_cont2("KF1")
  s <- boundary_continuous2(des, 10, list(c(20, 5), seq(-5, 5, by = 1)),
                            kf1_cov, NULL)
  expect_identical(s$cond$mean_t1, c(5, 20))
  expect_identical(s$cond$effect1, c(15, 30))
  expect_identical(s$range, c(-5, 5))
  expect_length(s$values[[1L]], 15L)
  p <- s$prob(1L, 2)
  dec <- qdm_decide(des, ref_data_cont2$KF1$t)
  expect_equal(c(p$go, p$nogo), unname(c(dec$prob_go, dec$prob_nogo)))
  expect_equal(s$derived(1L, 2)$effect, 3)
  expect_identical(s$axes$x, "mean_t1")
  expect_identical(s$axes$y, "cut")
  expect_match(s$titles[["effect"]], "fixed control mean -1")
})

test_that("boundary_continuous2() uses the control data and distributions", {
  des <- ref_design_cont2("K1")
  dc <- ref_data_cont2$K1$c
  s <- boundary_continuous2(des, ref_design_cont2("K1")$n_t,
                            list(0, c(-1, 1)), ref_data_cont2$K1$t$cov, dc)
  expect_identical(s$axes$x, "effect1")
  expect_identical(s$axes$offset, dc$mean[2L])
  expect_equal(s$cond$effect1, -dc$mean[1L])
  # The probabilities are those of qdm_decide() with the same data
  p <- s$prob(1L, 1)
  dec <- qdm_decide(des, data_continuous2(12, c(0, 1),
                                          ref_data_cont2$K1$t$cov), dc)
  expect_equal(c(p$go, p$nogo), unname(c(dec$prob_go, dec$prob_nogo)))
  # Hypothetical and external controls
  s4 <- boundary_continuous2(ref_design_cont2("K4"), 15, list(0, 0),
                             kf1_cov, NULL)
  expect_equal(s4$derived(1L, 2)$effect, 3)
  s5 <- boundary_continuous2(ref_design_cont2("K5"), 12, list(0, 0),
                             kf1_cov, NULL)
  expect_match(s5$titles[["effect"]], "posterior mean of the control mean")
  expect_identical(s5$range, c(0, 0))
})

test_that("boundary_continuous2() validates its arguments", {
  des <- ref_design_cont2("KF1")
  expect_error(boundary_continuous2(des, 10, c(1, 2), kf1_cov, NULL),
               "'mean_t' must be a list of two")
  expect_error(boundary_continuous2(des, 10, list(1, 2), NULL, NULL),
               "'cov_t' must be given")
  expect_error(boundary_continuous2(ref_design_cont2("K1"), 30, list(1, 2),
                                    kf1_cov, data_continuous(30, 1, 2)),
               "'data_c' must be created by data_continuous2\\(\\).")
})
