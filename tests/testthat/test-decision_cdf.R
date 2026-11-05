# References: Python_reference/ref_plot_decision.out (D1, D8, C1, C7, T2,
# H1) and helper-designs-tte.R (ref_prob_tte, ref_prob_hr). The
# distribution function of each decision is evaluated at the target value
# and at the minimum acceptable value.
cdf_at <- function(dec, q) {
  decision_cdf(dec$design, dec$posterior_t, dec$posterior_c)(q)
}

test_that("decision_cdf() gives the distribution of a binary effect", {
  p <- cdf_at(qdm_decide(ref_design("D1"), data_binary(12, 9),
                         data_binary(12, 4)), c(0.2, 0.1, -1, 1))
  expect_equal(1 - p[1L], 0.820732075666572, tolerance = 1e-7)
  expect_equal(p[2L], 0.0729266458410321, tolerance = 1e-7)
  expect_equal(p[3:4], c(0, 1))
  p <- cdf_at(qdm_decide(ref_design("D8"), data_binary(20, 9)), c(0.2, 0.05))
  expect_equal(1 - p[1L], 0.691441812114726, tolerance = 1e-9)
  expect_equal(p[2L], 0.0206296046344505, tolerance = 1e-9)
})

test_that("decision_cdf() gives the distribution of a continuous effect", {
  p <- cdf_at(qdm_decide(ref_design_cont("C1"), ref_data_cont$C1$t,
                         ref_data_cont$C1$c), c(30, 15))
  expect_equal(1 - p[1L], 0.5, tolerance = 1e-7)
  expect_equal(p[2L], 0.0878402872646904, tolerance = 1e-7)
  p <- cdf_at(qdm_decide(ref_design_cont("C7"), ref_data_cont$C7$t), 30)
  expect_equal(1 - p, 0.797515237746117, tolerance = 1e-9)
  # Hypothetical control: the probabilities of the decision (numerical
  # integration)
  dec <- qdm_decide(ref_design_cont("C3"), ref_data_cont$C3$t)
  p <- cdf_at(dec, c(30, 15))
  expect_equal(1 - p[1L], unname(dec$prob_go), tolerance = 1e-9)
  expect_equal(p[2L], unname(dec$prob_nogo), tolerance = 1e-9)
})

test_that("decision_cdf() gives the distribution of a time-to-event effect", {
  p <- cdf_at(qdm_decide(ref_design_tte("T2"), ref_data_tte$T2$t,
                         ref_data_tte$T2$c), c(0.15, 0))
  expect_equal(1 - p[1L], ref_prob_tte$T2[["go"]], tolerance = 1e-7)
  expect_equal(p[2L], ref_prob_tte$T2[["nogo"]], tolerance = 1e-7)
  p <- cdf_at(qdm_decide(ref_design_tte("T1"), ref_data_tte$T1$t), 0.15)
  expect_equal(1 - p, ref_prob_tte$T1[["go"]], tolerance = 1e-7)
  # Hazard ratio: Go below the target value 0.7
  p <- cdf_at(qdm_decide(ref_design_hr("H1"), ref_data_hr$t, ref_data_hr$c),
              c(0.7, 1))
  expect_equal(p[1L], ref_prob_hr$H1[["go"]], tolerance = 1e-7)
  expect_equal(1 - p[2L], ref_prob_hr$H1[["nogo"]], tolerance = 1e-7)
  p <- cdf_at(qdm_decide(ref_design_hr("H2"), ref_data_hr$t), c(0.7, -1))
  expect_equal(p[1L], ref_prob_hr$H2[["go"]], tolerance = 1e-7)
  expect_identical(p[2L], 0)
})
