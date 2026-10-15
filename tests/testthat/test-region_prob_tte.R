test_that("region_prob_tte() matches the references (posterior)", {
  # Python_reference/ref_fixed_tte_qdm.py, designs T1 and T2
  p1 <- region_prob_tte(ref_design_tte("T1"), list(shape = 6.1, rate = 151),
                        NULL)
  expect_equal(c(p1$go, p1$nogo), unname(ref_prob_tte$T1), tolerance = 1e-10)
  p2 <- region_prob_tte(ref_design_tte("T2"),
                        list(shape = c(5.1, 5.1), rate = c(41, 41)),
                        list(shape = 9.1, rate = 36))
  expect_equal(p2$go, rep(ref_prob_tte$T2[["go"]], 2), tolerance = 1e-10)
  expect_equal(p2$nogo, rep(ref_prob_tte$T2[["nogo"]], 2), tolerance = 1e-10)
})

test_that("region_prob_tte() estimates predictive probabilities", {
  # Complete follow-up of the future study: closed forms
  # (Python_reference/ref_fixed_tte_qdm.py, designs T3 and T4)
  n_mc <- 20000L
  set.seed(1)
  p3 <- region_prob_tte(ref_design_tte("T3", nMC = n_mc),
                        list(shape = 6.1, rate = 71), NULL)
  ref3 <- ref_prob_tte$T3[["go"]]
  expect_lt(abs(p3$go - ref3), 4.5 * sqrt(ref3 * (1 - ref3) / n_mc))
  expect_identical(p3$nogo, 1 - p3$go)
  set.seed(2)
  p4 <- region_prob_tte(ref_design_tte("T4", nMC = n_mc),
                        list(shape = 6.1, rate = 71),
                        list(shape = 9.1, rate = 61))
  ref4 <- ref_prob_tte$T4[["go"]]
  expect_lt(abs(p4$go - ref4), 4.5 * sqrt(ref4 * (1 - ref4) / n_mc))
  # Several distributions; the control parameters are recycled
  set.seed(3)
  p5 <- region_prob_tte(ref_design_tte("T4", nMC = 200L),
                        list(shape = c(6.1, 30.1), rate = c(71, 71)),
                        list(shape = 9.1, rate = 61))
  expect_length(p5$go, 2L)
  expect_gt(p5$go[1], p5$go[2])
})
