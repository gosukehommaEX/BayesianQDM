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

test_that("region_prob_tte() reverses the regions for direction = \"less\"", {
  # Thresholds exchanged: the Go region of direction = "less" is the NoGo
  # region of the reference criterion and conversely (designs T1 and T2)
  less <- criterion_posterior(0, 0.15, 0.6, 0.8, direction = "less")
  d1 <- ref_design_tte("T1")
  d1$criterion <- less
  p1 <- region_prob_tte(d1, list(shape = 6.1, rate = 151), NULL)
  expect_equal(c(p1$go, p1$nogo), unname(rev(ref_prob_tte$T1)),
               tolerance = 1e-10)
  d2 <- ref_design_tte("T2")
  d2$criterion <- less
  p2 <- region_prob_tte(d2, list(shape = 5.1, rate = 41),
                        list(shape = 9.1, rate = 36))
  expect_equal(c(p2$go, p2$nogo), unname(rev(ref_prob_tte$T2)),
               tolerance = 1e-10)
  # Predictive: P(future effect < NULL) is the complement of the reference
  d3 <- ref_design_tte("T3", nMC = 500L)
  d3l <- d3
  d3l$criterion$direction <- "less"
  set.seed(5)
  a <- region_prob_tte(d3, list(shape = 6.1, rate = 71), NULL)
  set.seed(5)
  b <- region_prob_tte(d3l, list(shape = 6.1, rate = 71), NULL)
  expect_identical(b$go, a$nogo)
  expect_identical(b$nogo, a$go)
})

test_that("region_prob_tte() computes probabilities of the hazard ratio", {
  # Python_reference/ref_hazard_ratio.py (integration over the control
  # hazard)
  p1 <- region_prob_tte(ref_design_hr("H1"), list(shape = 12.1, rate = 201),
                        list(shape = 20.1, rate = 181))
  expect_equal(c(p1$go, p1$nogo), unname(ref_prob_hr$H1), tolerance = 1e-10)
  p2 <- region_prob_tte(ref_design_hr("H2"), list(shape = 12.1, rate = 201),
                        NULL)
  expect_equal(c(p2$go, p2$nogo), unname(ref_prob_hr$H2), tolerance = 1e-10)
  # Direction "greater": Go if P(HR > 1.3) and NoGo if P(HR <= 1)
  d <- ref_design_hr("H1")
  d$criterion <- criterion_posterior(1.3, 1, 0.6, 0.7)
  p3 <- region_prob_tte(d, list(shape = 12.1, rate = 201),
                        list(shape = 20.1, rate = 181))
  expect_equal(c(p3$go, p3$nogo), c(0.00684327358731, 0.9581768777711),
               tolerance = 1e-10)
})

test_that("region_prob_tte() estimates predictive probabilities of the hazard ratio", {
  # Complete follow-up of the future study: beta distributions of the
  # future total observed times (Python_reference/ref_hazard_ratio.py)
  n_mc <- 20000L
  set.seed(6)
  p1 <- region_prob_tte(ref_design_hr("P1", nMC = n_mc),
                        list(shape = 12.1, rate = 201), NULL)
  ref1 <- ref_prob_hr$P1[["go"]]
  expect_lt(abs(p1$go - ref1), 4.5 * sqrt(ref1 * (1 - ref1) / n_mc))
  expect_identical(p1$nogo, 1 - p1$go)
  set.seed(7)
  p2 <- region_prob_tte(ref_design_hr("P2", nMC = n_mc),
                        list(shape = 12.1, rate = 201),
                        list(shape = 20.1, rate = 181))
  ref2 <- ref_prob_hr$P2[["go"]]
  expect_lt(abs(p2$go - ref2), 4.5 * sqrt(ref2 * (1 - ref2) / n_mc))
})
