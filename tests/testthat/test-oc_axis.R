test_that("oc_axis() uses the treatment effect with a concurrent control", {
  ax <- oc_axis(ref_design("D1"), scen_c)
  expect_equal(ax$x, scen_c$pi_t - scen_c$pi_c)
  expect_equal(ax$lines, c(0.1, 0.2))
  expect_identical(ax$labels, c("theta[MAV]", "theta[TV]"))
  expect_null(ax$caption)
  ax5 <- oc_axis(ref_design("D5"), scen_c)
  expect_equal(ax5$lines, 0.05)
  expect_identical(ax5$labels, "theta[NULL]")
})

test_that("oc_axis() shifts the thresholds without a concurrent control", {
  # Control distribution Beta(5, 9) in D3 and Beta(4, 8) in D4
  ax <- oc_axis(ref_design("D3"), scen_u)
  expect_equal(ax$x, scen_u$pi_t)
  expect_equal(ax$lines, 5 / 14 + c(0.1, 0.2))
  expect_match(ax$caption, "shifted by the mean of the control distribution")
  expect_match(ax$caption, "0.357")
  expect_equal(oc_axis(ref_design("D4"), scen_u)$lines, 1 / 3 + c(0.1, 0.2))
})

test_that("oc_axis() uses the mean difference for continuous designs", {
  ax <- oc_axis(ref_design_cont("C1"), scen_cont_c)
  expect_equal(ax$x, scen_cont_c$mu_t - scen_cont_c$mu_c)
  expect_equal(ax$xlab, expression(theta == mu[t] - mu[c]))
  expect_equal(ax$lines, c(15, 30))
  expect_null(ax$caption)
  ax5 <- oc_axis(ref_design_cont("C5"), scen_cont_c)
  expect_equal(ax5$lines, 10)
})

test_that("oc_axis() shifts the thresholds by the control location", {
  # Hypothetical control located at -10 (C3) and at 0 (C6); external control
  # data with mean -12 (C4)
  ax <- oc_axis(ref_design_cont("C3"), scen_cont_u)
  expect_equal(ax$x, scen_cont_u$mu_t)
  expect_equal(ax$xlab, expression(mu[t]))
  expect_equal(ax$lines, -10 + c(15, 30))
  expect_match(ax$caption,
               "shifted by the location of the control distribution \\(-10\\)")
  expect_equal(oc_axis(ref_design_cont("C4"), scen_cont_u)$lines,
               -12 + c(15, 30))
  expect_equal(oc_axis(ref_design_cont("C6"), scen_cont_u)$lines, 5)
})

test_that("oc_axis() shifts the thresholds by a fixed control value", {
  ax <- oc_axis(ref_design("D8"), scen_u)
  expect_equal(ax$x, scen_u$pi_t)
  expect_equal(ax$lines, 0.2 + c(0.05, 0.2))
  expect_match(ax$caption, "shifted by the fixed control value \\(0.2\\)")
})

test_that("oc_axis() uses milestone survival for time-to-event designs", {
  sc <- data.frame(median_t = c(6, 12), median_c = 6)
  ax <- oc_axis(ref_design_tte("T2"), sc)
  s <- function(m) exp(-log(2) * 3 / m)
  expect_equal(ax$x, s(c(6, 12)) - s(6))
  expect_equal(ax$xlab, expression(theta == S[t](t[0]) - S[c](t[0])))
  expect_null(ax$caption)
  ax1 <- oc_axis(ref_design_tte("T1"), sc[, 1, drop = FALSE])
  expect_equal(ax1$x, s(c(6, 12)))
  expect_equal(ax1$xlab, expression(S[t](t[0])))
  expect_equal(ax1$lines, 0.72 + c(0, 0.15))
  ax5 <- oc_axis(ref_design_tte("T5"), sc[, 1, drop = FALSE])
  expect_match(ax5$caption, "shifted by the mean of the control distribution")
})
