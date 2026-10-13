test_that("plot.qdm_oc() draws the decision probabilities", {
  oc <- qdm_oc(ref_design("D1"),
               data.frame(pi_t = c(0.6, 0.2, 0.4), pi_c = 0.4))
  p <- plot(oc)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  # Layers: threshold lines, their labels, curves, points
  expect_equal(sort(unique(b$data[[1]]$xintercept)), c(0.1, 0.2))
  lines <- b$data[[3]]
  expect_equal(nrow(lines), 9L)
  expect_equal(sort(unique(lines$x)), c(-0.2, 0, 0.2))
  # Curves for Go, NoGo and Gray, in increasing order of the effect
  expect_identical(sort(unique(lines$group)), 1:3)
  expect_equal(lines$y[lines$group == 1], oc$Go[order(oc$pi_t)])
  expect_equal(lines$y[lines$group == 3], oc$Gray[order(oc$pi_t)])
  expect_identical(unique(lines$colour[lines$group == 1]), "#658D1B")
})

test_that("plot.qdm_oc() adds Miss decisions and shifted thresholds", {
  expect_warning(oc <- qdm_oc(ref_design("D6"), scen_c), "Miss")
  lines <- ggplot2::ggplot_build(plot(oc))$data[[3]]
  expect_identical(sort(unique(lines$group)), 1:4)
  expect_equal(lines$y[lines$group == 4], oc$Miss)
  expect_identical(unique(lines$colour[lines$group == 4]), "#E69F00")
  oc3 <- qdm_oc(ref_design("D3"), scen_u)
  b3 <- ggplot2::ggplot_build(plot(oc3))
  expect_equal(sort(unique(b3$data[[1]]$xintercept)), 5 / 14 + c(0.1, 0.2))
  expect_equal(sort(unique(b3$data[[3]]$x)), scen_u$pi_t)
})

test_that("plot.qdm_oc() needs the design attribute", {
  oc <- qdm_oc(ref_design("D1"), scen_c)
  attr(oc, "design") <- NULL
  expect_error(plot(oc), "'x' must be an unmodified result of qdm_oc().")
})

test_that("plot.qdm_oc() draws continuous operating characteristics", {
  oc <- qdm_oc(ref_design_cont("C1", "MM"), scen_cont_c, nsim = 200,
               seed = 1)
  b <- ggplot2::ggplot_build(plot(oc))
  expect_equal(sort(unique(b$data[[1]]$xintercept)), c(15, 30))
  expect_equal(sort(unique(b$data[[3]]$x)), c(0, 20, 40))
  oc3 <- qdm_oc(ref_design_cont("C3", "MM"), scen_cont_u, nsim = 200,
                seed = 1)
  b3 <- ggplot2::ggplot_build(plot(oc3))
  expect_equal(sort(unique(b3$data[[1]]$xintercept)), c(5, 20))
  expect_equal(sort(unique(b3$data[[3]]$x)), scen_cont_u$mu_t)
})
