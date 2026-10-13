test_that("plot.qdm_calibration() draws the calibration curves", {
  cal <- run_cal("D1")
  p <- plot(cal)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  # Layers: targets, chosen thresholds, curves
  expect_equal(sort(b$data[[1]]$yintercept), c(0.05, 0.2))
  expect_equal(sort(b$data[[2]]$xintercept), c(0.41, 0.7))
  steps <- b$data[[3]]
  expect_equal(nrow(steps), 2L * 99L)
  expect_equal(steps$y[steps$group == 1], cal$curve$prob_go)
  expect_equal(steps$y[steps$group == 2], cal$curve$prob_nogo)
  expect_identical(unique(steps$colour[steps$group == 2]), "#D91E49")
})

test_that("plot.qdm_calibration() omits thresholds that were not found", {
  a <- cal_args$D1
  cal <- suppressWarnings(
    qdm_calibrate(ref_design("D1"), a$scenario_go, a$scenario_nogo,
                  target_go = 1e-6, target_nogo = 0.2)
  )
  b <- ggplot2::ggplot_build(plot(cal))
  expect_equal(b$data[[2]]$xintercept, 0.41)
})
