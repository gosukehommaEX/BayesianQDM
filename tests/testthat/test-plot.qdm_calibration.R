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

test_that("plot.qdm_calibration() takes colors, line width and legend", {
  cal <- run_cal("D1")
  p <- plot(cal, col_go = "black", col_nogo = "blue", linewidth = 0.4,
            base_size = 9, legend_position = "top")
  b <- ggplot2::ggplot_build(p)
  steps <- b$data[[3]]
  expect_identical(unique(steps$colour[steps$group == 1]), "black")
  expect_identical(unique(steps$colour[steps$group == 2]), "blue")
  expect_identical(unique(steps$linewidth), 0.4)
  # The targets take the colors of their criteria
  expect_setequal(b$data[[1]]$colour, c("black", "blue"))
  th <- plot_theme(p)
  expect_identical(th[["legend.position"]], "top")
  expect_identical(element_size(th[["text"]]), 9)
  expect_identical(unique(ggplot2::ggplot_build(plot(cal))$data[[3]]$linewidth),
                   1)
  expect_error(plot(cal, col_nogo = NA),
               "'col_nogo' must be a single color name", fixed = TRUE)
  expect_error(plot(cal, base_size = -1), "'base_size' must be greater than 0.",
               fixed = TRUE)
  expect_warning(
    plot(cal, point_size = 2),
    "Ignoring unknown argument(s) of plot.qdm_calibration(): point_size.",
    fixed = TRUE
  )
})
