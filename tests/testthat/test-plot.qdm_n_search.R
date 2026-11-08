test_that("plot.qdm_n_search() draws the probabilities of the candidates", {
  res <- find_n_d1(n_t = seq(6, 20, by = 2))
  p <- plot(res)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  # Layers: targets, chosen sample size, lines, points. The two targets are
  # equal, so they are one gray line
  expect_equal(b$data[[1]]$yintercept, 0.1)
  expect_identical(b$data[[1]]$colour, "gray40")
  expect_equal(b$data[[2]]$xintercept, 14)
  lines <- b$data[[3]]
  expect_equal(lines$x[lines$group == 1], seq(6, 20, by = 2))
  expect_equal(lines$y[lines$group == 1], res$curve$prob_go)
  expect_equal(lines$y[lines$group == 2], res$curve$prob_nogo)
  expect_identical(unique(lines$colour[lines$group == 1]), "#658D1B")
  expect_identical(unique(lines$colour[lines$group == 2]), "#D91E49")
  expect_identical(plot_labels(p)$x,
                   "Number of subjects in the treatment group")
  # The legend names the false Go and the false NoGo, not the arguments
  expect_identical(levels(plot_data(p)$event), false_decision_labels()$legend)
  # The vertical axis runs from 0 to the largest probability
  expect_equal(b$layout$panel_scales_y[[1]]$get_limits(),
               c(0, max(res$curve$prob_go, res$curve$prob_nogo)))
})

test_that("plot.qdm_n_search() draws different targets in their colors", {
  a <- cal_args$D1
  res <- qdm_find_n(ref_design("D1"), a$scenario_go, a$scenario_nogo,
                    target_go = 0.3, target_nogo = 0.1, n_t = c(12, 14))
  b <- ggplot2::ggplot_build(plot(res))
  expect_equal(b$data[[1]]$yintercept, c(0.3, 0.1))
  expect_identical(b$data[[1]]$colour, c("#658D1B", "#D91E49"))
  # The axis contains the larger target
  expect_equal(b$layout$panel_scales_y[[1]]$get_limits(), c(0, 0.3))
})

test_that("plot.qdm_n_search() omits the choice when none was found", {
  none <- suppressWarnings(find_n_d1(n_t = c(6, 8)))
  expect_identical(nrow(ggplot2::ggplot_build(plot(none))$data[[2]]), 0L)
})

test_that("plot.qdm_n_search() takes colors, line width and legend", {
  res <- find_n_d1(n_t = seq(6, 20, by = 2))
  p <- plot(res, col_go = "black", col_nogo = "blue", linewidth = 0.4,
            point_size = 3, base_size = 9, legend_position = "top")
  b <- ggplot2::ggplot_build(p)
  lines <- b$data[[3]]
  expect_identical(unique(lines$colour[lines$group == 1]), "black")
  expect_identical(unique(lines$linewidth), 0.4)
  expect_identical(unique(b$data[[4]]$size), 3)
  th <- plot_theme(p)
  expect_identical(th[["legend.position"]], "top")
  expect_identical(element_size(th[["text"]]), 9)
  expect_error(plot(res, col_go = NA), "'col_go' must be a single color name",
               fixed = TRUE)
  expect_warning(plot(res, overlay = TRUE),
                 "Ignoring unknown argument(s) of plot.qdm_n_search(): overlay.",
                 fixed = TRUE)
})
