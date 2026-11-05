test_that("plot.qdm_comparison() draws one curve per design and decision", {
  cmp <- qdm_compare(list(a = ref_design("D1"), b = ref_design("D2")), scen_c)
  p <- plot(cmp)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  # Layers: thresholds, their labels, curves, points (the thresholds are
  # drawn in each panel)
  expect_equal(unique(b$data[[1]]$xintercept), c(0.1, 0.2))
  expect_identical(nrow(b$data[[1]]), 2L * 3L)
  expect_identical(nrow(b$data[[2]]), 2L * 3L)
  expect_equal(unique(b$data[[2]]$x), c(0.1, 0.2))
  expect_setequal(b$data[[2]]$label, c("theta[MAV]", "theta[TV]"))
  # The second label is lower, so that the labels do not overlap
  expect_equal(unique(b$data[[2]]$y), c(1, 0.93))
  lines <- b$data[[3]]
  expect_identical(nrow(lines), 2L * 3L * 4L)
  cols <- grDevices::hcl.colors(2, "Dark 3")
  expect_setequal(unique(lines$colour), cols)
  # Panels Go, NoGo and Gray
  expect_identical(nrow(b$layout$layout), 3L)
  go_b <- lines[lines$PANEL == 1 & lines$colour == cols[2], ]
  expect_equal(go_b$x, scen_c$pi_t - scen_c$pi_c)
  expect_equal(go_b$y, cmp$Go[cmp$design == "b"])
  gray_a <- lines[lines$PANEL == 3 & lines$colour == cols[1], ]
  expect_equal(gray_a$y, cmp$Gray[cmp$design == "a"])
})

test_that("plot.qdm_comparison() takes the decisions, colors and style", {
  cmp <- qdm_compare(list(a = ref_design("D1"), b = ref_design("D2")), scen_c)
  p <- plot(cmp, decisions = c("NoGo", "Go"), col_design = c("black", "blue"),
            linewidth = 0.5, point_size = 2, base_size = 9,
            legend_position = "top")
  b <- ggplot2::ggplot_build(p)
  expect_identical(nrow(b$layout$layout), 2L)
  lines <- b$data[[3]]
  expect_setequal(unique(lines$colour), c("black", "blue"))
  expect_identical(unique(lines$linewidth), 0.5)
  expect_identical(unique(b$data[[4]]$size), 2)
  # The first panel is Go (the decisions in their usual order)
  expect_equal(lines$y[lines$PANEL == 1 & lines$colour == "black"],
               cmp$Go[cmp$design == "a"])
  th <- plot_theme(p)
  expect_identical(th[["legend.position"]], "top")
  expect_identical(element_size(th[["text"]]), 9)
  expect_error(plot(cmp, decisions = "Win"),
               "'decisions' must contain some of", fixed = TRUE)
  expect_error(plot(cmp, decisions = character(0)),
               "'decisions' must contain some of", fixed = TRUE)
  expect_error(plot(cmp, col_design = "black"),
               "'col_design' must give one color for each of the 2 designs.",
               fixed = TRUE)
  expect_error(plot(cmp, col_design = c("black", "nocolor")),
               "'col_design' must be a single color name", fixed = TRUE)
  expect_error(plot(cmp, cell_values = FALSE),
               "'cell_values' applies only to the heat maps of two endpoints.",
               fixed = TRUE)
  expect_warning(plot(cmp, overlay = TRUE),
                 "Ignoring unknown argument(s) of plot.qdm_comparison(): overlay.",
                 fixed = TRUE)
  one <- qdm_compare(list(ref_design("D1")), scen_c)
  expect_error(plot(one, col_design = c("red", "blue")),
               "'col_design' must give one color for each of the 1 design.",
               fixed = TRUE)
})

test_that("plot.qdm_comparison() shows Miss and omits different thresholds", {
  cmp <- suppressWarnings(
    qdm_compare(list(a = ref_design("D1"), b = ref_design("D6")), scen_c)
  )
  expect_true(any(cmp$Miss > 0))
  expect_identical(nrow(ggplot2::ggplot_build(plot(cmp))$layout$layout), 4L)
  # D7 has other thresholds, so no threshold is drawn
  cmp2 <- qdm_compare(list(a = ref_design("D1"), b = ref_design("D7")), scen_c)
  b <- ggplot2::ggplot_build(plot(cmp2))
  expect_length(b$data, 2L)
})

test_that("plot.qdm_comparison() needs the same horizontal axis", {
  cmp <- qdm_compare(list(a = ref_design("D1"), b = ref_design("D3")), scen_c)
  expect_error(plot(cmp),
               "The designs do not have the same horizontal axis",
               fixed = TRUE)
  dup <- qdm_compare(list(a = ref_design("D1")),
                     data.frame(pi_t = c(0.3, 0.4), pi_c = c(0.3, 0.4)))
  expect_error(plot(dup),
               "The scenarios must have distinct values on the horizontal axis",
               fixed = TRUE)
  bad <- cmp
  attr(bad, "designs") <- NULL
  expect_error(plot(bad), "'x' must be an unmodified result of qdm_compare().",
               fixed = TRUE)
})
