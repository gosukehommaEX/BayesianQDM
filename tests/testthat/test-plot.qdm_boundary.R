test_that("plot.qdm_boundary() draws the decision intervals", {
  b <- qdm_boundary(ref_design("D8"), n_t = 10:20)
  p <- plot(b)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  rects <- built$data[[1]]
  expect_identical(nrow(rects), nrow(b))
  # Rectangles of whole numbers: n_t +- 0.5 and the responders +- 0.5
  expect_equal(rects$xmin, b$n_t - 0.5)
  expect_equal(rects$ymax, b$to + 0.5)
  cols <- c(Go = "#658D1B", NoGo = "#D91E49", Gray = "#939597",
            Miss = "#E69F00")
  expect_identical(unname(rects$fill), unname(cols[as.character(b$decision)]))
  labs <- plot_labels(p)
  expect_identical(labs$x, "Treatment group size")
  expect_identical(labs$y, "Responders in the treatment group")
  expect_identical(labs$fill, "Decision")
  # The boundary values: the smallest Go count and the largest NoGo count
  # for every sample size (11 values of n_t)
  text <- built$data[[2]]
  go <- b[b$decision == "Go", ]
  nogo <- b[b$decision == "NoGo", ]
  expect_identical(nrow(text), 22L)
  expect_setequal(paste(text$x, text$label),
                  c(paste(go$n_t, go$from), paste(nogo$n_t, nogo$to)))
})

test_that("plot.qdm_boundary() labels some or none of the conditions", {
  b <- qdm_boundary(ref_design("D8"), n_t = 10:40)
  # More than 16 values: about six of them
  text <- ggplot2::ggplot_build(plot(b))$data[[2]]
  expect_true(all(unique(text$x) %in% c(10, 15, 20, 25, 30, 35, 40)))
  expect_length(ggplot2::ggplot_build(plot(b, labels = FALSE))$data, 1L)
  text2 <- ggplot2::ggplot_build(plot(b, labels = c(12, 30)))$data[[2]]
  expect_setequal(unique(text2$x), c(12, 30))
  expect_error(plot(b, labels = 50),
               "'labels' must be TRUE, FALSE or values of the horizontal")
  expect_error(plot(b, labels = "a"),
               "'labels' must be TRUE, FALSE or values of the horizontal")
})

test_that("plot.qdm_boundary() shows the response rate on a secondary axis", {
  p <- plot(qdm_boundary(ref_design("D1")))
  sec <- plot_scales(p)$get_scales("y")$secondary.axis
  expect_false(is.null(sec))
  expect_identical(sec$name, "Observed treatment response rate")
  labs <- plot_labels(p)
  expect_identical(labs$x, "Responders in the control group (n = 12)")
})

test_that("plot.qdm_boundary() draws observed means against a continuous condition", {
  b <- qdm_boundary(ref_design_cont("C1", method = "MM"),
                    sd_t = c(20, 25, 30), data_c = ref_data_cont$C1$c)
  p <- plot(b)
  built <- ggplot2::ggplot_build(p)
  rects <- built$data[[1]]
  # Midpoints between the values of the standard deviation
  expect_equal(sort(unique(rects$xmin)), c(17.5, 22.5, 27.5))
  # The open ends are drawn to the edges of the panel
  expect_true(any(is.infinite(c(b$effect_from, b$effect_to))))
  text <- built$data[[2]]
  go <- b[b$decision == "Go", ]
  expect_equal(sort(text$y[text$label %in%
                             vapply(go$effect_from, format, "", digits = 3)]),
               sort(go$effect_from))
  labs <- plot_labels(p)
  expect_identical(labs$y, "Observed difference in means")
  expect_match(labs$caption, "^Control data: ")
  # The vertical range covers the boundaries
  lim <- built$layout$panel_params[[1]]$y.range
  expect_true(lim[1L] < min(b$effect_to[b$decision == "NoGo"]))
  expect_true(lim[2L] > max(go$effect_from))
})

test_that("plot.qdm_boundary() handles a single condition and a constant decision", {
  b <- qdm_boundary(ref_design_cont("C7"), sd_t = 20,
                    mean_t = seq(-10, 40, by = 5))
  rects <- ggplot2::ggplot_build(plot(b))$data[[1]]
  expect_equal(unique(rects$xmin), 11.5)
  # Without any finite boundary, the grid gives the range (and there are
  # no boundary values to write)
  b2 <- b
  b2$decision[] <- "Gray"
  b2$from <- -Inf
  b2$to <- Inf
  built <- ggplot2::ggplot_build(plot(b2))
  lim <- built$layout$panel_params[[1]]$y.range
  expect_true(lim[1L] <= -10 && lim[2L] >= 40)
  expect_length(built$data, 1L)
})

test_that("plot.qdm_boundary() needs the axes attribute", {
  b <- qdm_boundary(ref_design("D8"))
  attr(b, "axes") <- NULL
  expect_error(plot(b),
               "'x' must be an unmodified result of qdm_boundary\\(\\).")
})

test_that("plot.qdm_boundary() takes colors, text sizes and legend", {
  b <- qdm_boundary(ref_design("D8"), n_t = 10:20)
  p <- plot(b, col_go = "black", col_nogo = "blue", col_gray = "white",
            base_size = 8, legend_position = "none", label_size = 2)
  built <- ggplot2::ggplot_build(p)
  cols <- c(Go = "black", NoGo = "blue", Gray = "white", Miss = "#E69F00")
  expect_identical(unname(built$data[[1]]$fill),
                   unname(cols[as.character(b$decision)]))
  expect_identical(unique(built$data[[2]]$size), 2)
  th <- plot_theme(p)
  expect_identical(th[["legend.position"]], "none")
  expect_identical(element_size(th[["text"]]), 8)
  expect_error(plot(b, col_gray = "grey500"),
               "'col_gray' must be a single color name", fixed = TRUE)
  expect_warning(
    plot(b, linewidth = 2),
    "Ignoring unknown argument(s) of plot.qdm_boundary(): linewidth.",
    fixed = TRUE
  )
})
