test_that("plot_style() returns the default style", {
  s <- plot_style()
  expect_identical(s$linewidth, 1)
  expect_identical(s$point_size, 1.5)
  expect_equal(s$label_size, 0.8 * 11 / ggplot2::.pt)
  expect_identical(s$theme[["legend.position"]], "bottom")
  expect_identical(element_size(s$theme[["text"]]), 11)
})

test_that("plot_style() takes the given values", {
  s <- plot_style(linewidth = 0.5, point_size = 3, label_size = 4,
                  base_size = 16, legend_position = "none")
  expect_identical(s$linewidth, 0.5)
  expect_identical(s$point_size, 3)
  expect_identical(s$label_size, 4)
  expect_identical(s$theme[["legend.position"]], "none")
  expect_identical(element_size(s$theme[["text"]]), 16)
  # The default label size follows the base size
  expect_equal(plot_style(base_size = 20)$label_size, 0.8 * 20 / ggplot2::.pt)
  for (pos in c("top", "left", "right")) {
    expect_identical(plot_style(legend_position = pos)$theme[[
      "legend.position"]], pos)
  }
})

test_that("plot_style() rejects invalid values", {
  expect_error(plot_style(linewidth = -1), "'linewidth' must be at least 0.",
               fixed = TRUE)
  expect_error(plot_style(point_size = "a"), "'point_size' must be numeric.",
               fixed = TRUE)
  expect_error(plot_style(label_size = c(1, 2)),
               "'label_size' must have length 1.", fixed = TRUE)
  expect_error(plot_style(base_size = 0), "'base_size' must be greater than 0.",
               fixed = TRUE)
  msg <- "'legend_position' must be one of"
  expect_error(plot_style(legend_position = "inside"), msg, fixed = TRUE)
  expect_error(plot_style(legend_position = c("top", "left")), msg,
               fixed = TRUE)
  expect_error(plot_style(legend_position = 1), msg, fixed = TRUE)
  expect_error(plot_style(legend_position = NA_character_), msg, fixed = TRUE)
})
