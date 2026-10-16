# Expected colors: Python_reference/ref_plot_overlay.py (shades of the
# decision colors with weights 0.2, 0.4, 0.6, 0.8 and 1 for the five steps)
overlay_x <- function(miss = 0) {
  eff <- expand.grid(th1 = c(0, 0.2), th2 = c(0, 0.2, 0.4))
  data.frame(pi_t1 = 0.3 + eff$th1, pi_t2 = 0.3 + eff$th2, pi_c1 = 0.3,
             pi_c2 = 0.3,
             Go = c(0.62, 0.10, 0.05, NA, 0.72, 0.2),
             Gray = c(0.30, 0.45, 0.15, NA, 0.2, 0.3),
             NoGo = c(0.08, 0.45, 0.80, NA, 0.08, 0.499),
             Miss = miss)
}
overlay_cols <- decision_colors(col_go = "#658D1B", col_nogo = "#D91E49",
                                col_gray = "#939597", col_miss = "#E69F00")
# Rows of a built layer in the order of the scenarios (endpoint 1 varies
# fastest)
by_cell <- function(d) d[order(d$y, d$x), , drop = FALSE]

test_that("plot_oc_overlay() colors each cell by the decision with the highest probability", {
  p <- plot_oc_overlay(overlay_x(), ref_design_bin2("B2"), overlay_cols,
                       plot_style())
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  # Go step 3, Gray step 1, NoGo step 5, not evaluable, Go step 4, NoGo
  # step 1
  expect_identical(by_cell(b$data[[1]])$fill,
                   c("#A3BB76", "#E9EAEA", "#D91E49", "#FFFFFF", "#84A449",
                     "#F7D2DB"))
  expect_identical(by_cell(b$data[[3]])$label,
                   c("0.62", "0.45", "0.80", "NE", "0.72", "0.50"))
  # The legend layer is invisible
  expect_identical(unique(b$data[[2]]$alpha), 0)
  fill <- plot_scales(p)$get_scales("fill")
  expect_identical(fill$name, "Decision")
  expect_identical(unname(fill$get_labels()), c("Go", "Gray", "NoGo"))
  colour <- plot_scales(p)$get_scales("colour")
  expect_identical(colour$name, "Probability")
  expect_identical(colour$get_limits(),
                   c("< 0.50", "0.50-0.60", "0.60-0.70", "0.70-0.80",
                     ">= 0.80"))
  # The steps of the legend are shades of the gray #808080
  expect_identical(unname(colour$map(colour$get_limits())),
                   c("#E6E6E6", "#CCCCCC", "#B3B3B3", "#999999", "#808080"))
  # Every decision has a key, although only NoGo reaches the highest step
  layers <- if (is.null(attr(p, "layers"))) p$layers else attr(p, "layers")
  expect_true(isTRUE(layers[[1L]]$show.legend[["fill"]]))
  expect_identical(plot_labels(p)$x, "Treatment effect on endpoint 1")
})

test_that("plot_oc_overlay() leaves out the probabilities and adds Miss", {
  p <- plot_oc_overlay(overlay_x(), ref_design_bin2("B2"), overlay_cols,
                       plot_style(), cell_values = FALSE)
  expect_identical(by_cell(ggplot2::ggplot_build(p)$data[[3]])$label,
                   c("", "", "", "NE", "", ""))
  # Miss with probability 0.55 in the second scenario (Go 0.10, Gray 0.45,
  # NoGo 0.45): Miss step 2; the first scenario stays Go step 3
  p2 <- plot_oc_overlay(overlay_x(c(0, 0.55, 0, NA, 0, 0)),
                        ref_design_bin2("B2"), overlay_cols, plot_style())
  b2 <- ggplot2::ggplot_build(p2)
  expect_identical(by_cell(b2$data[[1]])$fill[1:2], c("#A3BB76", "#F5D999"))
  expect_identical(unname(plot_scales(p2)$get_scales("fill")$get_labels()),
                   c("Go", "Gray", "NoGo", "Miss"))
  # Other colors: with black, the shade of weight w is the gray level
  # 255 (1 - w)
  p3 <- plot_oc_overlay(overlay_x(), ref_design_bin2("B2"),
                        decision_colors(col_go = "black",
                                        col_nogo = "darkgreen",
                                        col_gray = "gray", col_miss = "red"),
                        plot_style(label_size = 4))
  b3 <- ggplot2::ggplot_build(p3)
  expect_identical(by_cell(b3$data[[1]])$fill[c(1L, 5L)],
                   c("#666666", "#333333"))
  expect_identical(unique(b3$data[[3]]$size), 4)
})
