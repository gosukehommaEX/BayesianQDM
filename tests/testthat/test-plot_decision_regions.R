test_that("plot_decision_regions() draws the regions of two endpoints", {
  dec <- qdm_decide(ref_design_cont2("K1"), ref_data_cont2$K1$t,
                    ref_data_cont2$K1$c)
  p <- plot(dec)
  b <- ggplot2::ggplot_build(p)
  tiles <- b$data[[1]]
  text <- b$data[[2]]
  expect_identical(nrow(tiles), 9L)
  prob <- function(r) {
    formatC(unname(dec$prob_regions[r]), format = "f", digits = 3L)
  }
  # Region 1 (both effects beyond their target values) at the top right,
  # region 9 at the bottom left and region 3 at the bottom right
  at <- function(d, x, y) d[abs(d$x - x) < 1e-9 & abs(d$y - y) < 1e-9, ]
  expect_identical(at(text, 3, 3)$label, paste0("R1\n", prob("R1")))
  expect_identical(at(text, 1, 1)$label, paste0("R9\n", prob("R9")))
  expect_identical(at(text, 3, 1)$label, paste0("R3\n", prob("R3")))
  expect_identical(at(text, 1, 3)$label, paste0("R7\n", prob("R7")))
  expect_setequal(text$label, paste0("R", 1:9, "\n", prob(paste0("R", 1:9))))
  # The Go region (R1) and the NoGo region (R9) have their colors
  shade <- function(r, col) {
    shade_color(0.15 + 0.85 * unname(dec$prob_regions[r]), col)
  }
  expect_identical(at(tiles, 3, 3)$fill, shade("R1", "#658D1B"))
  expect_identical(at(tiles, 1, 1)$fill, shade("R9", "#D91E49"))
  expect_identical(at(tiles, 2, 2)$fill, shade("R5", "#939597"))
  # Axis labels with the thresholds, the favorable side last
  sc <- plot_scales(p)
  expect_identical(sc$get_scales("x")$labels,
                   c("theta1 <= 5", "5 < theta1 <= 10", "theta1 > 10"))
  expect_identical(sc$get_scales("y")$labels,
                   c("theta2 <= 1", "1 < theta2 <= 2", "theta2 > 2"))
  expect_identical(plot_labels(p)$subtitle, decision_subtitle(dec))
})

test_that("plot_decision_regions() covers marginal and predictive criteria", {
  # Marginal rule: no region has the color of a decision. The tile at x and
  # y is in column 4 - x and row 4 - y of qdm_regions()
  dec <- qdm_decide(ref_design_cont2("K5"), ref_data_cont2$K5$t)
  tiles <- ggplot2::ggplot_build(plot(dec))$data[[1]]
  expect_identical(nrow(tiles), 9L)
  region <- (3 - tiles$x) * 3 + (4 - tiles$y)
  expect_identical(tiles$fill, shade_color(
    0.15 + 0.85 * unname(dec$prob_regions[paste0("R", region)]), "#939597"
  ))
  # Predictive criterion: four regions
  pred <- qdm_decide(ref_design_cont2("K3"), ref_data_cont2$K3$t,
                     ref_data_cont2$K3$c)
  p <- plot(pred)
  b <- ggplot2::ggplot_build(p)
  expect_identical(nrow(b$data[[1]]), 4L)
  expect_identical(plot_scales(p)$get_scales("x")$labels,
                   c("theta1 <= 1.984", "theta1 > 1.984"))
})
