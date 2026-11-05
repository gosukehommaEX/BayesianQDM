test_that("plot_decision_bars() draws the predictive probabilities", {
  dec <- qdm_decide(ref_design("D5"), data_binary(10, 8), data_binary(10, 3))
  p <- plot(dec)
  b <- ggplot2::ggplot_build(p)
  # Layers: bars, thresholds, values
  bars <- b$data[[1]]
  probs <- unname(c(dec$prob_go, dec$prob_nogo))
  expect_equal(bars$ymax, probs)
  expect_identical(bars$fill, c("#658D1B", "#D91E49"))
  expect_equal(b$data[[2]]$y, c(0.9, 0.7))
  expect_identical(b$data[[3]]$label, formatC(probs, format = "f",
                                              digits = 3L))
  expect_identical(plot_labels(p)$title,
                   sprintf("Decision: %s", dec$decision))
  expect_identical(plot_labels(p)$subtitle, decision_subtitle(dec))
})

test_that("plot_decision_bars() is used for a predictive criterion only", {
  dec <- qdm_decide(ref_design_tte("T3"), ref_data_tte$T3$t, seed = 1)
  b <- ggplot2::ggplot_build(plot(dec))
  expect_equal(b$data[[1]]$ymax, unname(c(dec$prob_go, dec$prob_nogo)))
  expect_equal(b$data[[2]]$y, c(0.8, 0.3))
})
