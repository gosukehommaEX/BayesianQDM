test_that("plot.qdm_decision() chooses the plot of the decision", {
  dec <- qdm_decide(ref_design("D1"), data_binary(12, 9), data_binary(12, 4))
  p <- plot(dec)
  expect_s3_class(p, "ggplot")
  expect_identical(plot_labels(p)$y, "Posterior density")
  pred <- qdm_decide(ref_design("D5"), data_binary(10, 8), data_binary(10, 3))
  expect_identical(plot_labels(plot(pred))$y,
                   "Posterior predictive probability")
  two <- qdm_decide(ref_design_cont2("K1"), ref_data_cont2$K1$t,
                    ref_data_cont2$K1$c)
  expect_identical(plot_labels(plot(two))$x, "Endpoint 1")
})

test_that("plot.qdm_decision() takes colors and style", {
  dec <- qdm_decide(ref_design("D1"), data_binary(12, 9), data_binary(12, 4))
  p <- plot(dec, col_go = "black", col_nogo = "blue", linewidth = 0.5,
            label_size = 3, base_size = 9, legend_position = "top")
  b <- ggplot2::ggplot_build(p)
  expect_setequal(unique(b$data[[1]]$fill), c("black", "blue"))
  expect_identical(unique(b$data[[2]]$linewidth), 0.5)
  expect_identical(unique(b$data[[4]]$size), 3)
  th <- plot_theme(p)
  expect_identical(th[["legend.position"]], "top")
  expect_identical(element_size(th[["text"]]), 9)
  expect_error(plot(dec, col_gray = NA),
               "'col_gray' must be a single color name", fixed = TRUE)
  expect_warning(plot(dec, overlay = TRUE),
                 "Ignoring unknown argument(s) of plot.qdm_decision(): overlay.",
                 fixed = TRUE)
  bad <- dec
  bad$design <- NULL
  expect_error(plot(bad), "'x' must be an unmodified result of qdm_decide().",
               fixed = TRUE)
})
