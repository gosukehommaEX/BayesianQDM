plot_labels <- function(p) {
  if (exists("get_labs", envir = asNamespace("ggplot2"), inherits = FALSE)) {
    ggplot2::get_labs(p)
  } else {
    p$labels
  }
}

test_that("plot.qdm_boundary() draws the decision of each grid point", {
  b <- qdm_boundary(ref_design("D8"), n_t = 10:20)
  p <- plot(b)
  expect_s3_class(p, "ggplot")
  tiles <- ggplot2::ggplot_build(p)$data[[1]]
  expect_identical(nrow(tiles), nrow(b))
  idx <- match(paste(tiles$x, tiles$y), paste(b$n_t, b$y_t))
  expect_false(anyNA(idx))
  cols <- c(Go = "#658D1B", NoGo = "#D91E49", Gray = "#939597",
            Miss = "#E69F00")
  expect_identical(unname(tiles$fill), unname(cols[as.character(
    b$decision[idx])]))
  labs <- plot_labels(p)
  expect_identical(labs$x, "Treatment group size")
  expect_identical(labs$y, "Responders in the treatment group")
  expect_identical(labs$fill, "Decision")
})

test_that("plot.qdm_boundary() shows the control data in a caption", {
  b <- qdm_boundary(ref_design_tte("T2"), events_t = 0:20,
                    exposure_t = seq(20, 200, by = 20),
                    data_c = data_tte(9, 35))
  p <- plot(b)
  labs <- plot_labels(p)
  expect_match(labs$caption, "^Control data: ")
  expect_identical(labs$x, "Total observed time of the treatment group")
  tiles <- ggplot2::ggplot_build(p)$data[[1]]
  expect_setequal(unique(tiles$fill),
                  c(Go = "#658D1B", NoGo = "#D91E49", Gray = "#939597",
                    Miss = "#E69F00")[unique(as.character(b$decision))])
})

test_that("plot.qdm_boundary() needs the axes attribute", {
  b <- qdm_boundary(ref_design("D8"))
  attr(b, "axes") <- NULL
  expect_error(plot(b),
               "'x' must be an unmodified result of qdm_boundary\\(\\).")
})

test_that("plot.qdm_boundary() takes colors, text size and legend", {
  b <- qdm_boundary(ref_design("D8"), n_t = 10:20)
  p <- plot(b, col_go = "black", col_nogo = "blue", col_gray = "white",
            base_size = 8, legend_position = "none")
  tiles <- ggplot2::ggplot_build(p)$data[[1]]
  idx <- match(paste(tiles$x, tiles$y), paste(b$n_t, b$y_t))
  cols <- c(Go = "black", NoGo = "blue", Gray = "white", Miss = "#E69F00")
  expect_identical(unname(tiles$fill), unname(cols[as.character(
    b$decision[idx])]))
  expect_true(all(c("Go", "NoGo", "Gray") %in% as.character(b$decision)))
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
