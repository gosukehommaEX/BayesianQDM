#' Line widths, sizes, text and legend of a plot
#'
#' Validates the style arguments of a \code{plot()} method and builds the
#' theme of the plot (\code{ggplot2::theme_bw()} with the given base size and
#' legend position).
#'
#' @param linewidth Width of the lines of the probabilities, in the units of
#'   the \code{linewidth} aesthetic of \pkg{ggplot2}.
#' @param point_size Size of the points, in millimeters.
#' @param label_size Size of the text drawn in the plot, in millimeters, or
#'   \code{NULL} for the size of the axis tick labels,
#'   \code{0.8 * base_size / .pt}.
#' @param base_size Base font size of the theme, in points.
#' @param legend_position One of \code{"bottom"}, \code{"top"},
#'   \code{"left"}, \code{"right"} and \code{"none"}.
#'
#' @return A list with elements \code{linewidth}, \code{point_size},
#'   \code{label_size} and \code{theme}.
#' @noRd
plot_style <- function(linewidth = 1, point_size = 1.5, label_size = NULL,
                       base_size = 11, legend_position = "bottom") {
  check_numeric(linewidth, "linewidth", lower = 0)
  check_numeric(point_size, "point_size", lower = 0)
  check_numeric(base_size, "base_size", lower = 0, lower_open = TRUE)
  if (is.null(label_size)) {
    label_size <- 0.8 * base_size / ggplot2::.pt
  } else {
    check_numeric(label_size, "label_size", lower = 0)
  }
  positions <- c("bottom", "top", "left", "right", "none")
  if (!is.character(legend_position) || length(legend_position) != 1L ||
      !(legend_position %in% positions)) {
    stop(paste0("'legend_position' must be one of \"bottom\", \"top\", ",
                "\"left\", \"right\" and \"none\"; use ggplot2::theme() ",
                "to place the legend inside the plot."), call. = FALSE)
  }
  theme <- ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.position = legend_position)
  list(linewidth = linewidth, point_size = point_size,
       label_size = label_size, theme = theme)
}
