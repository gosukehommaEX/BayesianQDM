#' Plot a search of the sample size
#'
#' Draws the probabilities that the Go criterion is met under
#' \code{scenario_go} and that the NoGo criterion is met under
#' \code{scenario_nogo} against the number of subjects in the treatment
#' group, with the targets (dashed lines, one gray line when the two targets
#' are equal) and the smallest candidate that meets them (dotted line). The
#' vertical axis runs from 0 to the largest probability or target, so that
#' small probabilities can be read.
#'
#' @param x An object of class \code{"qdm_n_search"}.
#' @param col_go,col_nogo Colors of the probabilities and targets of the Go
#'   and NoGo criteria.
#' @param linewidth Width of the lines, in the units of the \code{linewidth}
#'   aesthetic of ggplot2.
#' @param point_size Size of the points, in millimeters.
#' @param base_size Base font size of \code{ggplot2::theme_bw()}, in points.
#' @param legend_position Position of the legend: \code{"bottom"} (the
#'   default), \code{"top"}, \code{"left"}, \code{"right"} or \code{"none"}.
#' @param ... Not used; arguments given here are ignored with a warning.
#'
#' @return A ggplot object.
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' res <- qdm_find_n(des, data.frame(pi_t = 0.4, pi_c = 0.4),
#'                   data.frame(pi_t = 0.7, pi_c = 0.4), target_go = 0.1,
#'                   target_nogo = 0.1, n_t = seq(6, 30, by = 2))
#' plot(res)
#' @export
plot.qdm_n_search <- function(x, col_go = "#658D1B", col_nogo = "#D91E49",
                              linewidth = 1, point_size = 1.5,
                              base_size = 11, legend_position = "bottom",
                              ...) {
  warn_unused_dots(list(...), "plot.qdm_n_search")
  cols <- decision_colors(col_go = col_go, col_nogo = col_nogo)
  style <- plot_style(linewidth = linewidth, point_size = point_size,
                      base_size = base_size,
                      legend_position = legend_position)
  labels <- c("Go criterion met under scenario_go",
              "NoGo criterion met under scenario_nogo")
  curve <- x$curve
  long <- data.frame(
    n = rep(curve$n_t, 2L),
    prob = c(curve$prob_go, curve$prob_nogo),
    event = factor(rep(labels, each = nrow(curve)), levels = labels)
  )
  cols <- stats::setNames(unname(cols), labels)
  targets <- data.frame(y = c(x$target_go, x$target_nogo),
                        event = factor(labels, levels = labels))
  chosen <- data.frame(x = x$n_t)
  chosen <- chosen[!is.na(chosen$x), , drop = FALSE]
  # Equal targets are one line, which the colors of the two would hide
  target_line <- if (identical(x$target_go, x$target_nogo)) {
    ggplot2::geom_hline(yintercept = x$target_go, colour = "gray40",
                        linetype = "dashed")
  } else {
    ggplot2::geom_hline(data = targets,
                        ggplot2::aes(yintercept = .data$y,
                                     colour = .data$event),
                        linetype = "dashed")
  }
  ggplot2::ggplot(long, ggplot2::aes(x = .data$n, y = .data$prob,
                                     colour = .data$event)) +
    target_line +
    ggplot2::geom_vline(data = chosen,
                        ggplot2::aes(xintercept = .data$x),
                        colour = "gray40", linetype = "dotted") +
    ggplot2::geom_line(linewidth = style$linewidth) +
    ggplot2::geom_point(size = style$point_size) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::labs(x = "Number of subjects in the treatment group",
                  y = "Probability", colour = NULL) +
    style$theme +
    ggplot2::theme(legend.direction = "vertical",
                   panel.grid.minor = ggplot2::element_blank())
}
