#' Plot a calibration
#'
#' Draws the probability that the Go criterion is met under
#' \code{scenario_go} and the probability that the NoGo criterion is met under
#' \code{scenario_nogo} against the threshold, with the targets (horizontal
#' lines) and the chosen thresholds (vertical lines).
#'
#' The plot is a \code{ggplot} object, so it can be changed further by adding
#' \pkg{ggplot2} functions, for example \code{ggplot2::labs()} and
#' \code{ggplot2::theme()}.
#'
#' @inheritParams plot.qdm_oc
#' @param x An object of class \code{"qdm_calibration"}.
#' @param col_go,col_nogo Colors of the curves, targets and thresholds of the
#'   Go and NoGo criteria: color names or hexadecimal codes. The defaults
#'   are green and red.
#' @param linewidth Width of the calibration curves, in the units of the
#'   \code{linewidth} aesthetic of \pkg{ggplot2}.
#' @param base_size Base font size of \code{ggplot2::theme_bw()}, in points,
#'   which sets the size of the axis labels, tick labels and legend.
#' @param ... Not used; arguments given here are ignored with a warning.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' cal <- qdm_calibrate(des, data.frame(pi_t = 0.4, pi_c = 0.4),
#'                      data.frame(pi_t = 0.7, pi_c = 0.4),
#'                      target_go = 0.05, target_nogo = 0.2)
#' plot(cal)
#' plot(cal, linewidth = 0.6, base_size = 14, legend_position = "top")
#' @export
plot.qdm_calibration <- function(x, col_go = "#658D1B", col_nogo = "#D91E49",
                                 linewidth = 1, base_size = 11,
                                 legend_position = "bottom", ...) {
  warn_unused_dots(list(...), "plot.qdm_calibration")
  cols <- decision_colors(col_go = col_go, col_nogo = col_nogo)
  style <- plot_style(linewidth = linewidth, base_size = base_size,
                      legend_position = legend_position)
  labels <- c("Go criterion met under scenario_go",
              "NoGo criterion met under scenario_nogo")
  long <- data.frame(
    gamma = rep(x$curve$gamma, 2L),
    prob = c(x$curve$prob_go, x$curve$prob_nogo),
    event = factor(rep(labels, each = nrow(x$curve)), levels = labels)
  )
  cols <- stats::setNames(unname(cols), labels)
  targets <- data.frame(y = c(x$target_go, x$target_nogo),
                        event = factor(labels, levels = labels))
  chosen <- data.frame(x = c(x$gamma_go, x$gamma_nogo),
                       event = factor(labels, levels = labels))
  chosen <- chosen[!is.na(chosen$x), , drop = FALSE]
  ggplot2::ggplot(long, ggplot2::aes(x = .data$gamma, y = .data$prob,
                                     colour = .data$event)) +
    ggplot2::geom_hline(data = targets,
                        ggplot2::aes(yintercept = .data$y,
                                     colour = .data$event),
                        linetype = "dashed") +
    ggplot2::geom_vline(data = chosen,
                        ggplot2::aes(xintercept = .data$x,
                                     colour = .data$event),
                        linetype = "dotted") +
    ggplot2::geom_step(linewidth = style$linewidth) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                breaks = seq(0, 1, by = 0.2)) +
    ggplot2::labs(x = expression(gamma), y = "Probability", colour = NULL) +
    style$theme +
    ggplot2::theme(legend.direction = "vertical",
                   panel.grid.minor = ggplot2::element_blank())
}
