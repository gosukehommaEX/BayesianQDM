#' Predictive probabilities of a decision on one endpoint as bars
#'
#' @param x An object of class \code{"qdm_decision"} for one endpoint and a
#'   criterion based on posterior predictive probability.
#' @param cols Named colors of the decisions (\code{decision_colors()}).
#' @param style The style of the plot (\code{plot_style()}).
#'
#' @return A ggplot object with a bar for the predictive probability of the
#'   Go region and one for the NoGo region, and their thresholds.
#' @noRd
plot_decision_bars <- function(x, cols, style) {
  crit <- x$design$criterion
  labels <- event_labels(crit, 3L)
  bars <- data.frame(
    pos = 1:2,
    prob = unname(c(x$prob_go, x$prob_nogo)),
    gamma = c(crit$gamma_go, crit$gamma_nogo),
    event = factor(c("Go", "NoGo"), levels = c("Go", "NoGo"))
  )
  bars$text <- formatC(bars$prob, format = "f", digits = 3L)
  ggplot2::ggplot(bars) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$pos - 0.3,
                                    xmax = .data$pos + 0.3, ymin = 0,
                                    ymax = .data$prob, fill = .data$event),
                       alpha = 0.7) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$pos - 0.4,
                                       xend = .data$pos + 0.4,
                                       y = .data$gamma, yend = .data$gamma),
                          linewidth = style$linewidth, colour = "gray20") +
    ggplot2::geom_text(ggplot2::aes(x = .data$pos, y = .data$prob,
                                    label = .data$text),
                       vjust = -0.4, size = style$label_size) +
    ggplot2::scale_fill_manual(values = cols[c("Go", "NoGo")],
                               guide = "none") +
    ggplot2::scale_x_continuous(breaks = 1:2,
                                labels = c(labels$go, labels$nogo),
                                limits = c(0.4, 2.6)) +
    ggplot2::scale_y_continuous(limits = c(0, 1.08),
                                breaks = seq(0, 1, by = 0.2)) +
    ggplot2::labs(x = NULL, y = "Posterior predictive probability",
                  title = sprintf("Decision: %s", x$decision),
                  subtitle = decision_subtitle(x),
                  caption = "The horizontal segments mark the thresholds.") +
    style$theme +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank())
}
