#' Plot operating characteristics
#'
#' Draws the probabilities of Go, NoGo and Gray decisions (and Miss decisions
#' if they have positive probability) against the true treatment effect of
#' the scenarios. For a design without a concurrent control group, the
#' horizontal axis is the true response rate (binary endpoint) or mean
#' (continuous endpoint) of the treatment group, and the thresholds of the
#' criterion are shifted by the mean (binary) or location (continuous) of the
#' distribution of the control group.
#'
#' For two endpoints, the plot has one heat map per decision, with the true
#' treatment effects of the two endpoints (or the true parameters of the
#' treatment group without a concurrent control group) on the axes and the
#' probability of the decision in each cell; scenarios that are not evaluable
#' are marked "NE". The scenarios must have distinct values on the axes, so
#' scenarios with several correlations (or several parameters of the control
#' group) are plotted separately.
#'
#' @param x An object of class \code{"qdm_oc"}.
#' @param ... Not used.
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
#' oc <- qdm_oc(des, data.frame(pi_t = seq(0.2, 0.9, by = 0.05), pi_c = 0.4))
#' plot(oc)
#' @export
plot.qdm_oc <- function(x, ...) {
  design <- attr(x, "design")
  if (is.null(design)) {
    stop("'x' must be an unmodified result of qdm_oc().", call. = FALSE)
  }
  if (design$criterion$n_endpoints == 2L) {
    return(plot_oc_two(x, design))
  }
  axis <- oc_axis(design, x)
  decisions <- c("Go", "NoGo", "Gray", "Miss")
  if (all(x$Miss <= 0)) {
    decisions <- decisions[1:3]
  }
  ord <- order(axis$x)
  long <- data.frame(
    x = rep(axis$x[ord], length(decisions)),
    prob = unlist(lapply(decisions, function(d) x[[d]][ord]),
                  use.names = FALSE),
    decision = factor(rep(decisions, each = nrow(x)), levels = decisions)
  )
  cols <- c(Go = "#658D1B", NoGo = "#D91E49", Gray = "#939597",
            Miss = "#E69F00")[decisions]
  ltys <- c(Go = "solid", NoGo = "twodash", Gray = "dashed",
            Miss = "dotted")[decisions]
  ggplot2::ggplot(long, ggplot2::aes(x = .data$x, y = .data$prob,
                                     colour = .data$decision,
                                     linetype = .data$decision)) +
    ggplot2::geom_vline(xintercept = axis$lines, colour = "gray40",
                        linetype = "dotted") +
    ggplot2::annotate("text", x = axis$lines, y = 1, label = axis$labels,
                      parse = TRUE, hjust = -0.15, vjust = 1,
                      colour = "gray30") +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_linetype_manual(values = ltys) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                breaks = seq(0, 1, by = 0.2)) +
    ggplot2::labs(x = axis$xlab, y = "Probability", colour = "Decision",
                  linetype = "Decision", caption = axis$caption) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank())
}
