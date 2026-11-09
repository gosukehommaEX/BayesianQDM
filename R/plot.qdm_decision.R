#' Plot a decision
#'
#' Shows how the observed data lead to the decision, for example to explain
#' a decision to a clinical team.
#'
#' For one endpoint and a criterion based on posterior probability, the plot
#' shows the posterior density of the treatment effect \eqn{\theta}, with
#' the target value and the minimum acceptable value (dotted lines) and the
#' Go and NoGo regions shaded; the subtitle gives the posterior
#' probabilities of the regions, their thresholds and whether the criteria
#' are met. The horizontal axis is the treatment effect for every control
#' group (for a fixed control, \eqn{\theta} is the parameter of the
#' treatment group minus the fixed value). The density is computed from the
#' exact distribution function of the effect (numerical integration for a
#' continuous endpoint, whatever the method of the design), on a grid that
#' contains the thresholds and covers the central 99.9 percent of the
#' distribution and the thresholds; the probabilities in the subtitle are
#' those of the decision, computed by the method of the design.
#'
#' For one endpoint and a criterion based on posterior predictive
#' probability, the plot shows the predictive probabilities of the Go and
#' NoGo regions as bars, with their thresholds (horizontal segments).
#'
#' For two endpoints, the plot shows the regions of the plane of the two
#' effects (see \code{\link{qdm_regions}}; for a predictive criterion, the
#' effects observed in the future study, marked "future" in the labels of the
#' axes), with the favorable side of each endpoint at the right and at the
#' top, and the posterior (or predictive) probability of each region; the
#' shade of each region grows with its probability. With
#' \code{rule = "joint"}, the regions of the Go and NoGo
#' criteria have the colors of these decisions.
#'
#' @param x An object of class \code{"qdm_decision"}.
#' @param col_go,col_nogo,col_gray Colors of the Go region, the NoGo region
#'   and the other regions.
#' @param linewidth Width of the density curve, in the units of the
#'   \code{linewidth} aesthetic of ggplot2.
#' @param label_size Size of the text drawn in the plot, in millimeters, or
#'   \code{NULL} for 80 percent of \code{base_size}.
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
#' dec <- qdm_decide(des, data_t = data_binary(12, 9),
#'                   data_c = data_binary(12, 4))
#' plot(dec)
#'
#' # Two continuous endpoints
#' des2 <- design_continuous2(
#'   n_t = 12, n_c = 12, prior_t = prior_vague(), prior_c = prior_vague(),
#'   criterion = criterion_posterior(theta_tv = c(10, 2), theta_mav = c(5, 1),
#'                                   gamma_go = 0.5, gamma_nogo = 0.3)
#' )
#' vc <- function(s1, s2, r) matrix(c(s1^2, r * s1 * s2, r * s1 * s2, s2^2), 2)
#' dec2 <- qdm_decide(
#'   des2,
#'   data_t = data_continuous2(12, c(18, 2.5), vc(60, 5, 0.4)),
#'   data_c = data_continuous2(12, c(-10, -1), vc(70, 6, 0.2))
#' )
#' plot(dec2)
#' @export
plot.qdm_decision <- function(x, col_go = "#658D1B", col_nogo = "#D91E49",
                              col_gray = "#939597", linewidth = 1,
                              label_size = NULL, base_size = 11,
                              legend_position = "bottom", ...) {
  if (!inherits(x$design, "qdm_design")) {
    stop("'x' must be an unmodified result of qdm_decide().", call. = FALSE)
  }
  warn_unused_dots(list(...), "plot.qdm_decision")
  cols <- decision_colors(col_go = col_go, col_nogo = col_nogo,
                          col_gray = col_gray)
  style <- plot_style(linewidth = linewidth, label_size = label_size,
                      base_size = base_size,
                      legend_position = legend_position)
  crit <- x$design$criterion
  if (crit$n_endpoints == 2L) {
    return(plot_decision_regions(x, cols, style))
  }
  if (crit$type == "predictive") {
    return(plot_decision_bars(x, cols, style))
  }
  plot_decision_density(x, cols, style)
}
