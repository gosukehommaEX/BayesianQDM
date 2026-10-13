#' Horizontal axis of a plot of operating characteristics
#'
#' For a design with a concurrent control group, the axis is the true
#' treatment effect of each scenario and the thresholds of the criterion are
#' drawn at their values. Otherwise the axis is the true parameter of the
#' treatment group, and the thresholds are shifted by the mean (binary
#' endpoint) or the location (continuous endpoint) of the distribution of the
#' control group.
#'
#' @param design An object of class \code{"qdm_design"}.
#' @param x The scenarios (a data frame).
#'
#' @return A list with elements \code{x} (axis values), \code{xlab}
#'   (expression), \code{lines} (positions of the thresholds), \code{labels}
#'   (plotmath labels of the thresholds) and \code{caption}.
#' @noRd
oc_axis <- function(design, x) {
  crit <- design$criterion
  if (crit$type == "posterior") {
    lines <- c(crit$theta_mav, crit$theta_tv)
    labels <- c("theta[MAV]", "theta[TV]")
  } else {
    lines <- crit$theta_null
    labels <- "theta[NULL]"
  }
  binary <- design$endpoint == "binary"
  if (design$control == "concurrent") {
    if (binary) {
      return(list(x = x$pi_t - x$pi_c,
                  xlab = expression(theta == pi[t] - pi[c]),
                  lines = lines, labels = labels, caption = NULL))
    }
    return(list(x = x$mu_t - x$mu_c,
                xlab = expression(theta == mu[t] - mu[c]),
                lines = lines, labels = labels, caption = NULL))
  }
  center <- control_location(design)
  list(x = if (binary) x$pi_t else x$mu_t,
       xlab = if (binary) expression(pi[t]) else expression(mu[t]),
       lines = center + lines, labels = labels,
       caption = sprintf(paste("Thresholds shifted by the %s of the",
                               "control distribution (%s)"),
                         if (binary) "mean" else "location",
                         format_num(center, 3L)))
}
