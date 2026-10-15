#' Horizontal axis of a plot of operating characteristics
#'
#' For a design with a concurrent control group, the axis is the true
#' treatment effect of each scenario and the thresholds of the criterion are
#' drawn at their values. Otherwise the axis is the true parameter of the
#' treatment group, and the thresholds are shifted by the fixed control value,
#' or by the mean (binary and time-to-event endpoints) or the location
#' (continuous endpoint) of the distribution of the control group. For a
#' time-to-event endpoint the parameter is the milestone survival
#' probability \eqn{\exp(-\log(2) t_0 / M)} for the median survival time
#' \eqn{M} of the scenario.
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
  par <- switch(
    design$endpoint,
    binary = list(t = x$pi_t, c = x$pi_c,
                  xlab = expression(theta == pi[t] - pi[c]),
                  xlab_t = expression(pi[t]), center = "mean"),
    continuous = list(t = x$mu_t, c = x$mu_c,
                      xlab = expression(theta == mu[t] - mu[c]),
                      xlab_t = expression(mu[t]), center = "location"),
    tte = list(t = exp(-log(2) * design$milestone / x$median_t),
               c = exp(-log(2) * design$milestone / x$median_c),
               xlab = expression(theta == S[t](t[0]) - S[c](t[0])),
               xlab_t = expression(S[t](t[0])), center = "mean")
  )
  if (design$control == "concurrent") {
    return(list(x = par$t - par$c, xlab = par$xlab, lines = lines,
                labels = labels, caption = NULL))
  }
  center <- control_location(design)
  what <- if (design$control == "fixed") {
    "fixed control value"
  } else {
    paste(par$center, "of the control distribution")
  }
  list(x = par$t, xlab = par$xlab_t, lines = center + lines, labels = labels,
       caption = sprintf("Thresholds shifted by the %s (%s)", what,
                         format_num(center, 3L)))
}
