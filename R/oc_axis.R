#' Horizontal axis of a plot of operating characteristics
#'
#' For a design with a concurrent control group, the axis is the true
#' treatment effect of each scenario and the thresholds of the criterion are
#' drawn at their values. Otherwise the axis is the true parameter of the
#' treatment group, and the thresholds are shifted by the mean of the
#' distribution of the control group.
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
  if (design$control == "concurrent") {
    return(list(x = x$pi_t - x$pi_c,
                xlab = expression(theta == pi[t] - pi[c]),
                lines = lines, labels = labels, caption = NULL))
  }
  post_c <- posterior_binary(design$prior_c, 0, 0,
                             external = design$external_c,
                             hypothetical = design$hypothetical_c)
  center <- post_c$a / (post_c$a + post_c$b)
  list(x = x$pi_t, xlab = expression(pi[t]), lines = center + lines,
       labels = labels,
       caption = sprintf(paste("Thresholds shifted by the mean of the",
                               "control distribution (%s)"),
                         format_num(center, 3L)))
}
