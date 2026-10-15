#' Location of the distribution of the control group without a concurrent
#' control
#'
#' Used to shift the thresholds in plots of operating characteristics of a
#' design without a concurrent control group. For a fixed control it is the
#' fixed value. Otherwise, for a binary endpoint it is the mean of the Beta
#' distribution of the control response rate; for a continuous endpoint it
#' is the location of the t distribution of the control mean (the posterior
#' mean given the external control data, or \eqn{\mu_{0c}} of a hypothetical
#' control); for a time-to-event endpoint it is the mean of the milestone
#' survival probability, \eqn{(b / (b + t_0))^a} for a control hazard with
#' the Gamma distribution of shape \eqn{a} and rate \eqn{b}.
#'
#' @param design An object of class \code{"qdm_design"} for one endpoint
#'   without a concurrent control group.
#'
#' @return A number.
#' @noRd
control_location <- function(design) {
  if (design$control == "fixed") {
    return(design$fixed_c)
  }
  switch(
    design$endpoint,
    binary = {
      post <- posterior_binary(design$prior_c, 0, 0,
                               external = design$external_c,
                               hypothetical = design$hypothetical_c)
      post$a / (post$a + post$b)
    },
    continuous = if (design$control == "hypothetical") {
      design$hypothetical_c$mean
    } else {
      posterior_continuous(design$prior_c, 0L,
                           external = design$external_c)$mu
    },
    tte = {
      post <- posterior_tte(design$prior_c, 0, 0,
                            external = design$external_c,
                            hypothetical = design$hypothetical_c)
      (post$rate / (post$rate + design$milestone))^post$shape
    }
  )
}
