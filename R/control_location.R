#' Location of the distribution of the control group without a concurrent
#' control
#'
#' Used to shift the thresholds in plots of operating characteristics of a
#' design without a concurrent control group. For a binary endpoint it is the
#' mean of the Beta distribution of the control response rate; for a
#' continuous endpoint it is the location of the t distribution of the
#' control mean (the posterior mean given the external control data, or
#' \eqn{\mu_{0c}} of a hypothetical control).
#'
#' @param design An object of class \code{"qdm_design"} without a concurrent
#'   control group.
#'
#' @return A number.
#' @noRd
control_location <- function(design) {
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
    }
  )
}
