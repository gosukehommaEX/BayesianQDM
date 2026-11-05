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
#' the Gamma distribution of shape \eqn{a} and rate \eqn{b}, or the mean of
#' the hazard, \eqn{a / b}, for the hazard ratio (\eqn{\log(2)} divided by
#' the median survival time for a fixed control).
#'
#' For two endpoints, it is one value per endpoint: the fixed values, the
#' means of the response rates of the Dirichlet distribution of the control
#' cells (two binary endpoints), or the location of the distribution of the
#' control means (two continuous endpoints, as for one continuous endpoint).
#'
#' @param design An object of class \code{"qdm_design"} without a concurrent
#'   control group.
#'
#' @return A number, or two numbers for two endpoints.
#' @noRd
control_location <- function(design) {
  hr <- identical(design$effect, "hazard_ratio")
  if (design$control == "fixed") {
    return(if (hr) log(2) / design$fixed_c else design$fixed_c)
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
    binary2 = {
      alpha <- posterior_binary2(design$prior_c, c(0, 0, 0, 0),
                                 external = design$external_c,
                                 hypothetical = design$hypothetical_c)[1L, ]
      # The cells are "00", "01", "10" and "11" for endpoints 1 and 2
      unname(c(alpha["10"] + alpha["11"], alpha["01"] + alpha["11"]) /
               sum(alpha))
    },
    continuous2 = if (design$control == "hypothetical") {
      design$hypothetical_c$mean
    } else {
      as.numeric(posterior_continuous2(design$prior_c, 0L,
                                       external = design$external_c)$mu)
    },
    tte = {
      post <- posterior_tte(design$prior_c, 0, 0,
                            external = design$external_c,
                            hypothetical = design$hypothetical_c)
      if (hr) {
        post$shape / post$rate
      } else {
        (post$rate / (post$rate + design$milestone))^post$shape
      }
    }
  )
}
