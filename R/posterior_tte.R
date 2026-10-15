#' Posterior distribution of the hazard of a time-to-event endpoint
#'
#' With the Gamma prior of shape \eqn{a} and rate \eqn{b}, external data of
#' \eqn{d_e} events and total observed time \eqn{E_e} with power prior weight
#' \eqn{a_0}, hypothetical data \eqn{(d_h, E_h)} and current data
#' \eqn{(d, E)}, the posterior distribution of the hazard is the Gamma
#' distribution with shape \eqn{a + a_0 d_e + d_h + d} and rate
#' \eqn{b + a_0 E_e + E_h + E}.
#'
#' @param prior An object of class \code{"qdm_prior_gamma"}.
#' @param events,exposure Numeric vectors of current numbers of events and
#'   total observed times (0 for none).
#' @param external External data (\code{"qdm_borrow"}) or \code{NULL}.
#' @param hypothetical Hypothetical control (\code{"qdm_hypothetical"}) or
#'   \code{NULL}.
#'
#' @return A list with numeric vectors \code{shape} and \code{rate}.
#' @noRd
posterior_tte <- function(prior, events, exposure, external = NULL,
                          hypothetical = NULL) {
  shape <- prior$shape + events
  rate <- prior$rate + exposure
  if (!is.null(external)) {
    shape <- shape + external$a0 * external$data$events
    rate <- rate + external$a0 * external$data$exposure
  }
  if (!is.null(hypothetical)) {
    shape <- shape + hypothetical$data$events
    rate <- rate + hypothetical$data$exposure
  }
  list(shape = shape, rate = rate)
}
