#' Gamma prior for a hazard rate
#'
#' Creates a Gamma prior distribution for the hazard rate \eqn{\lambda} of
#' one group in an analysis of a time-to-event endpoint with exponential
#' event times. The milestone survival probability at time \eqn{t_0} is
#' \eqn{S(t_0) = \exp(-\lambda t_0)}.
#'
#' The prior is equivalent to \code{shape} events observed during a total
#' observed time of \code{rate}, so small values give a weakly informative
#' prior: \code{prior_gamma(0.1, 1)} carries the information of a tenth of an
#' event.
#'
#' @param shape,rate Positive shape and rate parameters (the prior mean of
#'   the hazard is \code{shape / rate}).
#'
#' @return An object of class \code{c("qdm_prior_gamma", "qdm_prior")} with
#'   elements \code{shape} and \code{rate}.
#'
#' @family priors
#' @examples
#' prior_gamma(0.1, 1)
#' @export
prior_gamma <- function(shape, rate) {
  check_numeric(shape, "shape", lower = 0, lower_open = TRUE)
  check_numeric(rate, "rate", lower = 0, lower_open = TRUE)
  structure(list(shape = shape, rate = rate),
            class = c("qdm_prior_gamma", "qdm_prior"))
}
