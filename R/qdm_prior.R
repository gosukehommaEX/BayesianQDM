#' Set the prior distributions of a design specification
#'
#' Sets the arguments \code{prior_t} and \code{prior_c} of a design
#' specification started by \code{\link{qdm_endpoint}}. They are checked
#' when the design is built, by the \code{design_*()} function of the
#' endpoint.
#'
#' @param x A design specification started by \code{\link{qdm_endpoint}}.
#' @param prior_t,prior_c Prior distributions of the treatment and control
#'   groups, created by the \code{prior_*()} function of the endpoint (see
#'   \code{\link{design_binary}}, \code{\link{design_continuous}},
#'   \code{\link{design_binary2}}, \code{\link{design_continuous2}} and
#'   \code{\link{design_tte}}). \code{prior_c} is \code{NULL} (the default)
#'   for a fixed control.
#'
#' @return \code{x} with the prior distributions set.
#'
#' @family design steps
#' @examples
#' qdm_endpoint("binary") |>
#'   qdm_prior(prior_t = prior_beta(0.5, 0.5),
#'             prior_c = prior_beta(0.5, 0.5))
#' @export
qdm_prior <- function(x, prior_t, prior_c = NULL) {
  spec_set(x, list(prior_t = prior_t, prior_c = prior_c), "qdm_prior")
}
