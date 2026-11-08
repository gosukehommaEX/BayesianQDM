#' Set a hypothetical or fixed control of a design specification
#'
#' Sets the arguments \code{hypothetical_c} and \code{fixed_c} of a design
#' specification started by \code{\link{qdm_endpoint}}, for a study without
#' a concurrent control group (\code{n_c = 0}, see
#' \code{\link{qdm_sample_size}}). A concurrent control group and external
#' control data only (see \code{\link{qdm_external}}) need no call of this
#' function. The arguments are checked when the design is built, by the
#' \code{design_*()} function of the endpoint. All the arguments of the step
#' are set together: an argument that is not given returns to the default, so
#' a hypothetical control replaces a fixed control and the other way round.
#'
#' @param x A design specification started by \code{\link{qdm_endpoint}},
#'   or a design, which is converted to the specification that builds it.
#' @param hypothetical_c Hypothetical control, created by
#'   \code{\link{hypothetical}}, or \code{NULL}.
#' @param fixed_c Fixed value of the control parameter, or \code{NULL}. Its
#'   meaning depends on the endpoint (see \code{\link{design_binary}},
#'   \code{\link{design_continuous}}, \code{\link{design_binary2}},
#'   \code{\link{design_continuous2}} and \code{\link{design_tte}}).
#'
#' @return A design specification: \code{x}, or the specification of the
#'   design \code{x}, with the control set.
#'
#' @family design steps
#' @examples
#' # Single-arm study compared with a fixed control response rate of 0.2
#' qdm_endpoint("binary") |>
#'   qdm_sample_size(n_t = 20, n_c = 0) |>
#'   qdm_prior(prior_t = prior_beta(1, 1)) |>
#'   qdm_criterion(criterion_posterior(theta_tv = 0.2, theta_mav = 0.05,
#'                                     gamma_go = 0.8, gamma_nogo = 0.8)) |>
#'   qdm_control(fixed_c = 0.2)
#' @export
qdm_control <- function(x, hypothetical_c = NULL, fixed_c = NULL) {
  spec_set(x, list(hypothetical_c = hypothetical_c, fixed_c = fixed_c),
           "qdm_control")
}
