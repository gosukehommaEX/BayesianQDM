#' Set the decision criterion of a design specification
#'
#' Sets the argument \code{criterion} of a design specification started by
#' \code{\link{qdm_endpoint}}. It is checked when the design is built, by the
#' \code{design_*()} function of the endpoint.
#'
#' @param x A design specification started by \code{\link{qdm_endpoint}},
#'   or a design, which is converted to the specification that builds it.
#' @param criterion Decision criterion created by
#'   \code{\link{criterion_posterior}} or \code{\link{criterion_predictive}},
#'   for one endpoint or for two endpoints as the endpoint requires.
#'
#' @return A design specification: \code{x}, or the specification of the
#'   design \code{x}, with the decision criterion set.
#'
#' @family design steps
#' @examples
#' qdm_endpoint("binary") |>
#'   qdm_criterion(criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                     gamma_go = 0.5, gamma_nogo = 0.8))
#' @export
qdm_criterion <- function(x, criterion) {
  spec_set(x, list(criterion = criterion), "qdm_criterion")
}
