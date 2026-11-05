#' Build a design from a design specification
#'
#' Creates the design specified step by step from \code{\link{qdm_endpoint}},
#' by calling the \code{design_*()} function of the endpoint
#' (\code{\link{design_binary}}, \code{\link{design_continuous}},
#' \code{\link{design_binary2}}, \code{\link{design_continuous2}} or
#' \code{\link{design_tte}}) with the arguments set. The design is identical
#' to the one created by calling that function directly with the same
#' arguments, and invalid arguments give the same errors. The sample sizes
#' (\code{\link{qdm_sample_size}}), the prior distribution of the treatment
#' group (\code{\link{qdm_prior}}) and the decision criterion
#' (\code{\link{qdm_criterion}}) must have been set.
#'
#' @param x A design specification started by \code{\link{qdm_endpoint}},
#'   or a design, for which the same design is returned.
#'
#' @return A design, an object of class \code{"qdm_design"}, as returned by
#'   the \code{design_*()} function of the endpoint.
#'
#' @family design steps
#' @examples
#' qdm_endpoint("binary") |>
#'   qdm_sample_size(n_t = 12, n_c = 12) |>
#'   qdm_prior(prior_t = prior_beta(0.5, 0.5),
#'             prior_c = prior_beta(0.5, 0.5)) |>
#'   qdm_criterion(criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                     gamma_go = 0.5, gamma_nogo = 0.8)) |>
#'   qdm_build()
#' @export
qdm_build <- function(x) {
  x <- check_spec(x)
  unset <- spec_missing(x)
  if (length(unset) > 0L) {
    stop(sprintf("The design specification is incomplete: set %s.",
                 paste(unset, collapse = ", ")), call. = FALSE)
  }
  do.call(design_function(x$endpoint), x$args)
}
