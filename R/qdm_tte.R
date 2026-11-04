#' Set the arguments of a time-to-event design specification
#'
#' Sets the arguments of \code{\link{design_tte}} that are specific to a
#' time-to-event endpoint: the effect measure, the milestone time and the
#' conduct of the study and of the future study. The specification must
#' have been started by \code{qdm_endpoint("tte")}. The arguments are checked
#' when the design is built, by \code{\link{design_tte}}.
#'
#' @param x A design specification started by \code{qdm_endpoint("tte")}.
#' @param milestone Milestone time \eqn{t_0} for \code{effect = "milestone"},
#'   or \code{NULL}.
#' @param effect Effect measure, \code{"milestone"} or
#'   \code{"hazard_ratio"}, or \code{NULL} for the default of
#'   \code{\link{design_tte}} (\code{"milestone"}).
#' @param conduct Conduct of the study, created by
#'   \code{\link{conduct_tte}}, or \code{NULL}.
#' @param conduct_future Conduct of the future study for a predictive
#'   criterion, created by \code{\link{conduct_tte}}, or \code{NULL} for the
#'   default of \code{\link{design_tte}} (\code{conduct}).
#'
#' @return \code{x} with the arguments set.
#'
#' @family design steps
#' @examples
#' qdm_endpoint("tte") |>
#'   qdm_sample_size(n_t = 20, n_c = 20) |>
#'   qdm_prior(prior_t = prior_gamma(0.1, 1),
#'             prior_c = prior_gamma(0.1, 1)) |>
#'   qdm_criterion(criterion_posterior(theta_tv = 0.15, theta_mav = 0,
#'                                     gamma_go = 0.6, gamma_nogo = 0.8)) |>
#'   qdm_tte(milestone = 3,
#'           conduct = conduct_tte(accrual_time = c(0, 5, 10, 15),
#'                                 accrual_rate = c(0.4, 2, 4.5),
#'                                 dropout_hazard = 0.00226,
#'                                 cut_time = 15.44)) |>
#'   qdm_build()
#' @export
qdm_tte <- function(x, milestone = NULL, effect = NULL, conduct = NULL,
                    conduct_future = NULL) {
  check_spec(x)
  if (x$endpoint != "tte") {
    stop("qdm_tte() is used only with qdm_endpoint(\"tte\").",
         call. = FALSE)
  }
  spec_set(x, list(milestone = milestone, effect = effect, conduct = conduct,
                   conduct_future = conduct_future), "qdm_tte")
}
