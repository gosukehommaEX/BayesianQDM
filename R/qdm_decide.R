#' Go/NoGo decision for observed data
#'
#' Computes the probabilities of the Go and NoGo regions given the observed
#' proof-of-concept data, and makes the decision of a design.
#'
#' The decision is Go if only the Go criterion is met, NoGo if only the NoGo
#' criterion is met, Gray if neither is met, and Miss if both are met (see
#' \code{\link{criterion_posterior}} and \code{\link{criterion_predictive}}).
#' A probability that is below a threshold by less than \eqn{10^{-10}} is
#' treated as equal to the threshold, so that ties that arise in symmetric
#' settings are not decided by rounding errors.
#'
#' The sample sizes of the observed data may differ from those of the design,
#' which are used only by \code{\link{qdm_oc}} and \code{\link{qdm_calibrate}}.
#'
#' @param design A design created by \code{\link{design_binary}}.
#' @param data_t Observed data of the treatment group, created by
#'   \code{\link{data_binary}}.
#' @param data_c Observed data of the concurrent control group, created by
#'   \code{\link{data_binary}}; \code{NULL} (the default) for a design without
#'   a concurrent control group.
#'
#' @return An object of class \code{"qdm_decision"}: a list with elements
#'   \describe{
#'     \item{decision}{\code{"Go"}, \code{"NoGo"}, \code{"Gray"} or
#'       \code{"Miss"}.}
#'     \item{prob_go, prob_nogo}{Probabilities of the Go and NoGo regions.}
#'     \item{met_go, met_nogo}{Whether the Go and NoGo criteria are met.}
#'     \item{posterior_t, posterior_c}{Parameters of the posterior
#'       distributions of the treatment and control groups (for a binary
#'       endpoint, the Beta parameters \code{a} and \code{b}).}
#'     \item{design, data_t, data_c}{The arguments.}
#'   }
#'
#' @family operations
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' qdm_decide(des, data_t = data_binary(12, 9), data_c = data_binary(12, 4))
#' @export
qdm_decide <- function(design, data_t, data_c = NULL) {
  check_design(design)
  res <- switch(design$endpoint,
                binary = decide_binary(design, data_t, data_c))
  crit <- design$criterion
  met_go <- meets_threshold(res$prob_go, crit$gamma_go)
  met_nogo <- meets_threshold(res$prob_nogo, crit$gamma_nogo)
  structure(list(decision = decision_label(met_go, met_nogo),
                 prob_go = res$prob_go, prob_nogo = res$prob_nogo,
                 met_go = met_go, met_nogo = met_nogo,
                 posterior_t = res$posterior_t,
                 posterior_c = res$posterior_c, design = design,
                 data_t = data_t, data_c = data_c),
            class = "qdm_decision")
}
