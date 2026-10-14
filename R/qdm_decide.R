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
#' @param design A design created by \code{\link{design_binary}},
#'   \code{\link{design_continuous}} or \code{\link{design_binary2}}.
#' @param data_t Observed data of the treatment group, created by
#'   \code{\link{data_binary}}, \code{\link{data_continuous}} or
#'   \code{\link{data_binary2}} (matching the endpoint of the design).
#' @param data_c Observed data of the concurrent control group, created by
#'   the same function as \code{data_t}; \code{NULL} (the default) for a
#'   design without a concurrent control group.
#' @param seed Random seed for probabilities computed by Monte Carlo
#'   simulation (a continuous design with \code{method = "MC"}, or posterior
#'   probabilities for two binary endpoints), or \code{NULL} to use the
#'   current random number stream. If given, the previous state of the random
#'   number generator is restored on exit. Not used otherwise.
#'
#' @return An object of class \code{"qdm_decision"}: a list with elements
#'   \describe{
#'     \item{decision}{\code{"Go"}, \code{"NoGo"}, \code{"Gray"} or
#'       \code{"Miss"}.}
#'     \item{prob_go, prob_nogo}{Probabilities of the Go and NoGo regions
#'       (for two endpoints with a marginal rule, one probability per
#'       endpoint).}
#'     \item{prob_regions}{For two endpoints, the probabilities of all
#'       regions (see \code{\link{qdm_regions}}); \code{NULL} otherwise.}
#'     \item{met_go, met_nogo}{Whether the Go and NoGo criteria are met.}
#'     \item{posterior_t, posterior_c}{Parameters of the posterior
#'       distributions of the treatment and control groups: for a binary
#'       endpoint, the Beta parameters \code{a} and \code{b}; for a
#'       continuous endpoint, the Normal-Inverse-chi-squared parameters
#'       \code{mu}, \code{kappa}, \code{nu} and \code{sigma} (on the standard
#'       deviation scale), so that the mean has the distribution
#'       \eqn{t_{\nu}(\mu, \sigma^2 / \kappa)}. For a hypothetical control,
#'       \code{posterior_c} gives the distribution of the control mean in the
#'       same form (see \code{\link{design_continuous}}). For two binary
#'       endpoints, the Dirichlet parameters of the cells \code{"00"},
#'       \code{"01"}, \code{"10"} and \code{"11"}.}
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
#'
#' # Continuous endpoint
#' des2 <- design_continuous(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_vague(), prior_c = prior_vague(),
#'   criterion = criterion_posterior(theta_tv = 30, theta_mav = 15,
#'                                   gamma_go = 0.6, gamma_nogo = 0.7)
#' )
#' qdm_decide(des2, data_t = data_continuous(12, 18, 25),
#'            data_c = data_continuous(12, -12, 25))
#' @export
qdm_decide <- function(design, data_t, data_c = NULL, seed = NULL) {
  check_design(design)
  if (!is.null(seed)) {
    check_numeric(seed, "seed", integer = TRUE)
  }
  res <- with_seed(seed, switch(
    design$endpoint,
    binary = decide_binary(design, data_t, data_c),
    continuous = decide_continuous(design, data_t, data_c),
    binary2 = decide_binary2(design, data_t, data_c)
  ))
  crit <- design$criterion
  met_go <- meets_criterion(res$prob_go, crit$gamma_go, crit$rule)
  met_nogo <- meets_criterion(res$prob_nogo, crit$gamma_nogo, crit$rule)
  structure(list(decision = decision_label(met_go, met_nogo),
                 prob_go = res$prob_go, prob_nogo = res$prob_nogo,
                 met_go = met_go, met_nogo = met_nogo,
                 prob_regions = res$prob_regions,
                 posterior_t = res$posterior_t,
                 posterior_c = res$posterior_c, design = design,
                 data_t = data_t, data_c = data_c),
            class = "qdm_decision")
}
