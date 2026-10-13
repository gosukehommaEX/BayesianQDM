#' Operating characteristics of a design
#'
#' Computes the probabilities of Go, NoGo, Gray and Miss decisions of a design
#' under scenarios for the true parameters.
#'
#' For a binary endpoint, all possible outcomes of the proof-of-concept study
#' are enumerated and the probabilities are exact (Kang et al., 2026, Section
#' 4.3); \code{nsim} and \code{seed} are not used.
#'
#' The decision for each outcome is made as in \code{\link{qdm_decide}}. A
#' warning is given when Miss decisions (both criteria met) have positive
#' probability in some scenario.
#'
#' @param design A design created by \code{\link{design_binary}}.
#' @param scenarios A data frame with one row per scenario. For a binary
#'   endpoint the columns \code{pi_t} and \code{pi_c} give the true response
#'   rates of the treatment and control groups; \code{pi_c} is not needed for
#'   a design without a concurrent control group. Other columns, such as
#'   labels, are kept in the result.
#' @param nsim Number of simulated data sets for endpoints whose operating
#'   characteristics are computed by simulation.
#' @param seed Random seed for the simulation, or \code{NULL} to use the
#'   current random number stream. If given, the previous state of the random
#'   number generator is restored on exit.
#'
#' @return An object of class \code{c("qdm_oc", "data.frame")}: the scenarios
#'   with the columns \code{Go}, \code{NoGo}, \code{Gray} and \code{Miss}
#'   added. The attributes \code{design} and \code{method} record the design
#'   and the method of computation.
#'
#' @references
#' Kang, R., Yamaguchi, Y. and Han, C. (2026). Bayesian quantitative
#' decision-making framework in cell and gene therapy development for rare
#' diseases. \emph{Journal of Biopharmaceutical Statistics}.
#' \doi{10.1080/10543406.2026.2655410}
#'
#' @family operations
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' oc <- qdm_oc(des, data.frame(pi_t = seq(0.2, 0.8, by = 0.1), pi_c = 0.4))
#' oc
#' plot(oc)
#' @export
qdm_oc <- function(design, scenarios, nsim = 10000L, seed = NULL) {
  check_design(design)
  scenarios <- check_scenarios(design, scenarios)
  check_numeric(nsim, "nsim", lower = 1, integer = TRUE)
  if (!is.null(seed)) {
    check_numeric(seed, "seed", integer = TRUE)
  }
  sample <- with_seed(seed, oc_sample(design, scenarios, nsim))
  probs <- oc_from_sample(sample, design$criterion$gamma_go,
                          design$criterion$gamma_nogo)
  out <- cbind(scenarios, as.data.frame(probs))
  rownames(out) <- NULL
  warn_miss(out$Miss)
  structure(out, class = c("qdm_oc", "data.frame"), design = design,
            method = sample$method)
}
