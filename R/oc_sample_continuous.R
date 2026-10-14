#' Simulated outcomes of a continuous design and their probabilities
#'
#' Simulates the sufficient statistics of \code{nsim} proof-of-concept
#' studies under each scenario and computes the probabilities of the Go and
#' NoGo regions for each simulated study. For a group with \eqn{n} subjects,
#' true mean \eqn{\mu} and standard deviation \eqn{\sigma}, the sample mean
#' and standard deviation are
#' \eqn{\bar y = \mu + \sigma Z / \sqrt{n}} and
#' \eqn{s = \sigma \sqrt{X / (n - 1)}} with independent
#' \eqn{Z \sim N(0, 1)} and \eqn{X \sim \chi^2_{n - 1}}. The same draws of
#' \eqn{Z} and \eqn{X} are used for all scenarios (common random numbers), so
#' that differences between scenarios are not blurred by simulation noise.
#' The draws are made in the order \code{rnorm(nsim)},
#' \code{rchisq(nsim, n_t - 1)} for the treatment group and then the same for
#' the concurrent control group (if any).
#'
#' @param design An object of class \code{"qdm_design_continuous"}.
#' @param scenarios Data frame of scenarios (checked).
#' @param nsim Number of simulated studies.
#'
#' @return A list with elements \code{p_go} and \code{p_nogo} (matrices with
#'   one row per simulated study and one column per scenario),
#'   \code{per_scenario} (\code{TRUE}), \code{weights}
#'   (matrix of the same size with all elements \code{1 / nsim}) and
#'   \code{method}.
#' @noRd
oc_sample_continuous <- function(design, scenarios, nsim) {
  nsim <- as.integer(round(nsim))
  n_t <- design$n_t
  z_t <- stats::rnorm(nsim)
  x_t <- stats::rchisq(nsim, n_t - 1)
  concurrent <- design$control == "concurrent"
  if (concurrent) {
    n_c <- design$n_c
    z_c <- stats::rnorm(nsim)
    x_c <- stats::rchisq(nsim, n_c - 1)
  }
  n_scen <- nrow(scenarios)
  p_go <- matrix(0, nsim, n_scen)
  p_nogo <- matrix(0, nsim, n_scen)
  for (j in seq_len(n_scen)) {
    sd_t <- scenarios$sd_t[j]
    post_t <- posterior_continuous(
      design$prior_t, n_t, scenarios$mu_t[j] + sd_t * z_t / sqrt(n_t),
      sd_t * sqrt(x_t / (n_t - 1)), external = design$external_t
    )
    if (concurrent) {
      sd_c <- scenarios$sd_c[j]
      post_c <- control_posterior_continuous(
        design, post_t, n_c, scenarios$mu_c[j] + sd_c * z_c / sqrt(n_c),
        sd_c * sqrt(x_c / (n_c - 1))
      )
    } else {
      post_c <- control_posterior_continuous(design, post_t)
    }
    p <- region_prob_continuous(design$criterion, post_t, post_c,
                                design$method, design$nMC)
    p_go[, j] <- p$go
    p_nogo[, j] <- p$nogo
  }
  list(p_go = p_go, p_nogo = p_nogo, per_scenario = TRUE,
       weights = matrix(1 / nsim, nsim, n_scen),
       method = sprintf("simulation of %d studies; probabilities by %s",
                        nsim, method_label(design$method, design$nMC)))
}
