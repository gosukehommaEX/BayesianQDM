#' Simulated studies of a design with two continuous endpoints and their
#' probabilities
#'
#' Simulates the sufficient statistics of \code{nsim} studies under each
#' scenario with the same random numbers for all scenarios (see
#' \code{draws_continuous2()} and \code{stats_continuous2()}); the random
#' numbers are drawn for the treatment group and then for the concurrent
#' control group. For each study, the probabilities of the Go and NoGo
#' regions are computed by the method of the design (with
#' \code{method = "MC"}, the draws are made after the random numbers of the
#' studies, separately for each study and scenario).
#'
#' @param design An object of class \code{"qdm_design_continuous2"}.
#' @param scenarios Data frame of scenarios (checked).
#' @param nsim Number of simulated studies.
#'
#' @return A list with elements \code{p_go} and \code{p_nogo} (matrices with
#'   one row per simulated study and one column per scenario, or lists of two
#'   such matrices for the marginal rules), \code{per_scenario}
#'   (\code{TRUE}), \code{weights} (matrix of the same size with all
#'   elements \code{1 / nsim}) and \code{method}.
#' @noRd
oc_sample_continuous2 <- function(design, scenarios, nsim) {
  nsim <- as.integer(round(nsim))
  crit <- design$criterion
  sets <- region_sets(crit)
  n_go <- length(sets$go)
  concurrent <- design$control == "concurrent"
  draws_t <- draws_continuous2(nsim, design$n_t)
  if (concurrent) {
    draws_c <- draws_continuous2(nsim, design$n_c)
  }
  n_scen <- nrow(scenarios)
  p <- matrix(0, nsim * n_scen, 2L * n_go)
  for (j in seq_len(n_scen)) {
    sc <- scenarios[j, , drop = FALSE]
    st <- stats_continuous2(draws_t, design$n_t, c(sc$mu_t1, sc$mu_t2),
                            c(sc$sd_t1, sc$sd_t2), sc$rho_t)
    post_t <- posterior_continuous2(design$prior_t, design$n_t, st$mean,
                                    st$cov, external = design$external_t)
    if (concurrent) {
      sc_c <- stats_continuous2(draws_c, design$n_c, c(sc$mu_c1, sc$mu_c2),
                                c(sc$sd_c1, sc$sd_c2), sc$rho_c)
      post_c <- control_posterior_continuous2(design, post_t, design$n_c,
                                              sc_c$mean, sc_c$cov)
    } else {
      post_c <- control_posterior_continuous2(design, post_t)
    }
    dist_t <- mean_dist_continuous2(post_t, crit$type, crit$m_t)
    dist_c <- if (is.null(post_c)) {
      NULL
    } else {
      mean_dist_continuous2(post_c, crit$type, crit$m_c)
    }
    p[(j - 1L) * nsim + seq_len(nsim), ] <- region_prob_continuous2(
      crit, dist_t, dist_c, c(sets$go, sets$nogo), design$method,
      design$nMC, fixed = design$fixed_c
    )
  }
  s <- split_sets(p, n_go, dim = c(nsim, n_scen))
  list(p_go = s$go, p_nogo = s$nogo, per_scenario = TRUE,
       weights = matrix(1 / nsim, nsim, n_scen),
       method = sprintf("simulation of %d studies; probabilities by %s",
                        nsim, method_label(design$method, design$nMC,
                                           design$control == "fixed")))
}
