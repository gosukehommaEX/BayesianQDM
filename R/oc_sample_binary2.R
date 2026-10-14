#' Outcomes of a design with two binary endpoints and their probabilities
#'
#' With \code{design$enumerate = TRUE}, all outcomes of the proof-of-concept
#' study are enumerated (Yamaguchi et al., submitted, equation A.5): the
#' probabilities of the Go and NoGo regions are computed for every pair of
#' treatment and control outcomes, and the probability of a pair under a
#' scenario is the product of the multinomial probabilities of the two
#' outcomes (kept as one matrix per group); enumeration is limited to 1e7
#' pairs of outcomes. Otherwise \code{nsim} studies are
#' simulated under each scenario with common random numbers (see
#' \code{simulate_counts_binary2()}); the uniform random numbers are drawn
#' for the treatment group and then for the concurrent control group, and the
#' probabilities of the regions are computed once for each distinct pair of
#' outcomes. A scenario whose correlation is infeasible for its response
#' rates is not evaluable: its probabilities are \code{NA} and a warning is
#' given.
#'
#' @param design An object of class \code{"qdm_design_binary2"}.
#' @param scenarios Data frame of scenarios (checked).
#' @param nsim Number of simulated studies (simulation only).
#'
#' @return A list with elements \code{p_go} and \code{p_nogo} (a matrix, or
#'   a list of two matrices for the marginal rules), \code{per_scenario}
#'   (whether the matrices have one column per scenario), either
#'   \code{weights} (simulation) or \code{weights_t} and \code{weights_c}
#'   (enumeration; one column per scenario), and \code{method}.
#' @noRd
oc_sample_binary2 <- function(design, scenarios, nsim) {
  crit <- design$criterion
  sets <- region_sets(crit)
  n_go <- length(sets$go)
  concurrent <- design$control == "concurrent"
  n_scen <- nrow(scenarios)
  cells_t <- cellprob_scenarios(scenarios$pi_t1, scenarios$pi_t2,
                                scenarios$rho_t)
  bad <- is.na(cells_t[, 1L])
  if (concurrent) {
    cells_c <- cellprob_scenarios(scenarios$pi_c1, scenarios$pi_c2,
                                  scenarios$rho_c)
    bad <- bad | is.na(cells_c[, 1L])
  }
  if (any(bad)) {
    warning(sprintf(paste("Scenario(s) %s are not evaluable: the",
                          "correlation is outside the feasible range for",
                          "the response rates. Their probabilities are NA."),
                    paste(which(bad), collapse = ", ")), call. = FALSE)
  }
  prob_label <- if (crit$type == "posterior") {
    sprintf("posterior probabilities by Monte Carlo simulation (nMC = %d)",
            design$nMC)
  } else {
    "predictive probabilities computed exactly"
  }
  fixed_c <- posterior_binary2(design$prior_c, c(0, 0, 0, 0),
                               external = design$external_c,
                               hypothetical = design$hypothetical_c)
  if (design$enumerate) {
    n_pairs <- choose(design$n_t + 3, 3) *
      (if (concurrent) choose(design$n_c + 3, 3) else 1)
    if (n_pairs > 1e7) {
      stop(sprintf(paste("Enumerating all outcomes would need %s pairs of",
                         "treatment and control outcomes (at most 1e7); use",
                         "enumerate = FALSE in design_binary2()."),
                   format(n_pairs, big.mark = ",")), call. = FALSE)
    }
    x_t <- outcomes_binary2(design$n_t)
    alpha_t <- posterior_binary2(design$prior_t, x_t, design$external_t)
    w_t <- vapply(seq_len(n_scen), function(s) {
      dmultinom_counts(x_t, cells_t[s, ])
    }, numeric(nrow(x_t)))
    if (concurrent) {
      x_c <- outcomes_binary2(design$n_c)
      alpha_c <- posterior_binary2(design$prior_c, x_c, design$external_c)
      w_c <- vapply(seq_len(n_scen), function(s) {
        dmultinom_counts(x_c, cells_c[s, ])
      }, numeric(nrow(x_c)))
    } else {
      alpha_c <- fixed_c
      w_c <- matrix(1, 1L, n_scen)
    }
    w_t <- matrix(w_t, ncol = n_scen)
    w_c <- matrix(w_c, ncol = n_scen)
    w_t[, bad] <- NA_real_
    n_t <- nrow(alpha_t)
    n_c <- nrow(alpha_c)
    p <- region_prob_binary2(crit, alpha_t, alpha_c,
                             rep(seq_len(n_t), times = n_c),
                             rep(seq_len(n_c), each = n_t),
                             c(sets$go, sets$nogo), design$nMC)
    s <- split_sets(p, n_go, dim = c(n_t, n_c))
    return(list(p_go = s$go, p_nogo = s$nogo, per_scenario = FALSE,
                weights_t = w_t, weights_c = w_c,
                method = paste("exact enumeration of the outcomes;",
                               prob_label)))
  }
  nsim <- as.integer(round(nsim))
  u_t <- matrix(stats::runif(nsim * design$n_t), nrow = nsim)
  if (concurrent) {
    u_c <- matrix(stats::runif(nsim * design$n_c), nrow = nsim)
  }
  # Counts of all simulated studies, scenarios stacked (NA for scenarios
  # that are not evaluable)
  sim_t <- matrix(NA_integer_, nsim * n_scen, 4L)
  sim_c <- matrix(0L, nsim * n_scen, 4L)
  for (s in which(!bad)) {
    rows <- (s - 1L) * nsim + seq_len(nsim)
    sim_t[rows, ] <- simulate_counts_binary2(u_t, cells_t[s, ])
    if (concurrent) {
      sim_c[rows, ] <- simulate_counts_binary2(u_c, cells_c[s, ])
    }
  }
  ok <- !is.na(sim_t[, 1L])
  key_t <- outcome_key(sim_t[ok, , drop = FALSE])
  key_c <- outcome_key(sim_c[ok, , drop = FALSE])
  ut <- unique(key_t)
  uc <- unique(key_c)
  it <- match(key_t, ut)
  ic <- match(key_c, uc)
  pair <- unique(cbind(it, ic))
  first_t <- match(ut, key_t)
  first_c <- match(uc, key_c)
  alpha_t <- posterior_binary2(design$prior_t,
                               sim_t[ok, , drop = FALSE][first_t, ,
                                                          drop = FALSE],
                               design$external_t)
  alpha_c <- if (concurrent) {
    posterior_binary2(design$prior_c,
                      sim_c[ok, , drop = FALSE][first_c, , drop = FALSE],
                      design$external_c)
  } else {
    fixed_c
  }
  p <- region_prob_binary2(crit, alpha_t, alpha_c, pair[, 1L], pair[, 2L],
                           c(sets$go, sets$nogo), design$nMC)
  idx <- match(paste(it, ic), paste(pair[, 1L], pair[, 2L]))
  full <- matrix(NA_real_, nsim * n_scen, ncol(p))
  full[ok, ] <- p[idx, , drop = FALSE]
  s <- split_sets(full, n_go, dim = c(nsim, n_scen))
  list(p_go = s$go, p_nogo = s$nogo, per_scenario = TRUE,
       weights = matrix(1 / nsim, nsim, n_scen),
       method = sprintf("simulation of %d studies; %s", nsim, prob_label))
}
