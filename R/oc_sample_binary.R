#' Outcomes of a binary design and their probabilities
#'
#' Enumerates all possible outcomes of the proof-of-concept study of a binary
#' design, computes the probabilities of the Go and NoGo regions for each
#' outcome, and the probability of each outcome under each scenario.
#'
#' @param design An object of class \code{"qdm_design_binary"}.
#' @param scenarios Data frame of scenarios (checked).
#'
#' @return A list with elements \code{p_go} and \code{p_nogo} (vectors over
#'   the outcomes), \code{per_scenario} (\code{FALSE}), \code{weights}
#'   (matrix of outcome probabilities with one column per scenario) and
#'   \code{method}.
#' @noRd
oc_sample_binary <- function(design, scenarios) {
  n_t <- design$n_t
  y_t <- 0:n_t
  post_t <- posterior_binary(design$prior_t, n_t, y_t,
                             external = design$external_t)
  # Binomial probabilities: one row per outcome, one column per scenario
  w_t <- vapply(scenarios$pi_t, function(p) stats::dbinom(y_t, n_t, p),
                numeric(n_t + 1L))
  if (design$control == "concurrent") {
    n_c <- design$n_c
    y_c <- 0:n_c
    post_c <- posterior_binary(design$prior_c, n_c, y_c,
                               external = design$external_c)
    w_c <- vapply(scenarios$pi_c, function(p) stats::dbinom(y_c, n_c, p),
                  numeric(n_c + 1L))
    # Pairs (y_t, y_c) with y_t varying fastest
    i_t <- rep(seq_along(y_t), times = n_c + 1L)
    i_c <- rep(seq_along(y_c), each = n_t + 1L)
    weights <- w_t[i_t, , drop = FALSE] * w_c[i_c, , drop = FALSE]
    p <- region_prob_binary(design$criterion, post_t$a[i_t], post_t$b[i_t],
                            post_c$a[i_c], post_c$b[i_c])
  } else {
    post_c <- if (design$control == "fixed") {
      NULL
    } else {
      posterior_binary(design$prior_c, 0, 0, external = design$external_c,
                       hypothetical = design$hypothetical_c)
    }
    weights <- w_t
    p <- region_prob_binary(design$criterion, post_t$a, post_t$b, post_c$a,
                            post_c$b, fixed = design$fixed_c)
  }
  list(p_go = p$go, p_nogo = p$nogo, per_scenario = FALSE,
       weights = weights, method = "exact enumeration of the outcomes")
}
