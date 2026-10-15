#' Simulated outcomes of a time-to-event design and their probabilities
#'
#' Simulates \code{nsim} studies under each scenario as specified by the
#' conduct of the design (\code{cpp_tte_oc_data()}), with exponential event
#' times whose hazards are \eqn{\log(2) / M} for the median survival times
#' \eqn{M} of the scenario, and computes the probabilities of the Go and NoGo
#' regions for each simulated study. The same random numbers (enrollment
#' times, standard exponential numbers for the event times and dropout times)
#' are used for all scenarios.
#'
#' @param design An object of class \code{"qdm_design_tte"}.
#' @param scenarios Data frame of scenarios (checked).
#' @param nsim Number of simulated studies.
#'
#' @return A list with elements \code{p_go} and \code{p_nogo} (matrices with
#'   one row per simulated study and one column per scenario),
#'   \code{per_scenario} (\code{TRUE}), \code{weights} (matrix of the same
#'   size with all elements \code{1 / nsim}) and \code{method}.
#' @noRd
oc_sample_tte <- function(design, scenarios, nsim) {
  if (is.null(design$conduct)) {
    stop(paste("The operating characteristics need the conduct of the",
               "study: give 'conduct' in design_tte()."), call. = FALSE)
  }
  nsim <- as.integer(round(nsim))
  n_scen <- nrow(scenarios)
  concurrent <- design$control == "concurrent"
  ca <- conduct_args(design$conduct)
  hazard_t <- log(2) / scenarios$median_t
  hazard_c <- if (concurrent) log(2) / scenarios$median_c else rep(1, n_scen)
  sim <- cpp_tte_oc_data(nsim, design$n_t, if (concurrent) design$n_c else 0L,
                         hazard_t, hazard_c, ca$acc_time, ca$acc_cum,
                         ca$drop_start, ca$drop_hazard, ca$cut_time,
                         ca$cut_events)
  post_t <- posterior_tte(design$prior_t, sim[, 1L], sim[, 2L],
                          external = design$external_t)
  post_c <- switch(
    design$control,
    concurrent = posterior_tte(design$prior_c, sim[, 3L], sim[, 4L],
                               external = design$external_c),
    fixed = NULL,
    posterior_tte(design$prior_c, 0, 0, external = design$external_c,
                  hypothetical = design$hypothetical_c)
  )
  p <- region_prob_tte(design, post_t, post_c)
  list(p_go = matrix(p$go, nsim, n_scen),
       p_nogo = matrix(p$nogo, nsim, n_scen), per_scenario = TRUE,
       weights = matrix(1 / nsim, nsim, n_scen),
       method = sprintf("simulation of %d studies; probabilities by %s",
                        nsim, tte_method_label(design)))
}
