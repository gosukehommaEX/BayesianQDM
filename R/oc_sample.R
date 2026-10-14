#' Outcomes of a design and their probabilities
#'
#' Dispatches to the function for the endpoint of the design. For one binary
#' endpoint all outcomes are enumerated and \code{nsim} is not used; for
#' continuous endpoints \code{nsim} studies are simulated; for two binary
#' endpoints the outcomes are enumerated or simulated as specified by the
#' design.
#'
#' @param design An object of class \code{"qdm_design"}.
#' @param scenarios Data frame of scenarios (checked).
#' @param nsim Number of simulated data sets (simulation-based endpoints).
#'
#' @return A list with elements \code{p_go} and \code{p_nogo} (vectors or
#'   matrices, or lists of two for the marginal rules of two endpoints),
#'   \code{per_scenario} (whether they have one column per scenario), the
#'   probabilities of the outcomes (\code{weights}, or \code{weights_t} and
#'   \code{weights_c} for enumerated pairs of outcomes; see
#'   \code{outcome_prob()}) and \code{method}.
#' @noRd
oc_sample <- function(design, scenarios, nsim) {
  switch(design$endpoint,
         binary = oc_sample_binary(design, scenarios),
         continuous = oc_sample_continuous(design, scenarios, nsim),
         binary2 = oc_sample_binary2(design, scenarios, nsim),
         continuous2 = oc_sample_continuous2(design, scenarios, nsim))
}
