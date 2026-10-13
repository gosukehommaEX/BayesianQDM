#' Outcomes of a design and their probabilities
#'
#' Dispatches to the function for the endpoint of the design. For binary
#' endpoints all outcomes are enumerated and \code{nsim} is not used.
#'
#' @param design An object of class \code{"qdm_design"}.
#' @param scenarios Data frame of scenarios (checked).
#' @param nsim Number of simulated data sets (simulation-based endpoints).
#'
#' @return A list with elements \code{p_go} and \code{p_nogo} (vectors over
#'   the outcomes, or matrices with one column per scenario), \code{weights}
#'   (matrix with one column per scenario) and \code{method}.
#' @noRd
oc_sample <- function(design, scenarios, nsim) {
  switch(design$endpoint,
         binary = oc_sample_binary(design, scenarios))
}
