#' Columns of the scenarios that set the treatment effect
#'
#' @param design An object of class \code{"qdm_design"}.
#'
#' @return The names of the columns of the treatment group that
#'   \code{qdm_scenarios()} computes from the effect \code{theta}.
#' @noRd
scenario_location <- function(design) {
  switch(design$endpoint,
         binary = "pi_t",
         continuous = "mu_t",
         tte = "median_t",
         binary2 = c("pi_t1", "pi_t2"),
         continuous2 = c("mu_t1", "mu_t2"))
}
