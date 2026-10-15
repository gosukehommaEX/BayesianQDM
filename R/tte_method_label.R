#' Description of the computation of the probabilities of a time-to-event
#' design
#'
#' @param design An object of class \code{"qdm_design_tte"}.
#'
#' @return A character string.
#' @noRd
tte_method_label <- function(design) {
  if (design$criterion$type == "predictive") {
    return(sprintf("simulation of future studies (nMC = %d)", design$nMC))
  }
  if (design$control == "fixed") {
    "exact computation (fixed control)"
  } else if (identical(design$effect, "hazard_ratio")) {
    "exact computation (F distribution)"
  } else {
    "numerical integration"
  }
}
