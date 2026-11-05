#' Range of the treatment effect of one endpoint
#'
#' @param design A design for one endpoint.
#'
#' @return The lower and upper limits of the treatment effect (possibly
#'   infinite).
#' @noRd
effect_support <- function(design) {
  fixed <- design$fixed_c
  if (design$endpoint == "continuous") {
    return(c(-Inf, Inf))
  }
  if (identical(design$effect, "hazard_ratio")) {
    return(c(0, Inf))
  }
  # Difference of probabilities, or a probability minus a fixed value
  if (is.null(fixed)) c(-1, 1) else c(-fixed, 1 - fixed)
}
