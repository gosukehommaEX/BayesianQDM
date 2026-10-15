#' Description of the method used to compute probabilities
#'
#' @param method One of \code{"NI"}, \code{"MM"} or \code{"MC"}.
#' @param nMC Number of Monte Carlo draws.
#' @param fixed Logical; if \code{TRUE} (fixed control), the probabilities
#'   are computed exactly whatever the method.
#'
#' @return A character string.
#' @noRd
method_label <- function(method, nMC, fixed = FALSE) {
  if (fixed) {
    return("exact computation (fixed control)")
  }
  switch(method,
         NI = "numerical integration",
         MM = "moment-matching approximation",
         MC = sprintf("Monte Carlo simulation (nMC = %d)", nMC))
}
