#' Description of the method used to compute probabilities
#'
#' @param method One of \code{"NI"}, \code{"MM"} or \code{"MC"}.
#' @param nMC Number of Monte Carlo draws.
#'
#' @return A character string.
#' @noRd
method_label <- function(method, nMC) {
  switch(method,
         NI = "numerical integration",
         MM = "moment-matching approximation",
         MC = sprintf("Monte Carlo simulation (nMC = %d)", nMC))
}
