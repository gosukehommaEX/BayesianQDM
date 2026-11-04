#' Name of the design function of an endpoint
#'
#' @param endpoint Type of endpoint: \code{"binary"}, \code{"continuous"},
#'   \code{"binary2"}, \code{"continuous2"} or \code{"tte"}.
#'
#' @return The name of the \code{design_*()} function of the endpoint.
#' @noRd
design_function <- function(endpoint) {
  switch(endpoint,
         binary = "design_binary",
         continuous = "design_continuous",
         binary2 = "design_binary2",
         continuous2 = "design_continuous2",
         tte = "design_tte",
         stop("Unknown endpoint.", call. = FALSE))
}
