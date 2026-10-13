#' Columns of the scenarios of a design
#'
#' @param design An object of class \code{"qdm_design"}.
#'
#' @return A list with elements \code{names} (the required columns),
#'   \code{lower} and \code{upper} (bounds of their values).
#' @noRd
scenario_columns <- function(design) {
  switch(design$endpoint,
         binary = list(names = if (design$control == "concurrent") {
           c("pi_t", "pi_c")
         } else {
           "pi_t"
         }, lower = 0, upper = 1))
}
