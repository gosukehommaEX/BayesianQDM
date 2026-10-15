#' Check a fixed control value
#'
#' @param x Object to check: \code{NULL} or the fixed value(s) of the
#'   parameter of the control group.
#' @param n_endpoints Number of endpoints (the required length of \code{x}).
#' @param probability Logical; if \code{TRUE}, the values are probabilities
#'   (response rates or milestone survival probabilities) and must be from 0
#'   to 1.
#'
#' @return \code{x} as a numeric vector, or \code{NULL}.
#' @noRd
check_fixed_c <- function(x, n_endpoints, probability) {
  if (is.null(x)) {
    return(NULL)
  }
  if (probability) {
    check_numeric(x, "fixed_c", len = n_endpoints, lower = 0, upper = 1)
  } else {
    check_numeric(x, "fixed_c", len = n_endpoints)
  }
  as.numeric(x)
}
