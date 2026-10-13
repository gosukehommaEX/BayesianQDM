#' External data borrowed through a power prior
#'
#' Specifies external (historical) data for one group together with the
#' power prior weight \eqn{a_0}. The likelihood of the external data is raised
#' to the power \eqn{a_0}, so that \eqn{a_0 = 0} ignores the external data and
#' \eqn{a_0 = 1} pools it with the current data. For the conjugate models used
#' in this package the resulting posterior distribution has a closed form, so
#' no Markov chain Monte Carlo sampling is needed.
#'
#' @param data External data created by \code{\link{data_binary}},
#'   \code{\link{data_continuous}}, \code{\link{data_binary2}} or
#'   \code{\link{data_continuous2}}.
#' @param a0 Power prior weight, a number from 0 to 1.
#'
#' @return An object of class \code{"qdm_borrow"} with elements \code{data} and
#'   \code{a0}.
#'
#' @examples
#' borrow(data_binary(n = 20, y = 8), a0 = 0.5)
#' @export
borrow <- function(data, a0) {
  if (!inherits(data, "qdm_data")) {
    stop("'data' must be created by a data_*() function, such as data_binary().",
         call. = FALSE)
  }
  check_numeric(a0, "a0", lower = 0, upper = 1)
  structure(list(data = data, a0 = a0), class = "qdm_borrow")
}
