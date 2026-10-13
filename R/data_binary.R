#' Binary endpoint data
#'
#' Summarizes the data of one group for a single binary endpoint as the number
#' of subjects and the number of responders. Used for observed proof-of-concept
#' data, external data and hypothetical control data.
#'
#' @param n Number of subjects (a whole number of at least 1).
#' @param y Number of responders (a whole number from 0 to \code{n}).
#'
#' @return An object of class \code{c("qdm_data_binary", "qdm_data")} with
#'   integer elements \code{n} and \code{y}.
#'
#' @family data
#' @examples
#' data_binary(n = 12, y = 7)
#' @export
data_binary <- function(n, y) {
  check_numeric(n, "n", lower = 1, integer = TRUE)
  check_numeric(y, "y", lower = 0, integer = TRUE)
  if (y > n) {
    stop("'y' must not exceed 'n'.", call. = FALSE)
  }
  structure(list(n = as.integer(round(n)), y = as.integer(round(y))),
            class = c("qdm_data_binary", "qdm_data"))
}
