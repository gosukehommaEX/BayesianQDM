#' Continuous endpoint data
#'
#' Summarizes the data of one group for a single continuous endpoint by its
#' sufficient statistics: the number of subjects, the sample mean and the
#' sample standard deviation (with divisor \code{n - 1}).
#'
#' @param n Number of subjects (a whole number of at least 2).
#' @param mean Sample mean.
#' @param sd Sample standard deviation (positive).
#'
#' @return An object of class \code{c("qdm_data_continuous", "qdm_data")} with
#'   elements \code{n}, \code{mean} and \code{sd}.
#'
#' @family data
#' @examples
#' data_continuous(n = 12, mean = 3.2, sd = 2.0)
#' @export
data_continuous <- function(n, mean, sd) {
  check_numeric(n, "n", lower = 2, integer = TRUE)
  check_numeric(mean, "mean")
  check_numeric(sd, "sd", lower = 0, lower_open = TRUE)
  structure(list(n = as.integer(round(n)), mean = mean, sd = sd),
            class = c("qdm_data_continuous", "qdm_data"))
}
