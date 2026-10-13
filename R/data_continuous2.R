#' Data for two continuous endpoints
#'
#' Summarises the data of one group for two continuous endpoints by their
#' sufficient statistics: the number of subjects, the vector of sample means
#' and the sample covariance matrix (with divisor \code{n - 1}).
#'
#' @param n Number of subjects (a whole number of at least 3).
#' @param mean Numeric vector of length 2, the sample means of endpoints 1 and
#'   2.
#' @param cov A 2 x 2 symmetric positive definite sample covariance matrix.
#'
#' @return An object of class \code{c("qdm_data_continuous2", "qdm_data")}
#'   with elements \code{n}, \code{mean} and \code{cov}.
#'
#' @family data
#' @examples
#' data_continuous2(n = 12, mean = c(3.2, 1.5),
#'                  cov = matrix(c(4, 1.2, 1.2, 2), nrow = 2))
#' @export
data_continuous2 <- function(n, mean, cov) {
  check_numeric(n, "n", lower = 3, integer = TRUE)
  check_numeric(mean, "mean", len = 2L)
  check_pd_matrix(cov, "cov")
  structure(list(n = as.integer(round(n)), mean = as.numeric(mean),
                 cov = unname(cov)),
            class = c("qdm_data_continuous2", "qdm_data"))
}
