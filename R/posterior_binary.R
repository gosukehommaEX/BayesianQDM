#' Posterior Beta parameters for a binary endpoint
#'
#' Combines a Beta prior, external data borrowed through a power prior,
#' hypothetical control data and the current data into the parameters of the
#' Beta posterior distribution of the response rate. With prior
#' Beta(a_p, b_p), external data (n_e, y_e) with weight a_0, hypothetical
#' data (n_h, y_h) and current data (n, y), the posterior is Beta(a, b) with
#' a = a_p + a_0 y_e + y_h + y and
#' b = b_p + a_0 (n_e - y_e) + (n_h - y_h) + (n - y).
#'
#' @param prior An object of class \code{"qdm_prior_beta"}.
#' @param n Number of subjects in the current data (0 if there are none).
#' @param y Numeric vector of numbers of responders in the current data.
#' @param external An object of class \code{"qdm_borrow"} or \code{NULL}.
#' @param hypothetical An object of class \code{"qdm_hypothetical"} with
#'   binary data, or \code{NULL}.
#'
#' @return A list with numeric vectors \code{a} and \code{b} of the same
#'   length as \code{y}.
#' @noRd
posterior_binary <- function(prior, n, y, external = NULL,
                             hypothetical = NULL) {
  a <- prior$a + y
  b <- prior$b + (n - y)
  if (!is.null(external)) {
    a <- a + external$a0 * external$data$y
    b <- b + external$a0 * (external$data$n - external$data$y)
  }
  if (!is.null(hypothetical)) {
    a <- a + hypothetical$data$y
    b <- b + (hypothetical$data$n - hypothetical$data$y)
  }
  list(a = a, b = b)
}
