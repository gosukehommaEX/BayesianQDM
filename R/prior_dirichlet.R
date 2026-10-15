#' Dirichlet prior for two binary endpoints
#'
#' Creates a Dirichlet prior distribution for the four joint response
#' probabilities of two binary endpoints in one group.
#'
#' The four cells are ordered as \code{"00"}, \code{"01"}, \code{"10"} and
#' \code{"11"}, where the first digit is the response on endpoint 1 and the
#' second digit is the response on endpoint 2 (1 = response, 0 = no
#' response). The marginal response rates are therefore
#' \eqn{\pi_1 = p_{10} + p_{11}} and \eqn{\pi_2 = p_{01} + p_{11}}.
#' With \code{alpha = rep(0.25, 4)}, each marginal response rate has a
#' Beta(0.5, 0.5) prior distribution.
#'
#' @param alpha A numeric vector of four positive concentration parameters in
#'   the cell order \code{"00"}, \code{"01"}, \code{"10"}, \code{"11"}.
#'
#' @return An object of class \code{c("qdm_prior_dirichlet", "qdm_prior")}
#'   with element \code{alpha}, a named numeric vector of length 4.
#'
#' @family priors
#' @examples
#' prior_dirichlet(rep(0.25, 4))
#' @export
prior_dirichlet <- function(alpha) {
  check_numeric(alpha, "alpha", len = 4L, lower = 0, lower_open = TRUE)
  alpha <- stats::setNames(as.numeric(alpha), c("00", "01", "10", "11"))
  structure(list(alpha = alpha),
            class = c("qdm_prior_dirichlet", "qdm_prior"))
}
