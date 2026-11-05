#' Print R code
#'
#' @param x An object of class \code{"qdm_code"}, from
#'   \code{\link{show_code}}.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' print(show_code(des))
#' @export
print.qdm_code <- function(x, ...) {
  cat(x, sep = "\n")
  invisible(x)
}
