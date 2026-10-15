#' Print a design
#'
#' @param x An object of class \code{"qdm_design"}.
#' @param ... Passed to \code{\link{format.qdm_design}}.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 0,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8),
#'   hypothetical_c = hypothetical(data = data_binary(n = 12, y = 4))
#' )
#' print(des)
#' @export
print.qdm_design <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}
