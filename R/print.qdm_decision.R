#' Print a decision
#'
#' @param x An object of class \code{"qdm_decision"}.
#' @param ... Passed to \code{\link{format.qdm_decision}}.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_predictive(theta_null = 0, m_t = 50,
#'                                    gamma_go = 0.9, gamma_nogo = 0.7)
#' )
#' print(qdm_decide(des, data_binary(12, 9), data_binary(12, 4)))
#' @export
print.qdm_decision <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}
