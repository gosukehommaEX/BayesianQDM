#' Print a decision criterion
#'
#' @param x An object of class \code{"qdm_criterion"}.
#' @param ... Passed to \code{\link{format.qdm_criterion}}.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' print(criterion_predictive(theta_null = 0, m_t = 50, gamma_go = 0.9,
#'                            gamma_nogo = 0.7))
#' @export
print.qdm_criterion <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}
