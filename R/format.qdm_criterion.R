#' Format a decision criterion
#'
#' @param x An object of class \code{"qdm_criterion"}.
#' @param digits Number of significant digits for numeric values.
#' @param ... Not used.
#'
#' @return A character vector with one element per line of the description.
#'
#' @examples
#' format(criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                            gamma_go = 0.5, gamma_nogo = 0.8))
#' @export
format.qdm_criterion <- function(x, digits = 4L, ...) {
  labels <- event_labels(x, digits)
  go <- sprintf("%s >= %s", labels[["go"]], format_num(x$gamma_go, digits))
  nogo <- sprintf("%s >= %s", labels[["nogo"]],
                  format_num(x$gamma_nogo, digits))
  extra <- character(0)
  if (x$type == "posterior") {
    header <- "Posterior probability criterion"
    if (x$n_endpoints == 2L) {
      extra <- sprintf("  TV = (%s), MAV = (%s)",
                       format_num(x$theta_tv, digits),
                       format_num(x$theta_mav, digits))
    }
  } else {
    header <- sprintf(
      "Posterior predictive probability criterion (m_t = %d, m_c = %d)",
      x$m_t, x$m_c
    )
    if (x$n_endpoints == 2L) {
      extra <- sprintf("  NULL = (%s)", format_num(x$theta_null, digits))
    }
  }
  c(header, paste0("  Go:   ", go), paste0("  NoGo: ", nogo), extra)
}
