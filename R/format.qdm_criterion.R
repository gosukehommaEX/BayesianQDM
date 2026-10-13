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
  gamma_go <- format_num(x$gamma_go, digits)
  gamma_nogo <- format_num(x$gamma_nogo, digits)
  if (x$n_endpoints == 2L) {
    go <- sprintf("P(regions {%s} | data) >= %s",
                  paste(x$go_regions, collapse = ", "), gamma_go)
    nogo <- sprintf("P(regions {%s} | data) >= %s",
                    paste(x$nogo_regions, collapse = ", "), gamma_nogo)
  }
  if (x$type == "posterior") {
    header <- "Posterior probability criterion"
    tv <- format_num(x$theta_tv, digits)
    mav <- format_num(x$theta_mav, digits)
    if (x$n_endpoints == 1L) {
      go <- sprintf("P(theta > %s | data) >= %s", tv, gamma_go)
      nogo <- sprintf("P(theta <= %s | data) >= %s", mav, gamma_nogo)
      extra <- character(0)
    } else {
      extra <- sprintf("  TV = (%s), MAV = (%s)", tv, mav)
    }
  } else {
    header <- sprintf(
      "Posterior predictive probability criterion (m_t = %d, m_c = %d)",
      x$m_t, x$m_c
    )
    null <- format_num(x$theta_null, digits)
    if (x$n_endpoints == 1L) {
      go <- sprintf("P(future theta > %s | data) >= %s", null, gamma_go)
      nogo <- sprintf("P(future theta <= %s | data) >= %s", null, gamma_nogo)
      extra <- character(0)
    } else {
      extra <- sprintf("  NULL = (%s)", null)
    }
  }
  c(header, paste0("  Go:   ", go), paste0("  NoGo: ", nogo), extra)
}
