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
  marginal <- x$n_endpoints == 2L && x$rule != "joint"
  join <- if (x$rule == "or") " or " else " and "
  criteria <- function(lab, gamma) {
    paste(sprintf("%s >= %s", lab,
                  vapply(gamma, format_num, character(1), digits = digits)),
          collapse = join)
  }
  extra <- character(0)
  if (x$type == "posterior") {
    header <- "Posterior probability criterion"
    if (x$n_endpoints == 2L && !marginal) {
      extra <- sprintf("  TV = (%s), MAV = (%s)",
                       format_num(x$theta_tv, digits),
                       format_num(x$theta_mav, digits))
    }
  } else {
    header <- sprintf(
      "Posterior predictive probability criterion (m_t = %d, m_c = %d)",
      x$m_t, x$m_c
    )
    if (x$n_endpoints == 2L && !marginal) {
      extra <- sprintf("  NULL = (%s)", format_num(x$theta_null, digits))
    }
  }
  if (marginal) {
    header <- sprintf("%s, marginal rule \"%s\"", header, x$rule)
  }
  if (x$n_endpoints == 2L && !marginal && any(x$direction == "less")) {
    # The regions are defined on the favorable side of each endpoint
    extra <- c(extra, sprintf("  Favorable direction: theta1 %s, theta2 %s",
                              x$direction[1L], x$direction[2L]))
  }
  c(header, paste0("  Go:   ", criteria(labels$go, x$gamma_go)),
    paste0("  NoGo: ", criteria(labels$nogo, x$gamma_nogo)), extra)
}
