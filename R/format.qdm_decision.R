#' Format a decision
#'
#' @param x An object of class \code{"qdm_decision"}.
#' @param digits Number of significant digits for numeric values.
#' @param ... Not used.
#'
#' @return A character vector with one element per line of the description.
#'
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' format(qdm_decide(des, data_binary(12, 9), data_binary(12, 4)))
#' @export
format.qdm_decision <- function(x, digits = 4L, ...) {
  crit <- x$design$criterion
  labels <- event_labels(crit, digits)
  met <- function(m) if (m) "met" else "not met"
  c(
    sprintf("Decision: %s", x$decision),
    sprintf("  Go:   %s = %s (threshold %s, %s)", labels[["go"]],
            format_num(x$prob_go, digits), format_num(crit$gamma_go, digits),
            met(x$met_go)),
    sprintf("  NoGo: %s = %s (threshold %s, %s)", labels[["nogo"]],
            format_num(x$prob_nogo, digits),
            format_num(crit$gamma_nogo, digits), met(x$met_nogo)),
    "Posterior distributions",
    sprintf("  Treatment: %s", format_posterior(x$design, x$posterior_t,
                                                digits)),
    sprintf("  Control:   %s", format_posterior(x$design, x$posterior_c,
                                                digits))
  )
}
