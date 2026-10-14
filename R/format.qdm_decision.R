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
  num <- function(v) format_num(v, digits)
  describe <- function(name, lab, prob, gamma, met_all) {
    if (length(prob) == 1L) {
      return(sprintf("  %s %s = %s (threshold %s, %s)", name, lab, num(prob),
                     num(gamma), met(met_all)))
    }
    c(sprintf("  %s %s (rule \"%s\")", name, met(met_all), crit$rule),
      sprintf("    %s = %s (threshold %s, %s)", lab, vapply(prob, num, ""),
              vapply(gamma, num, ""),
              vapply(meets_threshold(prob, gamma), met, "")))
  }
  out <- c(sprintf("Decision: %s", x$decision),
           describe("Go:  ", labels$go, x$prob_go, crit$gamma_go, x$met_go),
           describe("NoGo:", labels$nogo, x$prob_nogo, crit$gamma_nogo,
                    x$met_nogo))
  if (!is.null(x$prob_regions)) {
    out <- c(out, sprintf("Region probabilities: %s",
                          paste(names(x$prob_regions),
                                vapply(x$prob_regions, num, ""),
                                sep = " = ", collapse = ", ")))
  }
  c(out, "Posterior distributions",
    sprintf("  Treatment: %s", format_posterior(x$design, x$posterior_t,
                                                digits)),
    sprintf("  Control:   %s", format_posterior(x$design, x$posterior_c,
                                                digits)))
}
