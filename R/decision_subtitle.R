#' Probabilities of the criteria of a decision as text
#'
#' @param x An object of class \code{"qdm_decision"}.
#' @param digits Number of decimal places of the probabilities.
#'
#' @return A character string with one line for the Go criterion and one for
#'   the NoGo criterion: the probability, the threshold and whether the
#'   criterion is met (for each endpoint with a marginal rule).
#' @noRd
decision_subtitle <- function(x, digits = 3L) {
  crit <- x$design$criterion
  labels <- event_labels(crit, digits)
  prob <- function(p) formatC(p, format = "f", digits = digits)
  met <- function(m) if (m) "met" else "not met"
  line <- function(name, lab, p, gamma, m) {
    p <- unname(p)
    if (length(p) == 1L) {
      return(sprintf("%s: %s = %s (threshold %s, %s)", name, lab, prob(p),
                     format_num(gamma, digits), met(m)))
    }
    sprintf("%s (rule \"%s\", %s): %s", name, crit$rule, met(m),
            paste(sprintf("%s = %s (threshold %s)", lab, prob(p),
                          vapply(gamma, format_num, character(1),
                                 digits = digits)), collapse = "; "))
  }
  paste(line("Go", labels$go, x$prob_go, crit$gamma_go, x$met_go),
        line("NoGo", labels$nogo, x$prob_nogo, crit$gamma_nogo, x$met_nogo),
        sep = "\n")
}
