#' Format a search of the sample size
#'
#' @param x An object of class \code{"qdm_n_search"}.
#' @param digits Number of significant digits.
#' @param ... Not used.
#'
#' @return A character vector of lines.
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' res <- qdm_find_n(des, data.frame(pi_t = 0.4, pi_c = 0.4),
#'                   data.frame(pi_t = 0.7, pi_c = 0.4), target_go = 0.1,
#'                   target_nogo = 0.1, n_t = seq(6, 20, by = 2))
#' format(res)
#' @export
format.qdm_n_search <- function(x, digits = 4L, ...) {
  num <- function(v) format_num(v, digits)
  size <- function(nt, nc) sprintf("n_t = %d, n_c = %d", nt, nc)
  head <- c(sprintf("Search of the sample size (%s)", x$method),
            sprintf(paste("  Targets: P(Go criterion met | scenario_go) <= %s",
                          "and P(NoGo criterion met | scenario_nogo) <= %s"),
                    num(x$target_go), num(x$target_nogo)))
  if (is.na(x$n_t)) {
    return(c(head, sprintf("  No candidate from %s to %s meets the targets.",
                           size(x$curve$n_t[1L], x$curve$n_c[1L]),
                           size(x$curve$n_t[nrow(x$curve)],
                                x$curve$n_c[nrow(x$curve)]))))
  }
  out <- c(head,
           sprintf("  Smallest candidate that meets the targets: %s",
                   size(x$n_t, x$n_c)),
           sprintf("    P(Go criterion met | scenario_go) = %s",
                   num(x$prob_go)),
           sprintf("    P(NoGo criterion met | scenario_nogo) = %s",
                   num(x$prob_nogo)))
  if (is.na(x$n_t_stable)) {
    out <- c(out, paste("  The largest candidate does not meet the targets",
                        "(the probabilities are not monotone in the sample",
                        "size)."))
  } else if (x$n_t_stable != x$n_t) {
    out <- c(out, sprintf(paste("  Smallest candidate from which all larger",
                                "candidates meet the targets: %s"),
                          size(x$n_t_stable, x$n_c_stable)))
  }
  out
}
