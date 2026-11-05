#' Search of the sample size as a data frame
#'
#' Converts a result of \code{\link{qdm_find_n}} to a data frame with one
#' row: the chosen sample sizes, the probabilities of meeting the criteria
#' with them, the targets and the smallest sample sizes from which all
#' larger candidates meet the targets. The probabilities for all candidates
#' are in the element \code{curve} of the result, which is a data frame.
#'
#' @param x An object of class \code{"qdm_n_search"}.
#' @param row.names \code{NULL} (the default) or the name of the row.
#' @param optional Not used.
#' @param ... Not used.
#'
#' @return A data frame with one row and the columns \code{n_t}, \code{n_c},
#'   \code{prob_go}, \code{target_go}, \code{prob_nogo}, \code{target_nogo},
#'   \code{n_t_stable} and \code{n_c_stable} (see \code{\link{qdm_find_n}}).
#'
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
#' as.data.frame(res)
#' @export
as.data.frame.qdm_n_search <- function(x, row.names = NULL, optional = FALSE,
                                       ...) {
  out <- data.frame(n_t = x$n_t, n_c = x$n_c, prob_go = x$prob_go,
                    target_go = x$target_go, prob_nogo = x$prob_nogo,
                    target_nogo = x$target_nogo, n_t_stable = x$n_t_stable,
                    n_c_stable = x$n_c_stable)
  if (!is.null(row.names)) {
    rownames(out) <- row.names
  }
  out
}
