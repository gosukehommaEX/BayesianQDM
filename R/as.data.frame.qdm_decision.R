#' Decision as a data frame
#'
#' Converts a result of \code{\link{qdm_decide}} to a data frame with one
#' row, so that the decisions for several data sets (for example for
#' subgroups, or for the prior distributions of a sensitivity analysis) can
#' be combined with \code{rbind()} and written to a file.
#'
#' The columns are \code{decision}; \code{prob_go}, \code{gamma_go} and
#' \code{met_go} (the probability of the Go criterion, its threshold and
#' whether the criterion is met); and \code{prob_nogo}, \code{gamma_nogo}
#' and \code{met_nogo} for the NoGo criterion. For two endpoints with a
#' marginal rule (\code{rule = "and"} or \code{"or"}), the probabilities and
#' the thresholds have one column per endpoint (\code{prob_go1},
#' \code{prob_go2}, \code{gamma_go1}, \code{gamma_go2}, and so on), and
#' \code{met_go} and \code{met_nogo} combine the endpoints by the rule. For
#' two endpoints, the probabilities of the regions follow, in the columns
#' \code{R1} to \code{R9} (\code{R1} to \code{R4} for a predictive
#' criterion; see \code{\link{qdm_regions}}).
#'
#' @param x An object of class \code{"qdm_decision"}.
#' @param row.names \code{NULL} (the default) or the name of the row.
#' @param optional Not used.
#' @param ... Not used.
#' @param posterior \code{TRUE} to add the parameters of the posterior
#'   distributions (see \code{\link{qdm_decide}}) in columns named
#'   \code{post_t_} and \code{post_c_} followed by the names of the
#'   parameters, for example \code{post_t_a} and \code{post_t_b}.
#'
#' @return A data frame with one row.
#'
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' dec1 <- qdm_decide(des, data_t = data_binary(12, 9),
#'                    data_c = data_binary(12, 4))
#' dec2 <- qdm_decide(des, data_t = data_binary(12, 6),
#'                    data_c = data_binary(12, 5))
#' rbind(as.data.frame(dec1), as.data.frame(dec2))
#' as.data.frame(dec1, posterior = TRUE)
#' @export
as.data.frame.qdm_decision <- function(x, row.names = NULL, optional = FALSE,
                                       ..., posterior = FALSE) {
  check_flag(posterior, "posterior")
  crit <- x$design$criterion
  # One column for a single value, or one per endpoint
  columns <- function(name, v) {
    v <- unname(v)
    nms <- if (length(v) == 1L) name else paste0(name, seq_along(v))
    stats::setNames(as.list(v), nms)
  }
  cols <- c(list(decision = unname(x$decision)),
            columns("prob_go", x$prob_go), columns("gamma_go", crit$gamma_go),
            list(met_go = unname(x$met_go)),
            columns("prob_nogo", x$prob_nogo),
            columns("gamma_nogo", crit$gamma_nogo),
            list(met_nogo = unname(x$met_nogo)))
  if (!is.null(x$prob_regions)) {
    cols <- c(cols, as.list(x$prob_regions))
  }
  if (posterior) {
    cols <- c(cols,
              stats::setNames(as.list(unname(x$posterior_t)),
                              paste0("post_t_", names(x$posterior_t))),
              stats::setNames(as.list(unname(x$posterior_c)),
                              paste0("post_c_", names(x$posterior_c))))
  }
  out <- data.frame(cols, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(row.names)) {
    rownames(out) <- row.names
  }
  out
}
