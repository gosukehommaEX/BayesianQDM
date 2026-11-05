#' Conduct of a time-to-event study as R code
#'
#' @param x An object of class \code{"qdm_conduct"}.
#'
#' @return A call from \code{code_node()} of \code{conduct_tte()}, without
#'   the arguments equal to their defaults.
#' @noRd
code_conduct <- function(x) {
  args <- list(accrual_time = code_num(x$accrual_time))
  if (!identical(x$accrual_rate, rep(1, length(x$accrual_time) - 1L))) {
    args$accrual_rate <- code_num(x$accrual_rate)
  }
  if (!identical(x$dropout_time, 0)) {
    args$dropout_time <- code_num(x$dropout_time)
  }
  if (!identical(x$dropout_hazard, 0)) {
    args$dropout_hazard <- code_num(x$dropout_hazard)
  }
  if (is.null(x$cut_time)) {
    args$cut_events <- code_num(x$cut_events, keep_type = FALSE)
  } else {
    args$cut_time <- code_num(x$cut_time)
  }
  code_node("conduct_tte", args)
}
