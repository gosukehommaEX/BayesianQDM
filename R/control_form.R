#' Form of the control group of a design
#'
#' Checks the combination of the control group arguments of a design and
#' returns the form of the control group.
#'
#' @param n_c Number of subjects in the concurrent control group.
#' @param external_c External control data or \code{NULL}.
#' @param hypothetical_c Hypothetical control or \code{NULL}.
#' @param fixed_c Fixed control value(s) or \code{NULL}.
#'
#' @return \code{"concurrent"}, \code{"external"}, \code{"hypothetical"} or
#'   \code{"fixed"}.
#' @noRd
control_form <- function(n_c, external_c, hypothetical_c, fixed_c = NULL) {
  if (n_c > 0) {
    if (!is.null(hypothetical_c)) {
      stop("'hypothetical_c' can be used only when 'n_c' is 0.",
           call. = FALSE)
    }
    if (!is.null(fixed_c)) {
      stop("'fixed_c' can be used only when 'n_c' is 0.", call. = FALSE)
    }
    return("concurrent")
  }
  given <- c(external = !is.null(external_c),
             hypothetical = !is.null(hypothetical_c),
             fixed = !is.null(fixed_c))
  if (!any(given)) {
    stop(paste("When 'n_c' is 0, specify 'external_c' (external control",
               "data), 'hypothetical_c' (hypothetical control) or 'fixed_c'",
               "(fixed control value)."), call. = FALSE)
  }
  if (sum(given) > 1L) {
    stop(paste("When 'n_c' is 0, specify only one of 'external_c',",
               "'hypothetical_c' and 'fixed_c'."), call. = FALSE)
  }
  names(given)[given]
}
