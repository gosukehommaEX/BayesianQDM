#' Form of the control group of a design
#'
#' Checks the combination of the control group arguments of a design and
#' returns the form of the control group.
#'
#' @param n_c Number of subjects in the concurrent control group.
#' @param external_c External control data or \code{NULL}.
#' @param hypothetical_c Hypothetical control or \code{NULL}.
#'
#' @return \code{"concurrent"}, \code{"external"} or \code{"hypothetical"}.
#' @noRd
control_form <- function(n_c, external_c, hypothetical_c) {
  if (n_c > 0) {
    if (!is.null(hypothetical_c)) {
      stop("'hypothetical_c' can be used only when 'n_c' is 0.",
           call. = FALSE)
    }
    return("concurrent")
  }
  if (is.null(external_c) && is.null(hypothetical_c)) {
    stop(paste("When 'n_c' is 0, specify 'external_c' (external control",
               "data) or 'hypothetical_c' (hypothetical control)."),
         call. = FALSE)
  }
  if (!is.null(external_c) && !is.null(hypothetical_c)) {
    stop(paste("When 'n_c' is 0, specify only one of 'external_c' and",
               "'hypothetical_c'."), call. = FALSE)
  }
  if (is.null(external_c)) "hypothetical" else "external"
}
