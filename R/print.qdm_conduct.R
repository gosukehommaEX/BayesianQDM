#' Print the conduct of a study
#'
#' @param x An object of class \code{"qdm_conduct"}.
#' @param ... Passed to \code{\link{format.qdm_conduct}}.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' conduct_tte(accrual_time = c(0, 12), cut_events = 20)
#' @export
print.qdm_conduct <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}
