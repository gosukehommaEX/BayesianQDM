#' Check the conduct of a study
#'
#' @param x Object to check: an object of class \code{"qdm_conduct"} (or
#'   \code{NULL} if \code{allow_null} is \code{TRUE}).
#' @param name Name of the argument, used in error messages.
#' @param n Number of subjects of the study; an analysis at an event must
#'   not need more events than subjects.
#' @param allow_null Logical; whether \code{NULL} is accepted.
#'
#' @return \code{x}, invisibly.
#' @noRd
check_conduct <- function(x, name, n, allow_null = FALSE) {
  if (is.null(x) && allow_null) {
    return(invisible(x))
  }
  if (!inherits(x, "qdm_conduct")) {
    stop(sprintf("'%s' must be created by conduct_tte().", name),
         call. = FALSE)
  }
  if (!is.null(x$cut_events) && x$cut_events > n) {
    stop(sprintf(paste("The analysis of '%s' is at event %d, but the study",
                       "has only %d subjects."), name, x$cut_events, n),
         call. = FALSE)
  }
  invisible(x)
}
