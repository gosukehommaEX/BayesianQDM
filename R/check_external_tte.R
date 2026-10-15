#' Check external data for a time-to-event endpoint
#'
#' @param x Object to check: \code{NULL} or an object of class
#'   \code{"qdm_borrow"} with time-to-event data.
#' @param name Name of the argument, used in the error message.
#'
#' @return \code{x}, invisibly.
#' @noRd
check_external_tte <- function(x, name) {
  if (!is.null(x) && (!inherits(x, "qdm_borrow") ||
                      !inherits(x$data, "qdm_data_tte"))) {
    stop(sprintf("'%s' must be NULL or created by borrow() with data_tte().",
                 name), call. = FALSE)
  }
  invisible(x)
}
