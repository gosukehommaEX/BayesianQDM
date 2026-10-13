#' Check external data for a binary endpoint
#'
#' @param x Object to check: \code{NULL} or an object of class
#'   \code{"qdm_borrow"} with binary data.
#' @param name Name of the argument, used in the error message.
#'
#' @return \code{x}, invisibly.
#' @noRd
check_external_binary <- function(x, name) {
  if (!is.null(x) && (!inherits(x, "qdm_borrow") ||
                      !inherits(x$data, "qdm_data_binary"))) {
    stop(sprintf("'%s' must be NULL or created by borrow() with data_binary().",
                 name), call. = FALSE)
  }
  invisible(x)
}
