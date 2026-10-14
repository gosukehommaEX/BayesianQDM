#' Check external data for two binary endpoints
#'
#' @param x Object to check: \code{NULL} or an object of class
#'   \code{"qdm_borrow"} with data for two binary endpoints.
#' @param name Name of the argument, used in the error message.
#'
#' @return \code{x}, invisibly.
#' @noRd
check_external_binary2 <- function(x, name) {
  if (!is.null(x) && (!inherits(x, "qdm_borrow") ||
                      !inherits(x$data, "qdm_data_binary2"))) {
    stop(sprintf(paste("'%s' must be NULL or created by borrow() with",
                       "data_binary2()."), name), call. = FALSE)
  }
  invisible(x)
}
