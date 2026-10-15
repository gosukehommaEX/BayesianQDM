#' Check external data for a continuous endpoint
#'
#' @param x Object to check: \code{NULL} or an object of class
#'   \code{"qdm_borrow"} with continuous data.
#' @param name Name of the argument, used in the error message.
#'
#' @return \code{x}, invisibly.
#' @noRd
check_external_continuous <- function(x, name) {
  if (!is.null(x) && (!inherits(x, "qdm_borrow") ||
                      !inherits(x$data, "qdm_data_continuous"))) {
    stop(sprintf(paste("'%s' must be NULL or created by borrow() with",
                       "data_continuous()."), name), call. = FALSE)
  }
  invisible(x)
}
