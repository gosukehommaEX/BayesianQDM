#' Check that an argument is a Beta prior
#'
#' @param x Object to check.
#' @param name Name of the argument, used in the error message.
#'
#' @return \code{x}, invisibly.
#' @noRd
check_prior_beta <- function(x, name) {
  if (!inherits(x, "qdm_prior_beta")) {
    stop(sprintf("'%s' must be created by prior_beta().", name),
         call. = FALSE)
  }
  invisible(x)
}
