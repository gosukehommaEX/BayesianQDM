#' Check that an argument is a Gamma prior
#'
#' @param x Object to check.
#' @param name Name of the argument, used in the error message.
#'
#' @return \code{x}, invisibly.
#' @noRd
check_prior_gamma <- function(x, name) {
  if (!inherits(x, "qdm_prior_gamma")) {
    stop(sprintf("'%s' must be created by prior_gamma().", name),
         call. = FALSE)
  }
  invisible(x)
}
