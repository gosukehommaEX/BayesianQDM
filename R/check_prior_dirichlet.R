#' Check that an argument is a Dirichlet prior
#'
#' @param x Object to check.
#' @param name Name of the argument, used in the error message.
#'
#' @return \code{x}, invisibly.
#' @noRd
check_prior_dirichlet <- function(x, name) {
  if (!inherits(x, "qdm_prior_dirichlet")) {
    stop(sprintf("'%s' must be created by prior_dirichlet().", name),
         call. = FALSE)
  }
  invisible(x)
}
