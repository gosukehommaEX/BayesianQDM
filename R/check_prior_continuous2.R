#' Check that an argument is a prior for two continuous endpoints
#'
#' @param x Object to check.
#' @param name Name of the argument, used in the error message.
#'
#' @return \code{x}, invisibly.
#' @noRd
check_prior_continuous2 <- function(x, name) {
  if (!inherits(x, c("qdm_prior_niw", "qdm_prior_vague"))) {
    stop(sprintf("'%s' must be created by prior_niw() or prior_vague().",
                 name), call. = FALSE)
  }
  invisible(x)
}
