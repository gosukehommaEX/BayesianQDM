#' Validate a logical flag
#'
#' Checks that an argument is a single non-missing logical value.
#'
#' @param x Object to check.
#' @param name Name of the argument, used in error messages.
#'
#' @return \code{x}, invisibly.
#' @noRd
check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("'%s' must be TRUE or FALSE.", name), call. = FALSE)
  }
  invisible(x)
}
