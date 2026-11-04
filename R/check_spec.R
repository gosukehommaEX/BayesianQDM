#' Check that an argument is a design specification
#'
#' @param x Object to check.
#'
#' @return \code{x}, invisibly.
#' @noRd
check_spec <- function(x) {
  if (!inherits(x, "qdm_spec")) {
    stop("'x' must be a design specification started by qdm_endpoint().",
         call. = FALSE)
  }
  invisible(x)
}
