#' Validate a color argument
#'
#' Checks that an argument is a single color that R recognizes: a color name
#' (see \code{grDevices::colors()}) or a hexadecimal code such as
#' \code{"#658D1B"}.
#'
#' @param x Object to check.
#' @param name Name of the argument, used in error messages.
#'
#' @return \code{x}, invisibly.
#' @noRd
check_color <- function(x, name) {
  ok <- is.character(x) && length(x) == 1L && !is.na(x) &&
    !inherits(tryCatch(grDevices::col2rgb(x), error = function(e) e),
              "error")
  if (!ok) {
    stop(sprintf(paste("'%s' must be a single color name or hexadecimal",
                       "color code."), name), call. = FALSE)
  }
  invisible(x)
}
