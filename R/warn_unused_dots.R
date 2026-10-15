#' Warn about arguments that a plot method does not use
#'
#' Arguments that a \code{plot()} method does not have are collected by
#' \code{...}. They are ignored with a warning that names them, as
#' \pkg{ggplot2} does for unknown parameters, so that a misspelled argument
#' or one that applies only to another plot is not dropped silently.
#'
#' @param dots The arguments in \code{...}, as a list.
#' @param fun Name of the method, used in the warning.
#'
#' @return \code{NULL}, invisibly.
#' @noRd
warn_unused_dots <- function(dots, fun) {
  if (length(dots) == 0L) {
    return(invisible(NULL))
  }
  nms <- names(dots)
  if (is.null(nms)) {
    nms <- rep("", length(dots))
  }
  nms[nms == ""] <- "an unnamed argument"
  warning(sprintf("Ignoring unknown argument(s) of %s(): %s.", fun,
                  paste(nms, collapse = ", ")), call. = FALSE)
  invisible(NULL)
}
