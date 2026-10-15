#' Check that an argument is a design
#'
#' @param design Object to check.
#'
#' @return \code{design}, invisibly.
#' @noRd
check_design <- function(design) {
  if (!inherits(design, "qdm_design")) {
    stop(paste("'design' must be created by a design_*() function, such as",
               "design_binary()."), call. = FALSE)
  }
  invisible(design)
}
