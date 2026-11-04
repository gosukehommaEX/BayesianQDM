#' Check that an argument is a design
#'
#' A design specification started by \code{qdm_endpoint()} is built by
#' \code{qdm_build()}.
#'
#' @param design Object to check.
#'
#' @return \code{design}, invisibly, or the design built from a design
#'   specification.
#' @noRd
check_design <- function(design) {
  if (inherits(design, "qdm_spec")) {
    return(qdm_build(design))
  }
  if (!inherits(design, "qdm_design")) {
    stop(paste("'design' must be created by a design_*() function, such as",
               "design_binary(), or by qdm_endpoint()."), call. = FALSE)
  }
  invisible(design)
}
