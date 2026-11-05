#' Check that an argument is a design specification or a design
#'
#' A design is converted to the specification that builds it
#' (\code{design_spec()}), so that the design steps accept a design.
#'
#' @param x Object to check.
#'
#' @return \code{x}, invisibly, or the specification of the design
#'   \code{x}.
#' @noRd
check_spec <- function(x) {
  if (inherits(x, "qdm_design")) {
    return(design_spec(x))
  }
  if (!inherits(x, "qdm_spec")) {
    stop(paste("'x' must be a design specification started by",
               "qdm_endpoint(), or a design created by a design_*()",
               "function."), call. = FALSE)
  }
  invisible(x)
}
