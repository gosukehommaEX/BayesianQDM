#' Design specification of a design
#'
#' Converts a design to the design specification that builds it, so that a
#' design can start a pipe of design steps (see \code{\link{qdm_endpoint}}).
#' The arguments of the \code{design_*()} function of the endpoint are taken
#' from the elements of the design, which keep them with the values checked
#' by that function. Arguments that are \code{NULL} are not set.
#'
#' @param design An object of class \code{"qdm_design"}.
#'
#' @return An object of class \code{"qdm_spec"} for which
#'   \code{qdm_build()} gives \code{design}.
#' @noRd
design_spec <- function(design) {
  fun <- get(design_function(design$endpoint), mode = "function")
  elements <- unclass(design)
  args <- elements[intersect(names(formals(fun)), names(elements))]
  args <- args[!vapply(args, is.null, logical(1))]
  structure(list(endpoint = design$endpoint, args = args), class = "qdm_spec")
}
