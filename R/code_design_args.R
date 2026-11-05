#' Arguments of a design specification that differ from their defaults
#'
#' @param spec An object of class \code{"qdm_spec"}.
#'
#' @return The arguments of \code{spec}, in the order of the formal arguments
#'   of the \code{design_*()} function of the endpoint, without those equal
#'   to their defaults (for \code{match.arg()}, the first choice; for a
#'   default that is another argument, the value of that argument).
#' @noRd
code_design_args <- function(spec) {
  f <- formals(get(design_function(spec$endpoint), mode = "function"))
  args <- spec$args[intersect(names(f), names(spec$args))]
  keep <- vapply(names(args), function(nm) {
    if (identical(f[[nm]], quote(expr = ))) {
      return(TRUE)
    }
    d <- f[[nm]]
    if (is.symbol(d)) {
      return(!identical(args[[nm]], args[[as.character(d)]]))
    }
    value <- eval(d, baseenv())
    if (is.character(value) && length(value) > 1L) {
      value <- value[1L]
    }
    !identical(args[[nm]], value)
  }, logical(1))
  args[keep]
}
