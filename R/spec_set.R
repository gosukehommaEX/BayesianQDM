#' Set arguments of a design specification
#'
#' Replaces arguments of the \code{design_*()} function stored in a design
#' specification. An argument given as \code{NULL} is removed, so that the
#' \code{design_*()} function uses its default.
#'
#' @param x A design specification (class \code{"qdm_spec"}).
#' @param values Named list of arguments of the \code{design_*()} function of
#'   the endpoint.
#' @param step Name of the function that sets the arguments, used in the
#'   error message.
#'
#' @return \code{x} with the arguments replaced.
#' @noRd
spec_set <- function(x, values, step) {
  check_spec(x)
  fun <- design_function(x$endpoint)
  used <- !vapply(values, is.null, logical(1))
  unknown <- setdiff(names(values)[used],
                     names(formals(get(fun, mode = "function"))))
  if (length(unknown) > 0L) {
    stop(sprintf("%s() cannot set %s, which %s() does not have.", step,
                 paste0("'", unknown, "'", collapse = ", "), fun),
         call. = FALSE)
  }
  keep <- setdiff(names(x$args), names(values))
  x$args <- c(x$args[keep], values[used])
  x
}
