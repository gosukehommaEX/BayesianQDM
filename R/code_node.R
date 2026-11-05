#' Call in the R code written by show_code()
#'
#' @param fun Name of the function.
#' @param args Named list of the arguments (an empty name for an argument
#'   given by position): character strings of R code or other calls from
#'   \code{code_node()}.
#'
#' @return A list with the elements \code{fun} and \code{args}.
#' @noRd
code_node <- function(fun, args = list()) {
  list(fun = fun, args = args)
}
