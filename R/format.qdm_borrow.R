#' Format an external data specification
#'
#' @param x An object of class \code{"qdm_borrow"}.
#' @param digits Number of significant digits for numeric values.
#' @param ... Not used.
#'
#' @return A character string describing the external data and the power prior
#'   weight.
#'
#' @examples
#' format(borrow(data_binary(n = 20, y = 8), a0 = 0.5))
#' @export
format.qdm_borrow <- function(x, digits = 4L, ...) {
  sprintf("%s, a0 = %s", format(x$data, digits = digits),
          format(x$a0, digits = digits))
}
