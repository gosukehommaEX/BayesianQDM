#' Print an external data specification
#'
#' @param x An object of class \code{"qdm_borrow"}.
#' @param ... Passed to \code{\link{format.qdm_borrow}}.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' print(borrow(data_binary(n = 20, y = 8), a0 = 0.5))
#' @export
print.qdm_borrow <- function(x, ...) {
  cat("External data:", format(x, ...), "\n")
  invisible(x)
}
