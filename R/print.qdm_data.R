#' Print a data specification
#'
#' @param x An object of class \code{"qdm_data"}.
#' @param ... Passed to \code{\link{format.qdm_data}}.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' print(data_continuous(n = 12, mean = 3.2, sd = 2))
#' @export
print.qdm_data <- function(x, ...) {
  cat("Data:", format(x, ...), "\n")
  invisible(x)
}
