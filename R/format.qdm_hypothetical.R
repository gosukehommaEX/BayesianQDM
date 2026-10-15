#' Format a hypothetical control specification
#'
#' @param x An object of class \code{"qdm_hypothetical"}.
#' @param digits Number of significant digits for numeric values.
#' @param ... Not used.
#'
#' @return A character string describing the hypothetical control.
#'
#' @examples
#' format(hypothetical(mean = 1, r = 1))
#' @export
format.qdm_hypothetical <- function(x, digits = 4L, ...) {
  if (x$type == "data") {
    sprintf("data (%s)", format(x$data, digits = digits))
  } else {
    sprintf("mean = (%s), r = %s",
            format_num(x$mean, digits),
            format(x$r, digits = digits))
  }
}
