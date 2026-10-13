#' Format a data specification
#'
#' @param x An object of class \code{"qdm_data"}.
#' @param digits Number of significant digits for numeric values.
#' @param ... Not used.
#'
#' @return A character string describing the data.
#'
#' @examples
#' format(data_binary(n = 12, y = 7))
#' @export
format.qdm_data <- function(x, digits = 4L, ...) {
  switch(
    class(x)[1L],
    qdm_data_binary = sprintf("binary, n = %d, y = %d", x$n, x$y),
    qdm_data_continuous = sprintf("continuous, n = %d, mean = %s, sd = %s",
                                  x$n, format(x$mean, digits = digits),
                                  format(x$sd, digits = digits)),
    qdm_data_binary2 = sprintf(
      "two binary endpoints, n = %d, counts (00, 01, 10, 11) = (%s)",
      x$n, paste(x$counts, collapse = ", ")
    ),
    qdm_data_continuous2 = sprintf(
      "two continuous endpoints, n = %d, mean = (%s), cov = (%s; %s)",
      x$n, format_num(x$mean, digits),
      format_num(x$cov[1L, ], digits),
      format_num(x$cov[2L, ], digits)
    ),
    stop("Unknown data class.", call. = FALSE)
  )
}
