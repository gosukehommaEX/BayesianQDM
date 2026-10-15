#' Data for two binary endpoints
#'
#' Summarizes the data of one group for two binary endpoints by the counts of
#' the four joint response patterns.
#'
#' The cells are ordered as \code{"00"}, \code{"01"}, \code{"10"} and
#' \code{"11"}, where the first digit is the response on endpoint 1 and the
#' second digit is the response on endpoint 2 (1 = response, 0 = no
#' response).
#'
#' @param counts A numeric vector of four non-negative whole numbers in the
#'   cell order \code{"00"}, \code{"01"}, \code{"10"}, \code{"11"}, with a
#'   total of at least 1.
#'
#' @return An object of class \code{c("qdm_data_binary2", "qdm_data")} with
#'   elements \code{n} (the total count) and \code{counts} (a named integer
#'   vector of length 4).
#'
#' @family data
#' @examples
#' data_binary2(c(4, 2, 3, 6))
#' @export
data_binary2 <- function(counts) {
  check_numeric(counts, "counts", len = 4L, lower = 0, integer = TRUE)
  counts <- stats::setNames(as.integer(round(counts)),
                            c("00", "01", "10", "11"))
  n <- sum(counts)
  if (n < 1L) {
    stop("'counts' must contain at least one subject.", call. = FALSE)
  }
  structure(list(n = n, counts = counts),
            class = c("qdm_data_binary2", "qdm_data"))
}
