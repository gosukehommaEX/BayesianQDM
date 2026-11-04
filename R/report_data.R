#' Summary of data in a report
#'
#' @param data Data created by a \code{data_*()} function.
#'
#' @return The description of \code{format()} without the type of endpoint,
#'   for example \code{"n = 20, y = 8"}, with "exposure" written as "total
#'   observed time".
#' @noRd
report_data <- function(data) {
  out <- sub("^[^,]*, ", "", format(data))
  sub("exposure", "total observed time", out, fixed = TRUE)
}
