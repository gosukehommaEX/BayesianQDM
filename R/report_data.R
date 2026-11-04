#' Summary of data in a report
#'
#' @param data Data created by a \code{data_*()} function.
#'
#' @return The description of \code{format()} without the type of endpoint,
#'   for example \code{"n = 20, y = 8"}, with "exposure" written as "total
#'   observed time" and each number formatted separately (see
#'   \code{report_num()}), so that the numbers of a vector do not share the
#'   number of decimals.
#' @noRd
report_data <- function(data) {
  num <- report_num
  switch(
    class(data)[1L],
    qdm_data_binary = sprintf("n = %d, y = %d", data$n, data$y),
    qdm_data_continuous = sprintf("n = %d, mean = %s, sd = %s", data$n,
                                  num(data$mean), num(data$sd)),
    qdm_data_binary2 = sprintf("n = %d, counts (00, 01, 10, 11) = (%s)",
                               data$n, paste(data$counts, collapse = ", ")),
    qdm_data_continuous2 = sprintf("n = %d, mean = (%s), cov = (%s; %s)",
                                   data$n, num(data$mean),
                                   num(data$cov[1L, ]), num(data$cov[2L, ])),
    qdm_data_tte = sprintf("events = %d, total observed time = %s",
                           data$events, num(data$exposure)),
    stop("Unknown data class.", call. = FALSE)
  )
}
