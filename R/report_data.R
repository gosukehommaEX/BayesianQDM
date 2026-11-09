#' Summary of data in a report
#'
#' @param data Data created by a \code{data_*()} function.
#'
#' @return A noun phrase that describes the data in a sentence of a report,
#'   for example \code{"9 responders among 12 subjects"}. Each number is
#'   formatted separately (see \code{report_num()}), so that the numbers of a
#'   vector do not share the number of decimals, and the noun after a count of
#'   one is singular.
#' @noRd
report_data <- function(data) {
  num <- report_num
  # A count of a noun, singular for one
  count <- function(k, noun) {
    sprintf("%d %s%s", k, noun, if (k == 1) "" else "s")
  }
  switch(
    class(data)[1L],
    qdm_data_binary = sprintf("%s among %s", count(data$y, "responder"),
                              count(data$n, "subject")),
    qdm_data_continuous = sprintf(
      "the outcomes of %s, with mean %s and standard deviation %s",
      count(data$n, "subject"), num(data$mean), num(data$sd)
    ),
    qdm_data_binary2 = sprintf(
      "%s with the counts (%s) of the responses 00, 01, 10 and 11",
      count(data$n, "subject"), paste(data$counts, collapse = ", ")
    ),
    qdm_data_continuous2 = sprintf(
      paste("the outcomes of %s, with means (%s) and covariance matrix",
            "$\\begin{pmatrix} %s & %s \\\\ %s & %s \\end{pmatrix}$"),
      count(data$n, "subject"), num(data$mean), num(data$cov[1L, 1L]),
      num(data$cov[1L, 2L]), num(data$cov[2L, 1L]), num(data$cov[2L, 2L])
    ),
    qdm_data_tte = sprintf("%s in a total observed time of %s",
                           count(data$events, "event"), num(data$exposure)),
    stop("Unknown data class.", call. = FALSE)
  )
}
