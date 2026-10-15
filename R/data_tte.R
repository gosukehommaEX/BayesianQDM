#' Time-to-event endpoint data
#'
#' Summarizes the data of one group for a time-to-event endpoint by the
#' sufficient statistics of the exponential model: the number of events and
#' the total observed time (the sum over subjects of the time from entry to
#' the event or censoring). Used for observed proof-of-concept data, external
#' data and hypothetical control data.
#'
#' With exponential event times and independent censoring, the likelihood of
#' the hazard \eqn{\lambda} is proportional to
#' \eqn{\lambda^d \exp(-\lambda E)} for \eqn{d} events and total observed
#' time \eqn{E}, so these two numbers are all that the analysis needs.
#'
#' @param events Number of events (a whole number of at least 0).
#' @param exposure Total observed time (positive), in the same unit as the
#'   milestone time of the design.
#'
#' @return An object of class \code{c("qdm_data_tte", "qdm_data")} with
#'   elements \code{events} and \code{exposure}.
#'
#' @family data
#' @examples
#' data_tte(events = 5, exposure = 40)
#' @export
data_tte <- function(events, exposure) {
  check_numeric(events, "events", lower = 0, integer = TRUE)
  check_numeric(exposure, "exposure", lower = 0, lower_open = TRUE)
  structure(list(events = as.integer(round(events)), exposure = exposure),
            class = c("qdm_data_tte", "qdm_data"))
}
