#' Conduct of a study with a time-to-event endpoint
#'
#' Specifies how subjects are enrolled, how they drop out and when the
#' analysis takes place. The operating characteristics of a design with a
#' time-to-event endpoint (\code{\link{design_tte}}) are computed by
#' simulating studies conducted in this way, and the posterior predictive
#' probabilities by simulating future studies.
#'
#' The subjects of all groups are enrolled at independent times with a
#' piecewise uniform distribution: the enrollment period is divided by
#' \code{accrual_time}, and the probability of an interval is proportional
#' to its rate in \code{accrual_rate} times its length (only the relative
#' rates matter). The dropout time of a subject, measured from enrollment,
#' has a piecewise exponential distribution with hazard
#' \code{dropout_hazard[j]} from time \code{dropout_time[j]} (the last hazard
#' applies until infinity), independently of the event time. The analysis
#' takes place at the calendar time \code{cut_time} after the start of
#' enrollment, or at the calendar time of the \code{cut_events}-th event
#' among all subjects of the study (if fewer events occur, all subjects are
#' followed until their event or dropout). A subject enrolled at time
#' \eqn{e} with event time \eqn{T} and dropout time \eqn{C} contributes the
#' observed time \eqn{\min(T, C, \textrm{cut} - e)} (0 if enrolled after the
#' analysis) and an event if \eqn{T \le C} and \eqn{e + T \le \textrm{cut}}.
#'
#' @param accrual_time Boundaries of the enrollment intervals: an increasing
#'   numeric vector that starts at 0, of length at least 2.
#' @param accrual_rate Non-negative relative enrollment rates of the
#'   intervals, one fewer than \code{accrual_time}, with a positive sum. The
#'   default is a constant rate (uniform enrollment).
#' @param dropout_time Start times of the pieces of the dropout hazard: an
#'   increasing numeric vector that starts at 0.
#' @param dropout_hazard Non-negative dropout hazards of the pieces, of the
#'   same length as \code{dropout_time}. The default \code{0} means no
#'   dropout.
#' @param cut_time Calendar time of the analysis after the start of
#'   enrollment (positive), or \code{NULL}.
#' @param cut_events Number of events that triggers the analysis (a whole
#'   number of at least 1), or \code{NULL}. Exactly one of \code{cut_time} and
#'   \code{cut_events} must be given.
#'
#' @return An object of class \code{"qdm_conduct"} with the arguments and the
#'   element \code{accrual_cum}, the cumulative probabilities of enrollment
#'   at the ends of the intervals.
#'
#' @seealso \code{\link{design_tte}}
#' @examples
#' # Enrollment over 15 months with increasing rates, dropout hazard 0.00226
#' # per month, analysis 15.44 months after the start of enrollment
#' conduct_tte(accrual_time = c(0, 5, 10, 15), accrual_rate = c(0.4, 2, 4.5),
#'             dropout_hazard = 0.00226, cut_time = 15.44)
#'
#' # Uniform enrollment over 12 months, analysis at the 20th event
#' conduct_tte(accrual_time = c(0, 12), cut_events = 20)
#' @export
conduct_tte <- function(accrual_time,
                        accrual_rate = rep(1, length(accrual_time) - 1L),
                        dropout_time = 0, dropout_hazard = 0,
                        cut_time = NULL, cut_events = NULL) {
  check_numeric(accrual_time, "accrual_time", len = NULL, lower = 0)
  if (length(accrual_time) < 2L || accrual_time[1L] != 0 ||
      any(diff(accrual_time) <= 0)) {
    stop(paste("'accrual_time' must be an increasing vector of at least two",
               "values that starts at 0."), call. = FALSE)
  }
  check_numeric(accrual_rate, "accrual_rate", len = length(accrual_time) - 1L,
                lower = 0)
  if (sum(accrual_rate) <= 0) {
    stop("'accrual_rate' must have a positive sum.", call. = FALSE)
  }
  check_numeric(dropout_time, "dropout_time", len = NULL, lower = 0)
  if (dropout_time[1L] != 0 || any(diff(dropout_time) <= 0)) {
    stop("'dropout_time' must be an increasing vector that starts at 0.",
         call. = FALSE)
  }
  check_numeric(dropout_hazard, "dropout_hazard", len = length(dropout_time),
                lower = 0)
  if (is.null(cut_time) == is.null(cut_events)) {
    stop("Specify exactly one of 'cut_time' and 'cut_events'.",
         call. = FALSE)
  }
  if (!is.null(cut_time)) {
    check_numeric(cut_time, "cut_time", lower = 0, lower_open = TRUE)
  } else {
    check_numeric(cut_events, "cut_events", lower = 1, integer = TRUE)
    cut_events <- as.integer(round(cut_events))
  }
  w <- accrual_rate * diff(accrual_time)
  cum <- cumsum(w) / sum(w)
  cum[length(cum)] <- 1
  structure(list(accrual_time = as.numeric(accrual_time),
                 accrual_rate = as.numeric(accrual_rate), accrual_cum = cum,
                 dropout_time = as.numeric(dropout_time),
                 dropout_hazard = as.numeric(dropout_hazard),
                 cut_time = cut_time, cut_events = cut_events),
            class = "qdm_conduct")
}
