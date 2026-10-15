#' Arguments of the C++ simulation functions for the conduct of a study
#'
#' @param conduct An object of class \code{"qdm_conduct"}.
#'
#' @return A list with elements \code{acc_time}, \code{acc_cum},
#'   \code{drop_start}, \code{drop_hazard}, \code{cut_time} (infinite for an
#'   analysis at an event) and \code{cut_events} (0 for an analysis at a
#'   time), in the order of the arguments of \code{cpp_tte_oc_data()}.
#' @noRd
conduct_args <- function(conduct) {
  list(acc_time = conduct$accrual_time, acc_cum = conduct$accrual_cum,
       drop_start = conduct$dropout_time,
       drop_hazard = conduct$dropout_hazard,
       cut_time = if (is.null(conduct$cut_time)) Inf else conduct$cut_time,
       cut_events = if (is.null(conduct$cut_events)) 0L else
         conduct$cut_events)
}
