#' Effects of no effect and of the thresholds of a criterion
#'
#' @param criterion An object of class \code{"qdm_criterion"}.
#' @param effect The effect of a time-to-event design (\code{"milestone"} or
#'   \code{"hazard_ratio"}), or \code{NULL}.
#'
#' @return A list with one numeric vector per endpoint: no effect (0, or 1
#'   for the hazard ratio), the minimum acceptable value and the target value
#'   (posterior criterion) or the null value (predictive criterion), with
#'   repeated values removed.
#' @noRd
threshold_effects <- function(criterion, effect = NULL) {
  none <- if (identical(effect, "hazard_ratio")) 1 else 0
  lapply(seq_len(criterion$n_endpoints), function(j) {
    v <- if (criterion$type == "posterior") {
      c(none, criterion$theta_mav[j], criterion$theta_tv[j])
    } else {
      c(none, criterion$theta_null[j])
    }
    unique(v)
  })
}
