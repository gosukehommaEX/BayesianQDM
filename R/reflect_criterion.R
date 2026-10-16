#' Criterion on the scale of the reflected effects
#'
#' For an endpoint on which smaller effects are favorable
#' (\code{direction = "less"}), the Go and NoGo regions of the effect
#' \eqn{\theta} are those of the reflected effect \eqn{-\theta} with the
#' thresholds \eqn{-\theta_{TV}}, \eqn{-\theta_{MAV}} and
#' \eqn{-\theta_{NULL}} and larger effects favorable: \eqn{\theta < \theta_{TV}}
#' if and only if \eqn{-\theta > -\theta_{TV}}, and \eqn{\theta \ge
#' \theta_{MAV}} if and only if \eqn{-\theta \le -\theta_{MAV}}. The functions
#' that compute the probabilities of the regions reflect the distributions of
#' the effects of these endpoints and use the criterion returned here.
#'
#' @param criterion An object of class \code{"qdm_criterion"}.
#'
#' @return A list with elements \code{criterion} (the criterion with the
#'   thresholds of the reflected endpoints negated and \code{direction} set
#'   to \code{"greater"}) and \code{flip} (logical vector, one element per
#'   endpoint, \code{TRUE} for the reflected endpoints).
#' @noRd
reflect_criterion <- function(criterion) {
  flip <- criterion$direction == "less"
  if (any(flip)) {
    for (nm in intersect(c("theta_tv", "theta_mav", "theta_null"),
                         names(criterion))) {
      criterion[[nm]][flip] <- -criterion[[nm]][flip]
    }
    criterion$direction <- rep("greater", criterion$n_endpoints)
  }
  list(criterion = criterion, flip = flip)
}
