#' Parameters of the treatment group for given effects
#'
#' Computes the columns of the treatment group that set the effect (see
#' \code{scenario_location()}) from the effects in the columns \code{theta}
#' (or \code{theta1} and \code{theta2}) and the control group: the columns of
#' the control group with a concurrent control group, or else the fixed
#' control value or the location of the distribution of the control group
#' (\code{control_location()}). This is the inverse of the horizontal axis of
#' the plots of operating characteristics.
#'
#' @param design An object of class \code{"qdm_design"}.
#' @param x A data frame of scenarios with the effects.
#'
#' @return \code{x} with the columns of the treatment group added.
#' @noRd
scenario_treatment <- function(design, x) {
  conc <- design$control == "concurrent"
  center <- if (conc) NULL else control_location(design)
  if (design$criterion$n_endpoints == 2L) {
    par <- if (design$endpoint == "binary2") "pi" else "mu"
    for (j in 1:2) {
      base <- if (conc) x[[sprintf("%s_c%d", par, j)]] else center[j]
      x[[sprintf("%s_t%d", par, j)]] <- base + x[[sprintf("theta%d", j)]]
    }
    return(x)
  }
  theta <- x$theta
  if (design$endpoint == "binary") {
    x$pi_t <- (if (conc) x$pi_c else center) + theta
  } else if (design$endpoint == "continuous") {
    x$mu_t <- (if (conc) x$mu_c else center) + theta
  } else if (identical(design$effect, "hazard_ratio")) {
    if (any(theta <= 0)) {
      stop("'theta' must be positive for the hazard ratio.", call. = FALSE)
    }
    # theta = lambda_t / lambda_c, with the hazard log(2) / M
    x$median_t <- if (conc) x$median_c / theta else log(2) / (theta * center)
  } else {
    t0 <- design$milestone
    s_c <- if (conc) exp(-log(2) * t0 / x$median_c) else center
    s_t <- s_c + theta
    if (any(s_t <= 0 | s_t >= 1)) {
      stop(paste("'theta' gives milestone survival probabilities of the",
                 "treatment group outside (0, 1)."), call. = FALSE)
    }
    x$median_t <- -log(2) * t0 / log(s_t)
  }
  x
}
