#' Posterior distributions and region probabilities for observed
#' time-to-event data
#'
#' @param design An object of class \code{"qdm_design_tte"}.
#' @param data_t,data_c Observed data of the treatment and control groups;
#'   \code{data_c} is \code{NULL} for a design without a concurrent control
#'   group.
#'
#' @return A list with elements \code{prob_go}, \code{prob_nogo},
#'   \code{posterior_t} and \code{posterior_c} (named vectors of the shape
#'   and rate of the Gamma distribution of the hazard; for a fixed control,
#'   \code{posterior_c} is the fixed value named \code{value}).
#' @noRd
decide_tte <- function(design, data_t, data_c) {
  if (!inherits(data_t, "qdm_data_tte")) {
    stop("'data_t' must be created by data_tte().", call. = FALSE)
  }
  if (design$control == "concurrent") {
    if (!inherits(data_c, "qdm_data_tte")) {
      stop(paste("'data_c' must be created by data_tte() for a design",
                 "with a concurrent control group."), call. = FALSE)
    }
    post_c <- posterior_tte(design$prior_c, data_c$events, data_c$exposure,
                            external = design$external_c)
  } else {
    if (!is.null(data_c)) {
      stop(paste("'data_c' must be NULL for a design without a concurrent",
                 "control group."), call. = FALSE)
    }
    post_c <- if (design$control == "fixed") {
      NULL
    } else {
      posterior_tte(design$prior_c, 0, 0, external = design$external_c,
                    hypothetical = design$hypothetical_c)
    }
  }
  post_t <- posterior_tte(design$prior_t, data_t$events, data_t$exposure,
                          external = design$external_t)
  p <- region_prob_tte(design, post_t, post_c)
  list(prob_go = p$go, prob_nogo = p$nogo,
       posterior_t = c(shape = post_t$shape, rate = post_t$rate),
       posterior_c = if (is.null(post_c)) {
         c(value = design$fixed_c)
       } else {
         c(shape = post_c$shape, rate = post_c$rate)
       })
}
