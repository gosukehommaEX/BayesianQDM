#' Decision boundary specification for one binary endpoint
#'
#' @param design An object of class \code{"qdm_design_binary"}.
#' @param n_t Numbers of subjects in the treatment group (a single number
#'   with a concurrent control group).
#'
#' @return A list with the elements used by \code{qdm_boundary()}:
#'   \code{cond} (data frame of the conditions), \code{values} (the numbers
#'   of responders of each condition), \code{prob}, \code{discrete},
#'   \code{derived} (function that computes the observed response rate and
#'   difference), \code{axes}, \code{titles} and \code{increasing}.
#' @noRd
boundary_binary <- function(design, n_t) {
  check_numeric(n_t, "n_t", len = NULL, lower = 1, integer = TRUE)
  n_t <- sort(unique(as.integer(round(n_t))))
  crit <- design$criterion
  concurrent <- design$control == "concurrent"
  if (concurrent) {
    if (length(n_t) != 1L) {
      stop(paste("'n_t' must be a single number for a design with a",
                 "concurrent control group."), call. = FALSE)
    }
    n_c <- design$n_c
    cond <- data.frame(n_t = n_t, n_c = n_c, y_c = 0:n_c)
    post_c <- posterior_binary(design$prior_c, n_c, cond$y_c,
                               external = design$external_c)
    control <- cond$y_c / n_c
  } else {
    cond <- data.frame(n_t = n_t)
    post_c <- if (design$control == "fixed") {
      NULL
    } else {
      posterior_binary(design$prior_c, 0, 0, external = design$external_c,
                       hypothetical = design$hypothetical_c)
    }
    control <- rep(control_location(design), nrow(cond))
  }
  prob <- function(idx, y) {
    post_t <- posterior_binary(design$prior_t, cond$n_t[idx], y,
                               external = design$external_t)
    if (concurrent) {
      region_prob_binary(crit, post_t$a, post_t$b, post_c$a[idx],
                         post_c$b[idx])
    } else {
      region_prob_binary(crit, post_t$a, post_t$b, post_c$a, post_c$b,
                         fixed = design$fixed_c)
    }
  }
  derived <- function(idx, y) {
    rate <- y / cond$n_t[idx]
    list(rate = rate, effect = rate - control[idx])
  }
  ref <- format_num(control[1L], 3L)
  effect_title <- switch(
    design$control,
    concurrent = paste("Observed difference in response rates",
                       "(y_t / n_t - y_c / n_c)"),
    fixed = sprintf(paste("Observed difference from the fixed control",
                          "rate (y_t / n_t - %s)"), ref),
    sprintf(paste("Observed difference from the mean of the control",
                  "distribution (y_t / n_t - %s)"), ref)
  )
  axes <- if (concurrent) {
    list(x = "y_c", y = "cut", discrete_x = TRUE, discrete_y = TRUE,
         xlab = sprintf("Responders in the control group (n = %d)", n_c),
         ylab = sprintf("Responders in the treatment group (n = %d)", n_t),
         sec_scale = n_t, caption = NULL)
  } else {
    list(x = "n_t", y = "cut", discrete_x = TRUE, discrete_y = TRUE,
         xlab = "Treatment group size",
         ylab = "Responders in the treatment group", sec_scale = NULL,
         caption = NULL)
  }
  list(cond = cond, values = lapply(cond$n_t, function(n) 0:n), prob = prob,
       discrete = TRUE, derived = derived, axes = axes,
       titles = c(cut = "Responders in the treatment group (y_t)",
                  rate = "Observed response rate (y_t / n_t)",
                  effect = effect_title),
       increasing = crit$direction == "greater")
}
