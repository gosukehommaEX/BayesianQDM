#' Decision boundary specification for one time-to-event endpoint
#'
#' The conditions are the total observed times of the treatment group and
#' the cut variable is its number of events.
#'
#' @param design An object of class \code{"qdm_design_tte"}.
#' @param n_t Number of subjects in the treatment group, used for the default
#'   numbers of events.
#' @param events_t Numbers of events of the treatment group (by default
#'   \code{0:n_t}).
#' @param exposure_t Total observed times of the treatment group.
#' @param data_c Observed control data (concurrent control only).
#'
#' @return A list with the elements used by \code{qdm_boundary()} (see
#'   \code{boundary_binary()}).
#' @noRd
boundary_tte <- function(design, n_t, events_t, exposure_t, data_c) {
  if (is.null(exposure_t)) {
    stop("'exposure_t' must be given for a time-to-event endpoint.",
         call. = FALSE)
  }
  check_numeric(n_t, "n_t", lower = 1, integer = TRUE)
  if (is.null(events_t)) {
    events_t <- 0:round(n_t)
  }
  check_numeric(events_t, "events_t", len = NULL, lower = 0, integer = TRUE)
  check_numeric(exposure_t, "exposure_t", len = NULL, lower = 0,
                lower_open = TRUE)
  concurrent <- design$control == "concurrent"
  if (concurrent && !inherits(data_c, "qdm_data_tte")) {
    stop("'data_c' must be created by data_tte().", call. = FALSE)
  }
  crit <- design$criterion
  hr <- identical(design$effect, "hazard_ratio")
  t0 <- design$milestone
  cond <- data.frame(exposure_t = sort(unique(exposure_t)))
  events <- sort(unique(round(events_t)))
  post_c <- switch(
    design$control,
    concurrent = posterior_tte(design$prior_c, data_c$events,
                               data_c$exposure, external = design$external_c),
    fixed = NULL,
    posterior_tte(design$prior_c, 0, 0, external = design$external_c,
                  hypothetical = design$hypothetical_c)
  )
  prob <- function(idx, d) {
    post_t <- posterior_tte(design$prior_t, d, cond$exposure_t[idx],
                            external = design$external_t)
    region_prob_tte(design, post_t, post_c)
  }
  # Control value of the observed effects: the observed hazard or milestone
  # survival of a concurrent control, otherwise the fixed value or the mean
  # of the control distribution (control_location())
  control <- if (concurrent) {
    haz_c <- data_c$events / data_c$exposure
    if (hr) haz_c else exp(-t0 * haz_c)
  } else {
    control_location(design)
  }
  ref <- format_num(control, 4L)
  derived <- function(idx, d) {
    haz <- d / cond$exposure_t[idx]
    if (hr) {
      # Without control events the observed hazard ratio is not defined
      list(effect = if (control > 0) haz / control else NA_real_ * haz)
    } else {
      surv <- exp(-t0 * haz)
      list(surv = surv, effect = surv - control)
    }
  }
  what <- switch(design$control,
                 concurrent = "observed in the control group",
                 fixed = "fixed control value",
                 "mean of the control distribution")
  titles <- if (hr) {
    c(cut = "Events in the treatment group (d_t)",
      effect = sprintf(paste("Observed hazard ratio: d_t / E_t divided by",
                             "the control hazard %s (%s)"), ref, what))
  } else {
    c(cut = "Events in the treatment group (d_t)",
      surv = sprintf("Observed milestone survival exp(-%s d_t / E_t)",
                     format_num(t0, 4L)),
      effect = sprintf(paste("Observed difference in milestone survival",
                             "(control %s, %s)"), ref, what))
  }
  list(cond = cond, values = rep(list(events), nrow(cond)), prob = prob,
       discrete = TRUE, derived = derived,
       axes = list(x = "exposure_t", y = "cut", discrete_x = FALSE,
                   discrete_y = TRUE,
                   xlab = "Total observed time of the treatment group",
                   ylab = "Events in the treatment group", sec_scale = NULL,
                   caption = if (concurrent) {
                     sprintf("Control data: %s", format(data_c, digits = 4L))
                   }),
       # More events give a larger hazard: a smaller milestone survival but
       # a larger hazard ratio
       titles = titles, increasing = xor(hr, crit$direction == "less"))
}
