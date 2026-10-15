#' Decision boundary data for one time-to-event endpoint
#'
#' @param design An object of class \code{"qdm_design_tte"}.
#' @param n_t Number of subjects in the treatment group, used for the default
#'   numbers of events.
#' @param events_t Numbers of events of the treatment group (by default
#'   \code{0:n_t}).
#' @param exposure_t Total observed times of the treatment group.
#' @param data_c Observed control data (concurrent control only).
#'
#' @return A list with elements \code{grid} (data frame), \code{p_go},
#'   \code{p_nogo} and \code{axes} (see \code{qdm_boundary()}).
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
  grid <- expand.grid(events_t = sort(unique(round(events_t))),
                      exposure_t = sort(unique(exposure_t)),
                      KEEP.OUT.ATTRS = FALSE)
  post_t <- posterior_tte(design$prior_t, grid$events_t, grid$exposure_t,
                          external = design$external_t)
  post_c <- switch(
    design$control,
    concurrent = posterior_tte(design$prior_c, data_c$events,
                               data_c$exposure, external = design$external_c),
    fixed = NULL,
    posterior_tte(design$prior_c, 0, 0, external = design$external_c,
                  hypothetical = design$hypothetical_c)
  )
  p <- region_prob_tte(design, post_t, post_c)
  axes <- list(x = "exposure_t", y = "events_t",
               xlab = "Total observed time of the treatment group",
               ylab = "Events in the treatment group",
               caption = if (concurrent) {
                 sprintf("Control data: %s", format(data_c, digits = 4L))
               })
  list(grid = grid, p_go = p$go, p_nogo = p$nogo, axes = axes)
}
