#' Scenarios for the operating characteristics
#'
#' Creates the data frame of scenarios for the true parameters that
#' \code{\link{qdm_oc}}, \code{\link{qdm_compare}}, \code{\link{qdm_calibrate}}
#' and \code{\link{qdm_find_n}} use, with the columns that the design needs
#' (see \code{\link{qdm_oc}}). Missing, unknown and repeated columns are
#' errors.
#'
#' The values of the columns are given in \code{...}, for example
#' \code{pi_t = seq(0.2, 0.8, by = 0.1), pi_c = 0.4} for one binary
#' endpoint. A name without the letter of the group gives the same value to
#' both groups: \code{sd} for a continuous endpoint (\code{sd_t} and
#' \code{sd_c}), \code{rho} for two binary endpoints (\code{rho_t} and
#' \code{rho_c}), and \code{sd1}, \code{sd2} and \code{rho} for two
#' continuous endpoints. The columns of the control group are not used for a
#' design without a concurrent control group. The correlations of two
#' endpoints have no default value, because the operating characteristics
#' depend on them (for two binary endpoints, the feasible correlations also
#' depend on the response rates; see \code{\link{cellprob_binary2}}).
#'
#' Alternatively, the scenarios can be given by the treatment effect
#' \code{theta}, on the scale of the decision criterion and of the
#' horizontal axis of the plots of \code{\link{qdm_oc}}. The parameter of the
#' treatment group that sets the effect (\code{pi_t}, \code{mu_t} or
#' \code{median_t}; \code{pi_t1} and \code{pi_t2}; or \code{mu_t1} and
#' \code{mu_t2}) is then computed, and is not given in \code{...}. With a
#' concurrent control group, the effect is measured from the parameter of
#' the control group given in \code{...}; without one, from the fixed
#' control value, or from the mean (binary and time-to-event endpoints) or
#' the location (continuous endpoints) of the distribution of the control
#' group. For a time-to-event endpoint evaluated by the milestone survival,
#' the median survival time of the treatment group is
#' \eqn{M_t = -\log(2) t_0 / \log(S_c(t_0) + \theta)}, where \eqn{S_c(t_0)}
#' is the milestone survival probability of the control group; for the
#' hazard ratio \eqn{\theta = \lambda_t / \lambda_c}, it is
#' \eqn{M_t = M_c / \theta}, or \eqn{\log(2) / (\theta \lambda_c)} without a
#' concurrent control group, where \eqn{\lambda_c} is the fixed control
#' hazard or the mean of the distribution of the control hazard.
#'
#' \code{theta = "thresholds"} gives the effects that a study protocol
#' usually reports: no effect (0, or 1 for a hazard ratio), the minimum
#' acceptable value and the target value of the criterion (the null value of
#' a predictive criterion), for each endpoint (repeated values are given
#' once).
#'
#' @param design A design created by \code{\link{design_binary}},
#'   \code{\link{design_continuous}}, \code{\link{design_binary2}},
#'   \code{\link{design_continuous2}} or \code{\link{design_tte}}, or a design
#'   specification started by \code{\link{qdm_endpoint}}, which is built by
#'   \code{\link{qdm_build}}.
#' @param ... Named numeric vectors: the values of the columns of the
#'   scenarios.
#' @param theta \code{NULL} (the default) to give all the columns in
#'   \code{...}; the treatment effects, a numeric vector for one endpoint or
#'   a list of two numeric vectors (the effects on endpoints 1 and 2) for two
#'   endpoints; or \code{"thresholds"}.
#' @param grid \code{TRUE} (the default) for all combinations of the values,
#'   as \code{expand.grid()} (the effects, then the values in the order
#'   given, the first varying fastest), or \code{FALSE} to combine the values
#'   row by row, each value having length 1 or the length of the longest
#'   value.
#'
#' @return A data frame with one row per scenario: the column \code{theta}
#'   (or \code{theta1} and \code{theta2}) when \code{theta} is given,
#'   followed by the columns of the scenarios in the order given in
#'   \code{\link{qdm_oc}}.
#'
#' @family operations
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' qdm_scenarios(des, pi_t = seq(0.2, 0.8, by = 0.2), pi_c = 0.4)
#'
#' # Scenarios given by the effect
#' qdm_scenarios(des, theta = c(0, 0.1, 0.2, 0.3), pi_c = 0.4)
#'
#' # No effect, the minimum acceptable value and the target value
#' scen <- qdm_scenarios(des, theta = "thresholds", pi_c = 0.4)
#' qdm_oc(des, scen)
#'
#' # Two continuous endpoints with the same standard deviations and
#' # correlation in both groups
#' des2 <- design_continuous2(
#'   n_t = 12, n_c = 12, prior_t = prior_vague(), prior_c = prior_vague(),
#'   criterion = criterion_posterior(theta_tv = c(10, 2), theta_mav = c(5, 1),
#'                                   gamma_go = 0.5, gamma_nogo = 0.3)
#' )
#' qdm_scenarios(des2, theta = list(c(0, 5, 10), c(0, 1, 2)), mu_c1 = 0,
#'               mu_c2 = 0, sd1 = 20, sd2 = 4, rho = 0.3)
#' @export
qdm_scenarios <- function(design, ..., theta = NULL, grid = TRUE) {
  design <- check_design(design)
  check_flag(grid, "grid")
  values <- list(...)
  nms <- names(values)
  if (length(values) > 0L &&
      (is.null(nms) || anyNA(nms) || any(nms == ""))) {
    stop("The values of the columns in '...' must be named.", call. = FALSE)
  }
  if (anyDuplicated(nms) > 0L) {
    stop("Each name in '...' must be given only once.", call. = FALSE)
  }
  cols <- scenario_columns(design)$names
  shared <- scenario_shared(design$endpoint, cols)
  loc <- scenario_location(design)
  effect <- if (is.null(theta)) NULL else scenario_theta(design, theta)
  needed <- if (is.null(effect)) cols else setdiff(cols, loc)
  # Columns set by each name ("sd" sets "sd_t" and "sd_c")
  covered <- unlist(lapply(nms, function(nm) {
    if (nm %in% names(shared)) shared[[nm]] else nm
  }), use.names = FALSE)
  if (!is.null(effect) && any(loc %in% covered)) {
    stop(sprintf("%s is computed from 'theta' and must not be given in '...'.",
                 paste0("'", intersect(loc, covered), "'", collapse = " and ")),
         call. = FALSE)
  }
  unknown <- setdiff(nms, c(needed, names(shared)))
  if (length(unknown) > 0L) {
    stop(sprintf(paste("The scenarios of this design have no column %s; they",
                       "have the column(s) %s."),
                 paste0("'", unknown, "'", collapse = ", "),
                 paste(needed, collapse = ", ")), call. = FALSE)
  }
  repeated <- unique(covered[duplicated(covered)])
  if (length(repeated) > 0L) {
    stop(sprintf("The column(s) %s are given more than once in '...'.",
                 paste0("'", repeated, "'", collapse = ", ")), call. = FALSE)
  }
  absent <- setdiff(needed, covered)
  if (length(absent) > 0L) {
    stop(sprintf("The column(s) %s must be given in '...'.",
                 paste0("'", absent, "'", collapse = ", ")), call. = FALSE)
  }
  for (nm in nms) {
    check_numeric(values[[nm]], nm, len = NULL)
  }
  if (!is.null(effect)) {
    names(effect) <- if (length(effect) == 1L) "theta" else
      c("theta1", "theta2")
  }
  vals <- c(effect, values)
  if (grid) {
    out <- expand.grid(vals, KEEP.OUT.ATTRS = FALSE)
  } else {
    lens <- lengths(vals)
    n <- max(lens)
    if (any(lens != 1L & lens != n)) {
      stop(paste("With 'grid = FALSE', each value must have length 1 or the",
                 "length of the longest value."), call. = FALSE)
    }
    out <- data.frame(lapply(vals, rep_len, length.out = n))
  }
  for (nm in intersect(names(shared), nms)) {
    for (target in shared[[nm]]) {
      out[[target]] <- out[[nm]]
    }
    out[[nm]] <- NULL
  }
  if (!is.null(effect)) {
    out <- scenario_treatment(design, out)
  }
  out <- out[c(names(effect), cols)]
  rownames(out) <- NULL
  check_scenarios(design, out)
  out
}
