#' Calibration of the decision thresholds
#'
#' Chooses the thresholds \eqn{\gamma_{go}} and \eqn{\gamma_{nogo}} of the
#' criterion of a design so that the probabilities of a false Go decision and
#' a false NoGo decision are controlled.
#'
#' The threshold \eqn{\gamma_{go}} is the smallest value in \code{grid} for
#' which the probability that the Go criterion is met under
#' \code{scenario_go} (typically a scenario with no or a small treatment
#' effect) is at most \code{target_go}. The threshold \eqn{\gamma_{nogo}} is
#' the smallest value in \code{grid} for which the probability that the NoGo
#' criterion is met under \code{scenario_nogo} (typically a scenario with the
#' expected treatment effect) is at most \code{target_nogo}. Both
#' probabilities decrease as the threshold increases, so the chosen thresholds
#' are the least strict ones that meet the targets.
#'
#' The decision probabilities of the design with the chosen thresholds are
#' computed under both scenarios. If \eqn{\gamma_{go} + \gamma_{nogo} \le 1}
#' (and, for two endpoints with the marginal rule \code{"or"}, for any
#' thresholds), Miss decisions (both criteria met) can occur, and a warning is
#' given when they have positive probability.
#'
#' @param design A design created by \code{\link{design_binary}},
#'   \code{\link{design_continuous}}, \code{\link{design_binary2}},
#'   \code{\link{design_continuous2}} or \code{\link{design_tte}}. Its
#'   thresholds are replaced by the calibrated ones. For two endpoints with a
#'   marginal rule (\code{rule = "and"} or \code{"or"}), the same threshold is
#'   used for both endpoints.
#' @param scenario_go,scenario_nogo Data frames with one row each, giving the
#'   scenarios used to calibrate \eqn{\gamma_{go}} and \eqn{\gamma_{nogo}}
#'   (see \code{\link{qdm_oc}} for the columns). For two binary endpoints, the
#'   correlations must be feasible.
#' @param target_go,target_nogo Upper bounds for the probabilities that the
#'   Go criterion is met under \code{scenario_go} and that the NoGo criterion
#'   is met under \code{scenario_nogo}, strictly between 0 and 1.
#' @param grid Candidate thresholds, strictly between 0 and 1.
#' @param nsim,seed As in \code{\link{qdm_oc}}. The probabilities for all
#'   values of \code{grid} are computed from the same outcomes (simulated
#'   studies or enumerated outcomes) and the same Monte Carlo draws.
#'
#' @return An object of class \code{"qdm_calibration"}: a list with elements
#'   \describe{
#'     \item{gamma_go, gamma_nogo}{The chosen thresholds (\code{NA} if no
#'       value in \code{grid} meets the target).}
#'     \item{prob_go, prob_nogo}{The probabilities that the Go criterion is
#'       met under \code{scenario_go} and that the NoGo criterion is met under
#'       \code{scenario_nogo}, at the chosen thresholds.}
#'     \item{target_go, target_nogo}{The targets.}
#'     \item{curve}{A data frame with the columns \code{gamma} (the grid),
#'       \code{prob_go} and \code{prob_nogo} (the two probabilities at each
#'       threshold).}
#'     \item{oc}{A data frame of the decision probabilities with the chosen
#'       thresholds under the two scenarios (\code{NULL} if a threshold was
#'       not found).}
#'     \item{design}{The design with the chosen thresholds (\code{NULL} if a
#'       threshold was not found).}
#'     \item{scenario_go, scenario_nogo, method}{The scenarios and the method
#'       of computation.}
#'   }
#'
#' @family operations
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' cal <- qdm_calibrate(des,
#'                      scenario_go = data.frame(pi_t = 0.4, pi_c = 0.4),
#'                      scenario_nogo = data.frame(pi_t = 0.7, pi_c = 0.4),
#'                      target_go = 0.05, target_nogo = 0.2)
#' cal
#' plot(cal)
#' cal$design
#' @export
qdm_calibrate <- function(design, scenario_go, scenario_nogo, target_go,
                          target_nogo, grid = seq(0.01, 0.99, by = 0.01),
                          nsim = 10000L, seed = NULL) {
  check_design(design)
  scenario_go <- check_scenarios(design, scenario_go, "scenario_go",
                                 one_row = TRUE, feasible = TRUE)
  scenario_nogo <- check_scenarios(design, scenario_nogo, "scenario_nogo",
                                   one_row = TRUE, feasible = TRUE)
  check_numeric(target_go, "target_go", lower = 0, upper = 1,
                lower_open = TRUE, upper_open = TRUE)
  check_numeric(target_nogo, "target_nogo", lower = 0, upper = 1,
                lower_open = TRUE, upper_open = TRUE)
  check_numeric(grid, "grid", len = NULL, lower = 0, upper = 1,
                lower_open = TRUE, upper_open = TRUE)
  grid <- sort(unique(grid))
  check_numeric(nsim, "nsim", lower = 1, upper = .Machine$integer.max,
                integer = TRUE)
  if (!is.null(seed)) {
    check_numeric(seed, "seed", lower = -.Machine$integer.max,
                  upper = .Machine$integer.max, integer = TRUE)
  }
  cols <- scenario_columns(design)$names
  scenarios <- rbind(scenario_go[, cols, drop = FALSE],
                     scenario_nogo[, cols, drop = FALSE])
  sample <- with_seed(seed, oc_sample(design, scenarios, nsim))
  rule <- design$criterion$rule
  curve <- data.frame(
    gamma = grid,
    prob_go = met_curve(sample_column(sample$p_go, 1L, sample$per_scenario),
                        sample_weights(sample, 1L), grid, rule),
    prob_nogo = met_curve(sample_column(sample$p_nogo, 2L,
                                        sample$per_scenario),
                          sample_weights(sample, 2L), grid, rule)
  )
  gamma_go <- select_gamma(grid, curve$prob_go, target_go)
  gamma_nogo <- select_gamma(grid, curve$prob_nogo, target_nogo)
  if (is.na(gamma_go)) {
    warning(sprintf(paste("No value in 'grid' gives a probability of meeting",
                          "the Go criterion of at most %s under",
                          "'scenario_go'."), format(target_go)),
            call. = FALSE)
  }
  if (is.na(gamma_nogo)) {
    warning(sprintf(paste("No value in 'grid' gives a probability of meeting",
                          "the NoGo criterion of at most %s under",
                          "'scenario_nogo'."), format(target_nogo)),
            call. = FALSE)
  }
  oc <- NULL
  new_design <- NULL
  if (!is.na(gamma_go) && !is.na(gamma_nogo)) {
    new_design <- design
    k <- length(design$criterion$gamma_go)
    new_design$criterion$gamma_go <- rep(gamma_go, k)
    new_design$criterion$gamma_nogo <- rep(gamma_nogo, k)
    probs <- oc_from_sample(sample, new_design$criterion$gamma_go,
                            new_design$criterion$gamma_nogo, rule)
    oc <- cbind(data.frame(scenario = c("scenario_go", "scenario_nogo")),
                scenarios, as.data.frame(probs))
    rownames(oc) <- NULL
    warn_miss(oc$Miss, rule)
  }
  structure(list(gamma_go = gamma_go, gamma_nogo = gamma_nogo,
                 prob_go = curve$prob_go[match(gamma_go, grid)],
                 prob_nogo = curve$prob_nogo[match(gamma_nogo, grid)],
                 target_go = target_go, target_nogo = target_nogo,
                 curve = curve, oc = oc, design = new_design,
                 scenario_go = scenario_go, scenario_nogo = scenario_nogo,
                 method = sample$method),
            class = "qdm_calibration")
}
