#' Search for the sample size
#'
#' Finds the smallest sample size, among candidates, for which the
#' probabilities of a false Go decision and a false NoGo decision are
#' controlled, with the thresholds of the criterion of the design fixed. It
#' is the counterpart of \code{\link{qdm_calibrate}}, which fixes the sample
#' size and chooses the thresholds.
#'
#' For each candidate, the design with the sample sizes \code{n_t[k]} and
#' \code{n_c[k]} (set by \code{\link{qdm_sample_size}}) gives the
#' probability that the Go criterion is met under \code{scenario_go}, a
#' scenario under which a Go decision would be wrong (typically with no or a
#' small treatment effect), and the probability that the NoGo criterion is
#' met under \code{scenario_nogo}, a scenario under which a NoGo decision
#' would be wrong (typically with the expected treatment effect), computed
#' as in
#' \code{\link{qdm_calibrate}}. A candidate meets the targets when the first
#' probability is at most \code{target_go} and the second at most
#' \code{target_nogo}.
#'
#' The probabilities need not decrease as the sample size increases: for a
#' binary endpoint they change in steps, so that a larger sample size can
#' fail the targets that a smaller one meets. The result therefore gives both
#' the smallest candidate that meets the targets (\code{n_t} and \code{n_c})
#' and the smallest candidate from which every larger candidate also meets
#' them (\code{n_t_stable} and \code{n_c_stable}).
#'
#' With a seed, the simulated studies of all candidates start from the same
#' random numbers, so the probabilities change smoothly with the sample size.
#' For a continuous or time-to-event endpoint, many candidates and the
#' default \code{nsim} can take long; \code{method = "MM"} (continuous
#' endpoints) or a smaller \code{nsim} gives a first answer quickly.
#'
#' @param design A design created by \code{\link{design_binary}},
#'   \code{\link{design_continuous}}, \code{\link{design_binary2}},
#'   \code{\link{design_continuous2}} or \code{\link{design_tte}}, or a design
#'   specification started by \code{\link{qdm_endpoint}}. Its sample sizes are
#'   replaced by the candidates.
#' @param scenario_go,scenario_nogo,target_go,target_nogo As in
#'   \code{\link{qdm_calibrate}}.
#' @param n_t Candidate numbers of subjects in the treatment group (distinct
#'   whole numbers, in any order; they are sorted).
#' @param n_c Candidate numbers of subjects in the concurrent control group,
#'   one per candidate of \code{n_t} or a single number for all, or
#'   \code{NULL} (the default) for as many as in the treatment group with a
#'   concurrent control group and 0 without one. A ratio of allocation is
#'   given by \code{n_c}, for example \code{n_c = n_t / 2} for even
#'   \code{n_t}.
#' @param nsim,seed As in \code{\link{qdm_oc}}.
#'
#' @return An object of class \code{"qdm_n_search"}: a list with elements
#'   \describe{
#'     \item{n_t, n_c}{The smallest candidate that meets the targets
#'       (\code{NA} if none does).}
#'     \item{n_t_stable, n_c_stable}{The smallest candidate from which every
#'       larger candidate also meets the targets (\code{NA} if the largest
#'       candidate does not meet them).}
#'     \item{prob_go, prob_nogo}{The probabilities that the Go criterion is
#'       met under \code{scenario_go} and that the NoGo criterion is met under
#'       \code{scenario_nogo}, for \code{n_t} and \code{n_c}.}
#'     \item{target_go, target_nogo}{The targets.}
#'     \item{curve}{A data frame with the columns \code{n_t}, \code{n_c},
#'       \code{prob_go}, \code{prob_nogo} and \code{met} (whether both targets
#'       are met), one row per candidate.}
#'     \item{oc}{A data frame of the decision probabilities with \code{n_t}
#'       and \code{n_c} under the two scenarios (\code{NULL} if no candidate
#'       meets the targets).}
#'     \item{design}{The design with \code{n_t} and \code{n_c} (\code{NULL}
#'       if no candidate meets the targets).}
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
#' res <- qdm_find_n(des,
#'                   scenario_go = data.frame(pi_t = 0.4, pi_c = 0.4),
#'                   scenario_nogo = data.frame(pi_t = 0.7, pi_c = 0.4),
#'                   target_go = 0.1, target_nogo = 0.1,
#'                   n_t = seq(6, 30, by = 2))
#' res
#' plot(res)
#' res$design
#' @export
qdm_find_n <- function(design, scenario_go, scenario_nogo, target_go,
                       target_nogo, n_t, n_c = NULL, nsim = 10000L,
                       seed = NULL) {
  design <- check_design(design)
  scenario_go <- check_scenarios(design, scenario_go, "scenario_go",
                                 one_row = TRUE, feasible = TRUE)
  scenario_nogo <- check_scenarios(design, scenario_nogo, "scenario_nogo",
                                   one_row = TRUE, feasible = TRUE)
  check_numeric(target_go, "target_go", lower = 0, upper = 1,
                lower_open = TRUE, upper_open = TRUE)
  check_numeric(target_nogo, "target_nogo", lower = 0, upper = 1,
                lower_open = TRUE, upper_open = TRUE)
  check_numeric(n_t, "n_t", len = NULL, lower = 1, integer = TRUE)
  n_t <- as.integer(round(n_t))
  if (anyDuplicated(n_t) > 0L) {
    stop("'n_t' must not contain repeated values.", call. = FALSE)
  }
  if (is.null(n_c)) {
    n_c <- if (design$control == "concurrent") n_t else 0L
  } else {
    check_numeric(n_c, "n_c", len = NULL, lower = 0, integer = TRUE)
    if (!(length(n_c) %in% c(1L, length(n_t)))) {
      stop("'n_c' must have length 1 or the length of 'n_t'.", call. = FALSE)
    }
  }
  n_c <- rep_len(as.integer(round(n_c)), length(n_t))
  ord <- order(n_t)
  n_t <- n_t[ord]
  n_c <- n_c[ord]
  check_numeric(nsim, "nsim", lower = 1, upper = .Machine$integer.max,
                integer = TRUE)
  if (!is.null(seed)) {
    check_numeric(seed, "seed", lower = -.Machine$integer.max,
                  upper = .Machine$integer.max, integer = TRUE)
  }
  cols <- scenario_columns(design)$names
  scenarios <- rbind(scenario_go[, cols, drop = FALSE],
                     scenario_nogo[, cols, drop = FALSE])
  rownames(scenarios) <- NULL
  crit <- design$criterion
  spec <- design_spec(design)
  candidates <- lapply(seq_along(n_t), function(k) {
    des_k <- qdm_build(qdm_sample_size(spec, n_t = n_t[k], n_c = n_c[k]))
    sample <- with_seed(seed, oc_sample(des_k, scenarios, nsim))
    go <- meets_criterion(sample$p_go, crit$gamma_go, crit$rule)
    nogo <- meets_criterion(sample$p_nogo, crit$gamma_nogo, crit$rule)
    list(design = des_k, method = sample$method,
         prob_go = outcome_prob(go, sample)[1L],
         prob_nogo = outcome_prob(nogo, sample)[2L],
         probs = oc_from_sample(sample, crit$gamma_go, crit$gamma_nogo,
                                crit$rule))
  })
  prob_go <- vapply(candidates, function(r) r$prob_go, numeric(1))
  prob_nogo <- vapply(candidates, function(r) r$prob_nogo, numeric(1))
  met <- prob_go <= target_go + 1e-12 & prob_nogo <= target_nogo + 1e-12
  first <- which(met)[1L]
  stable <- if (!met[length(met)]) {
    NA_integer_
  } else if (all(met)) {
    1L
  } else {
    max(which(!met)) + 1L
  }
  if (is.na(first)) {
    warning(sprintf(paste("No candidate sample size gives a probability of",
                          "meeting the Go criterion of at most %s under",
                          "'scenario_go' and of meeting the NoGo criterion",
                          "of at most %s under 'scenario_nogo'."),
                    format(target_go), format(target_nogo)), call. = FALSE)
  }
  oc <- NULL
  new_design <- NULL
  if (!is.na(first)) {
    oc <- cbind(data.frame(scenario = c("scenario_go", "scenario_nogo")),
                scenarios, as.data.frame(candidates[[first]]$probs))
    rownames(oc) <- NULL
    warn_miss(oc$Miss, crit$rule)
    new_design <- candidates[[first]]$design
  }
  structure(list(n_t = n_t[first], n_c = n_c[first],
                 n_t_stable = n_t[stable], n_c_stable = n_c[stable],
                 prob_go = prob_go[first], prob_nogo = prob_nogo[first],
                 target_go = target_go, target_nogo = target_nogo,
                 curve = data.frame(n_t = n_t, n_c = n_c, prob_go = prob_go,
                                    prob_nogo = prob_nogo, met = met),
                 oc = oc, design = new_design, scenario_go = scenario_go,
                 scenario_nogo = scenario_nogo,
                 method = candidates[[length(candidates)]]$method),
            class = "qdm_n_search")
}
