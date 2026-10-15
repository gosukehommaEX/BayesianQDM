# Validation against Huang, Talukder, Han and Kuan (2019), "Quantitative
# decision-making in randomized Phase II studies with a time-to-event
# endpoint", Journal of Biopharmaceutical Statistics 29(1), 189-202,
# Table 1: probability of a Go decision in the pancreatic cancer example
# (Section 4) for true hazard ratios 0.5 to 0.8.
#
# Example: 2:1 randomization of 168 subjects (112 treatment, 56 control),
# uniform enrollment over 13 months, control median overall survival 8.5
# months (exponential), analysis at 127 events, and exponential loss to
# follow-up. Go if the observed hazard ratio is at most TV = 0.67 and the
# upper bound of its one-sided 80% confidence interval is at most
# MAV = 0.80. The article computes P(Go) from the normal approximation of
# the log hazard ratio with variance 1 / (p (1 - p) 127), p being the
# expected proportion of events in the treatment group (Section 2).
#
# BayesianQDM uses the exponential model with Gamma(0.01, 0.01) priors for
# the hazards, simulates the studies (analysis at the 127th event) and
# applies the posterior criterion P(HR < 0.67 | data) >= 0.5, the Bayesian
# counterpart of an observed hazard ratio at most 0.67; the second
# condition of the article practically never matters in this example
# (Python_reference/ref_huang2019_table1.py). The loss to follow-up of
# BayesianQDM is the same in both groups, so the rows of Table 1 with equal
# hazards are validated (0.003, 0.01 and 0.05). Differences up to 0.01
# beyond rounding are EXPLAINED (reason "huang2019_normal"). With hazard
# 0.05, 127 events are not expected even with unlimited follow-up; the
# article keeps the variance of 127 events, whereas the simulated studies
# analyze the events observed, so these rows are INFO (reason
# "huang2019_events").

validate_huang2019_table1 <- function(nsim = 100000L, seed = 2019L) {
  study <- "Huang et al. (2019), Table 1"
  pub <- list("0.003" = c(0.946, 0.727, 0.406, 0.170),
              "0.01" = c(0.946, 0.726, 0.406, 0.170),
              "0.05" = c(0.948, 0.728, 0.406, 0.170))
  hr <- c(0.5, 0.6, 0.7, 0.8)
  median_c <- 8.5
  g <- prior_gamma(0.01, 0.01)
  crit <- criterion_posterior(theta_tv = 0.67, theta_mav = 1,
                              gamma_go = 0.5, gamma_nogo = 0.1,
                              direction = "less")
  rows <- list()
  for (d in names(pub)) {
    des <- design_tte(112, 56, g, g, crit, effect = "hazard_ratio",
                      conduct = conduct_tte(accrual_time = c(0, 13),
                                            dropout_hazard = as.numeric(d),
                                            cut_events = 127))
    oc <- qdm_oc(des, data.frame(median_t = median_c / hr,
                                 median_c = median_c),
                 nsim = nsim, seed = seed)
    info <- if (d == "0.05") "huang2019_events" else ""
    for (k in seq_along(hr)) {
      rows <- c(rows, list(number_row(
        study, sprintf("P(Go | HR = %s)", format(hr[k])),
        sprintf("loss-to-follow-up hazard %s in both groups", d),
        sprintf("simulation (nsim = %d)", nsim), pub[[d]][k], oc$Go[k], 3,
        explain_tol = 0.01, reason_explained = "huang2019_normal",
        info = info
      )))
    }
  }
  do.call(rbind, rows)
}
