# Validation against Yamaguchi, Homma, Takeda and Maruo (submitted), Table 3
# (case study 1) and Supplemental Tables S1 and S2 (case study 3, with and
# without external data): operating characteristics of two binary endpoints
# with a posterior criterion on the joint regions, for three correlations
# (rho_t = rho_c) and sixteen pairs of treatment effects. The published
# values (four decimals; NE for scenarios that are not evaluable) are in
# bivariate_binary2_tables.csv, extracted from the article and its
# supplement by Python_reference/extract_bivariate_tables.py in the work
# folder.
#
# Settings (article Section 3.1 and Supplemental Section S1.1): Dirichlet
# priors Dir(0.25, 0.25, 0.25, 0.25) in both groups;
#   Table 3   n_t = n_c = 10, control rates (0.15, 0.20), TV (0.15, 0.15),
#             MAV (0.10, 0.10), Go region 1, NoGo region 9;
#   S1, S2    n_t = n_c = 8, control rates (0.20, 0.15), TV (0.15, 0.25),
#             MAV (0.10, 0.15), Go regions 1 and 2, NoGo region 9; S1 adds
#             the external data (2, 1, 1, 2) with a0 = 0.7 (treatment) and
#             (14, 2, 3, 1) with a0 = 0.5 (control).
# The supplement states "15 patients per arm (nt = nc = 8)"; the tables are
# reproduced with 8 patients per group.
#
# Verdicts (one per scenario):
#   PASS       Go, Gray and NoGo agree within 0.00005 (rounding), or the
#              article reports the scenario as not evaluable (NE) and the
#              probabilities of BayesianQDM are NA.
#   EXPLAINED  the difference is at most the total probability of the
#              outcomes whose posterior probability of the Go or NoGo region
#              is within 0.01 of its threshold, plus 0.00005. The article and
#              BayesianQDM both estimate these posterior probabilities by
#              Monte Carlo simulation, so the decisions for these outcomes can
#              differ (reason "bivariate_mc").
#   FAIL       otherwise.
#
# This file defines validate_yamaguchi_binary2(); run_validation.R calls it.

yamaguchi_binary2_design <- function(table, gamma_go, gamma_nogo, nMC) {
  prior <- prior_dirichlet(rep(0.25, 4))
  if (table == "3") {
    criterion <- criterion_posterior(theta_tv = c(0.15, 0.15),
                                     theta_mav = c(0.10, 0.10),
                                     gamma_go = gamma_go,
                                     gamma_nogo = gamma_nogo,
                                     go_regions = 1L, nogo_regions = 9L)
    return(design_binary2(n_t = 10, n_c = 10, prior_t = prior,
                          prior_c = prior, criterion = criterion, nMC = nMC))
  }
  criterion <- criterion_posterior(theta_tv = c(0.15, 0.25),
                                   theta_mav = c(0.10, 0.15),
                                   gamma_go = gamma_go,
                                   gamma_nogo = gamma_nogo,
                                   go_regions = 1:2, nogo_regions = 9L)
  external_t <- NULL
  external_c <- NULL
  if (table == "S1") {
    external_t <- borrow(data_binary2(c(2, 1, 1, 2)), 0.7)
    external_c <- borrow(data_binary2(c(14, 2, 3, 1)), 0.5)
  }
  design_binary2(n_t = 8, n_c = 8, prior_t = prior, prior_c = prior,
                 criterion = criterion, external_t = external_t,
                 external_c = external_c, nMC = nMC)
}

yamaguchi_binary2_verdict <- function(res, rounding = 0.00005 + 1e-9) {
  diff <- pmax(abs(res$Go - res$pub_Go), abs(res$NoGo - res$pub_NoGo),
               abs(res$Gray - res$pub_Gray))
  ne <- is.na(res$pub_Go)
  verdict <- ifelse(diff <= rounding, "PASS",
                    ifelse(diff <= res$sensitivity + rounding, "EXPLAINED",
                           "FAIL"))
  verdict[ne] <- ifelse(is.na(res$Go[ne]), "PASS", "FAIL")
  verdict[!ne & is.na(diff)] <- "FAIL"
  data.frame(max_abs_diff = diff, verdict = verdict,
             reason = ifelse(verdict == "EXPLAINED", "bivariate_mc",
                             ifelse(ne & verdict == "PASS",
                                    "not evaluable in both", "")))
}

validate_yamaguchi_binary2 <- function(ref_file, nMC = 50000L, delta = 0.01,
                                       seed = 1L) {
  ref <- utils::read.csv(ref_file, stringsAsFactors = FALSE,
                         colClasses = c(table = "character"))
  control <- list("3" = c(0.15, 0.20), S1 = c(0.20, 0.15),
                  S2 = c(0.20, 0.15))
  res <- vector("list", length(control))
  for (g in seq_along(control)) {
    tab <- names(control)[g]
    d <- ref[ref$table == tab, ]
    pi_c <- control[[tab]]
    scen <- data.frame(pi_t1 = pi_c[1] + d$theta1,
                       pi_t2 = pi_c[2] + d$theta2, rho_t = d$rho,
                       pi_c1 = pi_c[1], pi_c2 = pi_c[2], rho_c = d$rho)
    # The posterior probabilities of the regions do not depend on the
    # thresholds or the scenario, so the outcomes are enumerated once per
    # table (with the thresholds of the first panel) and the operating
    # characteristics of each panel are computed from the same sample.
    des <- yamaguchi_binary2_design(tab, d$gamma_go[1], d$gamma_nogo[1],
                                    nMC)
    scen <- BayesianQDM:::check_scenarios(des, scen)
    s <- withCallingHandlers(
      BayesianQDM:::with_seed(seed, BayesianQDM:::oc_sample(des, scen, 1L)),
      warning = function(w) {
        if (grepl("not evaluable", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      })
    oc <- matrix(NA_real_, nrow(d), 4L)
    sens <- rep(NA_real_, nrow(d))
    for (k in unique(paste(d$gamma_go, d$gamma_nogo))) {
      rows <- which(paste(d$gamma_go, d$gamma_nogo) == k)
      gg <- d$gamma_go[rows[1]]
      gn <- d$gamma_nogo[rows[1]]
      oc[rows, ] <- BayesianQDM:::oc_from_sample(s, gg, gn)[rows, ]
      near <- abs(s$p_go - gg) < delta | abs(s$p_nogo - gn) < delta
      sens[rows] <- BayesianQDM:::outcome_prob(near, s)[rows]
    }
    res[[g]] <- data.frame(
      table = d$table, panel = d$panel, rho = d$rho,
      gamma_go = d$gamma_go, gamma_nogo = d$gamma_nogo,
      theta1 = d$theta1, theta2 = d$theta2,
      Go = oc[, 1L], NoGo = oc[, 2L], Gray = oc[, 3L], Miss = oc[, 4L],
      pub_Go = d$Go, pub_NoGo = d$NoGo, pub_Gray = d$Gray,
      sensitivity = sens, stringsAsFactors = FALSE
    )
  }
  res <- do.call(rbind, res)
  rownames(res) <- NULL
  cbind(res, yamaguchi_binary2_verdict(res))
}

# Self-check: shifting every published Go probability by 0.15 (about twice
# the largest sensitivity) must give FAIL for all evaluable scenarios, so the
# comparison is not vacuous.
yamaguchi_binary2_self_check <- function(res) {
  shifted <- res
  shifted$pub_Go <- shifted$pub_Go + 0.15
  v <- yamaguchi_binary2_verdict(shifted)$verdict
  ok_rows <- !is.na(res$pub_Go)
  list(n_evaluable = sum(ok_rows), n_fail = sum(v[ok_rows] == "FAIL"),
       ok = all(v[ok_rows] == "FAIL"),
       max_sensitivity = max(res$sensitivity[ok_rows]))
}
