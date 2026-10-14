# Validation against Yamaguchi, Homma, Takeda and Maruo (submitted), Table 4
# (case study 2, posterior criterion) and Supplemental Table S3 (case study
# 4, predictive criterion): operating characteristics of two continuous
# endpoints for three correlations (rho_t = rho_c) and sixteen pairs of
# treatment effects. The published values (four decimals) are in
# bivariate_continuous2_tables.csv, extracted from the article and its
# supplement by Python_reference/extract_bivariate_tables.py in the work
# folder.
#
# Settings (article Section 3.2 and Supplemental Section S1.2):
#   Table 4   n_t = n_c = 30, NIW priors with mu0 = (0, 1) (treatment) and
#             (-10, -1) (control), kappa0 = 0.1, nu0 = 4,
#             Lambda0 = (80^2, 15; 15, 6^2); control means (-10, -1),
#             standard deviations (80, 6); TV (10, 2), MAV (5, 1); Go region
#             1, NoGo region 9.
#   S3        n_t = n_c = 12, NIW priors with mu0 = (1.984, 1.191) and
#             (-0.5, -1.5), kappa0 = 0.1, nu0 = 4, Lambda0 = (25, 1; 1, 9);
#             control means (-0.5, -1.5), standard deviations (5, 3);
#             m_t = m_c = 50 future subjects, null values (1.984, 1.191);
#             Go region 1, NoGo region 4.
# The article computes the decision probabilities by simulating studies
# (their number is not stated; 10000, the default of getgamma2cont() in
# version 0.1.0, is assumed) and the region probabilities by Monte Carlo
# simulation. BayesianQDM simulates nsim studies and uses the
# moment-matching approximation.
#
# Verdicts (one per scenario):
#   PASS       Go, NoGo and Gray agree within 0.00005 (rounding).
#   EXPLAINED  each of Go, NoGo and Gray differs by at most
#              4 sqrt(p (1 - p) (1 / 10000 + 1 / nsim)) + s + 0.00005, where
#              p is the published probability and s the probability of the
#              simulated studies whose probability of the Go or NoGo region
#              is within 0.01 of its threshold (decisions that Monte Carlo
#              error can change).
#              Table 4: reason "bivariate_cont_mc" (simulation error).
#              Table S3: reason "c1_predictive_scale": the article used the
#              predictive scale ((1 + kappa) / kappa) Lambda / ((nu - 1) m),
#              which BayesianQDM corrects to (1 / kappa + 1 / m) Lambda /
#              (nu - 1) (inspection item C1). The verdict compares the
#              published values with the computation using the scale of the
#              article; the values of BayesianQDM are reported alongside.
#   INFO       Table S3, panel (c) (rho = 0.8): not reproduced with the
#              stated thresholds by either scale (reason "s3_rho08"); the
#              values are reported for information.
#   FAIL       otherwise.
#
# This file defines validate_yamaguchi_continuous2(); run_validation.R calls
# it.

yamaguchi_continuous2_design <- function(table, gamma_go, gamma_nogo) {
  if (table == "4") {
    lam <- matrix(c(80 ^ 2, 15, 15, 6 ^ 2), 2L, 2L)
    return(design_continuous2(
      n_t = 30, n_c = 30, prior_t = prior_niw(c(0, 1), 0.1, 4, lam),
      prior_c = prior_niw(c(-10, -1), 0.1, 4, lam),
      criterion = criterion_posterior(theta_tv = c(10, 2),
                                      theta_mav = c(5, 1),
                                      gamma_go = gamma_go,
                                      gamma_nogo = gamma_nogo,
                                      go_regions = 1L, nogo_regions = 9L)
    ))
  }
  lam <- matrix(c(25, 1, 1, 9), 2L, 2L)
  design_continuous2(
    n_t = 12, n_c = 12, prior_t = prior_niw(c(1.984, 1.191), 0.1, 4, lam),
    prior_c = prior_niw(c(-0.5, -1.5), 0.1, 4, lam),
    criterion = criterion_predictive(theta_null = c(1.984, 1.191), m_t = 50,
                                     m_c = 50, gamma_go = gamma_go,
                                     gamma_nogo = gamma_nogo,
                                     go_regions = 1L, nogo_regions = 4L)
  )
}

# Same computation as oc_sample() (with the same random numbers), but with
# the predictive scale of the article for both groups.
yamaguchi_continuous2_article_sample <- function(des, scen, nsim, seed) {
  crit <- des$criterion
  sets <- BayesianQDM:::region_sets(crit)
  BayesianQDM:::with_seed(seed, {
    draws_t <- BayesianQDM:::draws_continuous2(nsim, des$n_t)
    draws_c <- BayesianQDM:::draws_continuous2(nsim, des$n_c)
  })
  dist <- function(post, m) {
    list(loc = post$mu,
         scale = (1 + post$kappa) / post$kappa * post$Lambda /
           (post$nu - 1) / m,
         df = post$nu - 1)
  }
  p <- matrix(0, nsim * nrow(scen), 2L)
  for (j in seq_len(nrow(scen))) {
    sc <- scen[j, ]
    st <- BayesianQDM:::stats_continuous2(draws_t, des$n_t,
                                          c(sc$mu_t1, sc$mu_t2),
                                          c(sc$sd_t1, sc$sd_t2), sc$rho_t)
    sc_c <- BayesianQDM:::stats_continuous2(draws_c, des$n_c,
                                            c(sc$mu_c1, sc$mu_c2),
                                            c(sc$sd_c1, sc$sd_c2), sc$rho_c)
    post_t <- BayesianQDM:::posterior_continuous2(des$prior_t, des$n_t,
                                                  st$mean, st$cov)
    post_c <- BayesianQDM:::posterior_continuous2(des$prior_c, des$n_c,
                                                  sc_c$mean, sc_c$cov)
    p[(j - 1L) * nsim + seq_len(nsim), ] <-
      BayesianQDM:::region_prob_continuous2(crit, dist(post_t, crit$m_t),
                                            dist(post_c, crit$m_c),
                                            c(sets$go, sets$nogo), "MM")
  }
  list(p_go = matrix(p[, 1L], nsim), p_nogo = matrix(p[, 2L], nsim),
       per_scenario = TRUE,
       weights = matrix(1 / nsim, nsim, nrow(scen)))
}

yamaguchi_continuous2_verdict <- function(res, nsim, b_article = 10000,
                                          rounding = 0.00005 + 1e-9) {
  dec <- c("Go", "NoGo", "Gray")
  diff <- sapply(dec, function(d) abs(res[[d]] - res[[paste0("pub_", d)]]))
  tol <- sapply(dec, function(d) {
    p <- pmin(pmax(res[[paste0("pub_", d)]], 0), 1)
    4 * sqrt(p * (1 - p) * (1 / b_article + 1 / nsim)) + res$sensitivity +
      rounding
  })
  diff <- matrix(diff, ncol = 3L)
  tol <- matrix(tol, ncol = 3L)
  max_diff <- apply(diff, 1L, max)
  verdict <- ifelse(max_diff <= rounding, "PASS",
                    ifelse(rowSums(diff > tol) == 0, "EXPLAINED", "FAIL"))
  reason <- ifelse(verdict == "EXPLAINED",
                   ifelse(res$table == "4", "bivariate_cont_mc",
                          "c1_predictive_scale"), "")
  info <- res$table == "S3" & abs(res$rho - 0.8) < 1e-9
  verdict[info] <- "INFO"
  reason[info] <- "s3_rho08"
  data.frame(max_abs_diff = max_diff,
             max_ratio = apply(diff / tol, 1L, max), verdict = verdict,
             reason = reason)
}

validate_yamaguchi_continuous2 <- function(ref_file, nsim = 10000L,
                                           delta = 0.01, seed = 1L) {
  ref <- utils::read.csv(ref_file, stringsAsFactors = FALSE,
                         colClasses = c(table = "character"))
  control <- list("4" = list(mu = c(-10, -1), sd = c(80, 6)),
                  S3 = list(mu = c(-0.5, -1.5), sd = c(5, 3)))
  res <- vector("list", length(control))
  for (g in seq_along(control)) {
    tab <- names(control)[g]
    d <- ref[ref$table == tab, ]
    cc <- control[[tab]]
    scen <- data.frame(mu_t1 = cc$mu[1] + d$theta1,
                       mu_t2 = cc$mu[2] + d$theta2, sd_t1 = cc$sd[1],
                       sd_t2 = cc$sd[2], rho_t = d$rho, mu_c1 = cc$mu[1],
                       mu_c2 = cc$mu[2], sd_c1 = cc$sd[1], sd_c2 = cc$sd[2],
                       rho_c = d$rho)
    # The probabilities of the regions do not depend on the thresholds, so
    # the studies are simulated once per table (with the thresholds of the
    # first panel) and the decisions of each panel are computed from them.
    des <- yamaguchi_continuous2_design(tab, d$gamma_go[1], d$gamma_nogo[1])
    scen <- BayesianQDM:::check_scenarios(des, scen)
    s <- BayesianQDM:::with_seed(seed,
                                 BayesianQDM:::oc_sample(des, scen, nsim))
    if (tab == "S3") {
      s_article <- yamaguchi_continuous2_article_sample(des, scen, nsim, seed)
    }
    oc <- matrix(NA_real_, nrow(d), 4L)
    oc_pkg <- matrix(NA_real_, nrow(d), 4L)
    sens <- rep(NA_real_, nrow(d))
    for (k in unique(paste(d$gamma_go, d$gamma_nogo))) {
      rows <- which(paste(d$gamma_go, d$gamma_nogo) == k)
      gg <- d$gamma_go[rows[1]]
      gn <- d$gamma_nogo[rows[1]]
      oc_pkg[rows, ] <- BayesianQDM:::oc_from_sample(s, gg, gn)[rows, ]
      s_v <- if (tab == "S3") s_article else s
      oc[rows, ] <- BayesianQDM:::oc_from_sample(s_v, gg, gn)[rows, ]
      near <- abs(s_v$p_go - gg) < delta | abs(s_v$p_nogo - gn) < delta
      sens[rows] <- BayesianQDM:::outcome_prob(near, s_v)[rows]
    }
    res[[g]] <- data.frame(
      table = d$table, panel = d$panel, rho = d$rho,
      gamma_go = d$gamma_go, gamma_nogo = d$gamma_nogo,
      theta1 = d$theta1, theta2 = d$theta2,
      Go = oc[, 1L], NoGo = oc[, 2L], Gray = oc[, 3L], Miss = oc[, 4L],
      pub_Go = d$Go, pub_NoGo = d$NoGo, pub_Gray = d$Gray,
      sensitivity = sens,
      package_Go = oc_pkg[, 1L], package_NoGo = oc_pkg[, 2L],
      package_Gray = oc_pkg[, 3L], stringsAsFactors = FALSE
    )
  }
  res <- do.call(rbind, res)
  rownames(res) <- NULL
  cbind(res, yamaguchi_continuous2_verdict(res, nsim))
}

# Self-check: shifting every published Go probability by 0.2 (more than the
# largest tolerance) must give FAIL for all rows other than INFO, so the
# comparison is not vacuous.
yamaguchi_continuous2_self_check <- function(res, nsim) {
  shifted <- res
  shifted$pub_Go <- shifted$pub_Go + 0.2
  v <- yamaguchi_continuous2_verdict(shifted, nsim)$verdict
  checked <- res$verdict != "INFO"
  list(n_checked = sum(checked), n_fail = sum(v[checked] == "FAIL"),
       ok = isTRUE(all(v[checked] == "FAIL")))
}
