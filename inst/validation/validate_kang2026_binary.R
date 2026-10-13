# Validation against Kang, Yamaguchi and Han (2026), Supplemental Tables 2 to 5:
# operating characteristics of a single binary endpoint with posterior and
# posterior predictive criteria, several priors, sample sizes and power prior
# weights. The published values (three decimals) are in
# kang2026_supplement.csv, extracted from the supplementary Word file by
# Python_reference/extract_kang2026_supplement.py in the work folder.
#
# Verdicts:
#   PASS       all of Go, NoGo and Gray agree within 0.0005 (rounding).
#   EXPLAINED  (predictive criterion only) the difference is at most the total
#              probability of the outcomes whose exact predictive probability
#              is within 0.01 of a threshold, plus 0.0005. The article
#              computes the predictive probability by Monte Carlo
#              simulation, whose error can change the decisions for these
#              outcomes (reason "kang2026_mc").
#   FAIL       otherwise.
#
# This file defines validate_kang2026_binary(); run_validation.R calls it.

kang2026_design <- function(r) {
  criterion <- if (r$criterion == "posterior") {
    criterion_posterior(theta_tv = r$theta_tv, theta_mav = r$theta_mav,
                        gamma_go = r$gamma_go, gamma_nogo = r$gamma_nogo)
  } else {
    criterion_predictive(theta_null = r$theta_null, m_t = r$m_t,
                         m_c = r$m_c, gamma_go = r$gamma_go,
                         gamma_nogo = r$gamma_nogo)
  }
  external_t <- NULL
  external_c <- NULL
  if (r$n_e > 0) {
    external_t <- borrow(data_binary(r$n_e, r$y_et), r$a0_t)
    external_c <- borrow(data_binary(r$n_e, r$y_ec), r$a0_c)
  }
  design_binary(n_t = r$n_t, n_c = r$n_c,
                prior_t = prior_beta(r$a_t, r$b_t),
                prior_c = prior_beta(r$a_c, r$b_c),
                criterion = criterion, external_t = external_t,
                external_c = external_c)
}

kang2026_verdict <- function(res, rounding = 0.0005 + 1e-9) {
  diff <- pmax(abs(res$Go - res$pub_Go), abs(res$NoGo - res$pub_NoGo),
               abs(res$Gray - res$pub_Gray))
  verdict <- ifelse(diff <= rounding, "PASS",
                    ifelse(res$criterion == "predictive" &
                             diff <= res$sensitivity + rounding,
                           "EXPLAINED", "FAIL"))
  data.frame(max_abs_diff = diff, verdict = verdict,
             reason = ifelse(verdict == "EXPLAINED", "kang2026_mc", ""))
}

validate_kang2026_binary <- function(ref_file, delta = 0.01) {
  ref <- utils::read.csv(ref_file, stringsAsFactors = FALSE)
  key <- c("criterion", "a_t", "b_t", "a_c", "b_c", "n_t", "n_c", "m_t", "m_c",
           "n_e", "y_et", "y_ec", "a0_t", "a0_c")
  groups <- split(seq_len(nrow(ref)), do.call(paste, c(ref[key], sep = "|")))
  res <- vector("list", length(groups))
  for (g in seq_along(groups)) {
    idx <- groups[[g]]
    des <- kang2026_design(ref[idx[1], ])
    scen <- data.frame(pi_t = ref$pi_t[idx], pi_c = ref$pi_c[idx])
    oc <- qdm_oc(des, scen)
    # Probability of the outcomes whose region probability is close to a
    # threshold (sensitivity of the decisions to Monte Carlo error)
    s <- BayesianQDM:::oc_sample(des, scen, 1L)
    near <- abs(s$p_go - des$criterion$gamma_go) < delta |
      abs(s$p_nogo - des$criterion$gamma_nogo) < delta
    res[[g]] <- data.frame(
      table = ref$table[idx], row = ref$row[idx],
      criterion = ref$criterion[idx], setting = ref$setting[idx],
      a0_t = ref$a0_t[idx], a0_c = ref$a0_c[idx], effect = ref$effect[idx],
      Go = oc$Go, NoGo = oc$NoGo, Gray = oc$Gray, Miss = oc$Miss,
      pub_Go = ref$Go[idx], pub_NoGo = ref$NoGo[idx],
      pub_Gray = ref$Gray[idx],
      sensitivity = colSums(s$weights * near)
    )
  }
  res <- do.call(rbind, res)
  res <- res[order(res$table, res$row, res$criterion), ]
  rownames(res) <- NULL
  cbind(res, kang2026_verdict(res))
}

# Self-check: shifting every published Go probability by 0.02 must give FAIL
# for all rows of the posterior criterion, so the comparison is not vacuous.
kang2026_self_check <- function(res) {
  shifted <- res
  shifted$pub_Go <- shifted$pub_Go + 0.02
  v <- kang2026_verdict(shifted)$verdict
  post <- res$criterion == "posterior"
  list(n_posterior = sum(post), n_fail = sum(v[post] == "FAIL"),
       ok = all(v[post] == "FAIL"))
}
