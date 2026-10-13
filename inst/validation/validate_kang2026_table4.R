# Comparison with Kang, Yamaguchi and Han (2026), Table 4 (case study in
# Duchenne muscular dystrophy, Section 6). Inputs stated in Section 6:
# TV = 30 and MAV = 15; 12 subjects per group in the PoC expansion cohort
# with means 18 (treatment) and -12 (control) and standard deviation 25;
# external treatment data of 6 subjects (mean 20, SD 30, a0t = 0.5);
# external control data of 30 subjects (mean -12, SD 20, a0c = 0.2 to 0.8);
# vague priors.
#
# With these inputs, P(theta > TV) is 0.5 without external data (the
# observed difference equals TV), close to the published 0.499, but the
# other published values (P(theta > TV) of 0.598 to 0.640 and
# P(theta > MAV) of 0.999 to 1.000) are not reproduced; the article used
# Markov chain Monte Carlo sampling (Section 4.2.3). All rows are therefore
# reported as INFO (reason "kang2026_table4") until the inputs are
# clarified.

validate_kang2026_table4 <- function() {
  study <- "Kang et al. (2026), Table 4"
  analyses <- data.frame(
    setting = c("PoC data only", "External treatment data",
                sprintf("External treatment and control data, a0c = %s",
                        c(0.2, 0.4, 0.6, 0.8))),
    a0t = c(NA, 0.5, 0.5, 0.5, 0.5, 0.5),
    a0c = c(NA, NA, 0.2, 0.4, 0.6, 0.8),
    pub_tv = c(0.499, 0.598, 0.620, 0.628, 0.636, 0.640),
    pub_mav = c(0.999, 0.999, 1.000, 1.000, 1.000, 1.000),
    stringsAsFactors = FALSE
  )
  rows <- list()
  for (i in seq_len(nrow(analyses))) {
    a <- analyses[i, ]
    ext_t <- if (is.na(a$a0t)) NULL else
      borrow(data_continuous(6, 20, 30), a$a0t)
    ext_c <- if (is.na(a$a0c)) NULL else
      borrow(data_continuous(30, -12, 20), a$a0c)
    des <- design_continuous(12, 12, prior_vague(), prior_vague(),
                             criterion_posterior(theta_tv = 30,
                                                 theta_mav = 15,
                                                 gamma_go = 0.6,
                                                 gamma_nogo = 0.7),
                             external_t = ext_t, external_c = ext_c)
    res <- qdm_decide(des, data_continuous(12, 18, 25),
                      data_continuous(12, -12, 25))
    rows <- c(rows, list(
      number_row(study, "P(theta > TV | data)", a$setting, "NI", a$pub_tv,
                 res$prob_go, 3, info = "kang2026_table4"),
      number_row(study, "P(theta > MAV | data)", a$setting, "NI",
                 a$pub_mav, 1 - res$prob_nogo, 3, info = "kang2026_table4")
    ))
  }
  do.call(rbind, rows)
}
