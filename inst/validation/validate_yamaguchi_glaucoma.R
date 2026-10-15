# Validation against the glaucoma example of Yamaguchi, Homma, Takeda and
# Maruo (submitted), "Moment-matching approximation for t-differences in
# clinical trials with continuous endpoints", Section 5.2: treatment
# n = 36, mean -5.2, variance 15.62; control n = 31, mean -2.7, variance
# 13.28; vague priors. Published: P(Delta < 0 | data) = 0.9949 and
# P(Delta <= -2 | data) = 0.7011 (moment-matching approximation, four
# decimals). Both are checked with method = "NI" and method = "MM".

validate_yamaguchi_glaucoma <- function() {
  rows <- list()
  for (method in c("NI", "MM")) {
    des <- design_continuous(36, 31, prior_vague(), prior_vague(),
                             criterion_posterior(theta_tv = 0,
                                                 theta_mav = -2,
                                                 gamma_go = 0.8,
                                                 gamma_nogo = 0.6),
                             method = method)
    res <- qdm_decide(des, data_continuous(36, -5.2, sqrt(15.62)),
                      data_continuous(31, -2.7, sqrt(13.28)))
    study <- "Yamaguchi, Homma, Takeda and Maruo (submitted), Section 5.2"
    rows <- c(rows, list(
      number_row(study, "P(Delta < 0 | data)", "glaucoma example", method,
                 0.9949, 1 - res$prob_go, 4),
      number_row(study, "P(Delta <= -2 | data)", "glaucoma example", method,
                 0.7011, res$prob_nogo, 4)
    ))
  }
  do.call(rbind, rows)
}
