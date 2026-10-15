# Helpers for the validation of single published numbers (continuous
# endpoint). Each validate_*() function returns a data frame with one row per
# published number and the columns
#   study, quantity, setting, method, published, computed, decimals,
#   explain_tol, reason_explained, info
# and number_verdict() adds abs_diff, verdict and reason.
#
# Verdicts:
#   PASS       |computed - published| is at most half a unit of the last
#              published decimal (rounding).
#   EXPLAINED  the difference is at most rounding + explain_tol, where
#              explain_tol is set only for numbers that the article computed
#              by Monte Carlo simulation; the reason is registered in
#              run_validation.R.
#   FAIL       otherwise.
#   INFO       the published number could not be reproduced from the inputs
#              stated in the article; the row records both values (the
#              reason is registered in run_validation.R).

number_row <- function(study, quantity, setting, method, published, computed,
                       decimals, explain_tol = NA_real_,
                       reason_explained = "", info = "") {
  data.frame(study = study, quantity = quantity, setting = setting,
             method = method, published = published, computed = computed,
             decimals = decimals, explain_tol = explain_tol,
             reason_explained = reason_explained, info = info,
             stringsAsFactors = FALSE)
}

number_verdict <- function(res) {
  rounding <- 0.5 * 10 ^ (-res$decimals) + 1e-9
  diff <- abs(res$computed - res$published)
  explained <- !is.na(res$explain_tol) & diff <= rounding + res$explain_tol
  verdict <- ifelse(diff <= rounding, "PASS",
                    ifelse(explained, "EXPLAINED", "FAIL"))
  reason <- ifelse(verdict == "EXPLAINED", res$reason_explained, "")
  info <- nzchar(res$info)
  verdict[info] <- "INFO"
  reason[info] <- res$info[info]
  out <- res[, c("study", "quantity", "setting", "method", "published",
                 "computed")]
  out$abs_diff <- diff
  out$verdict <- verdict
  out$reason <- reason
  out
}

# Self-check: shifting every published number (except INFO rows) by 0.05
# must give FAIL, so that the comparison is not vacuous.
number_self_check <- function(raw) {
  checked <- !nzchar(raw$info)
  shifted <- raw
  shifted$published <- shifted$published + 0.05
  v <- number_verdict(shifted)$verdict[checked]
  list(n_checked = sum(checked), n_fail = sum(v == "FAIL"),
       ok = all(v == "FAIL"))
}
