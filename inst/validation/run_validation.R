# Runs the validation scripts of BayesianQDM and writes the results and a
# summary to inst/validation/results/. The summary is generated from the
# results; do not edit it by hand.
#
# Usage: from the package root, after devtools::load_all() (or with the
# package attached),
#   source("inst/validation/run_validation.R")

local({
  val_dir <- file.path("inst", "validation")
  if (!dir.exists(val_dir)) {
    stop("Run this script from the root of the BayesianQDM package.",
         call. = FALSE)
  }
  out_dir <- file.path(val_dir, "results")
  dir.create(out_dir, showWarnings = FALSE)
  env <- environment()
  for (f in c("validate_kang2026_binary.R", "validation_helpers.R",
              "validate_yamaguchi_glaucoma.R", "validate_huang_table3.R",
              "validate_kang2026_table4.R", "validate_yamaguchi_binary2.R",
              "validate_yamaguchi_continuous2.R")) {
    sys.source(file.path(val_dir, f), envir = env)
  }
  count <- function(x, v) sum(x == v)

  # One binary endpoint: operating characteristics
  res <- validate_kang2026_binary(file.path(val_dir,
                                            "kang2026_supplement.csv"))
  utils::write.csv(res, file.path(out_dir, "kang2026_binary.csv"),
                   row.names = FALSE)
  self <- kang2026_self_check(res)
  tab <- do.call(rbind, lapply(split(res, list(res$table, res$criterion),
                                     drop = TRUE), function(d) {
    data.frame(table = d$table[1], criterion = d$criterion[1],
               checks = nrow(d), PASS = count(d$verdict, "PASS"),
               EXPLAINED = count(d$verdict, "EXPLAINED"),
               FAIL = count(d$verdict, "FAIL"),
               max_abs_diff = max(d$max_abs_diff))
  }))
  tab <- tab[order(tab$table, tab$criterion), ]

  # Two binary endpoints: operating characteristics
  res2 <- validate_yamaguchi_binary2(file.path(val_dir,
                                               "bivariate_binary2_tables.csv"))
  utils::write.csv(res2, file.path(out_dir, "yamaguchi_binary2.csv"),
                   row.names = FALSE)
  self2 <- yamaguchi_binary2_self_check(res2)
  tab2 <- do.call(rbind, lapply(split(res2, factor(res2$table,
                                                   unique(res2$table))),
                                function(d) {
    ne <- is.na(d$pub_Go)
    data.frame(table = d$table[1], checks = nrow(d),
               PASS = count(d$verdict, "PASS"),
               EXPLAINED = count(d$verdict, "EXPLAINED"),
               FAIL = count(d$verdict, "FAIL"), NE = sum(ne),
               max_abs_diff = max(d$max_abs_diff[!ne], -Inf),
               max_sensitivity = max(d$sensitivity[!ne], -Inf))
  }))
  tab2$label <- ifelse(tab2$table == "3", "Table 3",
                       paste("Supplemental Table", tab2$table))

  # Two continuous endpoints: operating characteristics
  nsim2 <- 10000L
  res3 <- validate_yamaguchi_continuous2(
    file.path(val_dir, "bivariate_continuous2_tables.csv"), nsim = nsim2
  )
  utils::write.csv(res3, file.path(out_dir, "yamaguchi_continuous2.csv"),
                   row.names = FALSE)
  self3 <- yamaguchi_continuous2_self_check(res3, nsim2)
  tab3 <- do.call(rbind, lapply(split(res3, factor(res3$table,
                                                   unique(res3$table))),
                                function(d) {
    ok <- d$verdict != "INFO"
    data.frame(table = d$table[1], checks = nrow(d),
               PASS = count(d$verdict, "PASS"),
               EXPLAINED = count(d$verdict, "EXPLAINED"),
               FAIL = count(d$verdict, "FAIL"),
               INFO = count(d$verdict, "INFO"),
               max_abs_diff = max(d$max_abs_diff[ok], -Inf),
               max_ratio = max(d$max_ratio[ok], -Inf))
  }))
  tab3$label <- ifelse(tab3$table == "4", "Table 4",
                       paste("Supplemental Table", tab3$table))
  info3 <- res3[res3$verdict == "INFO", ]

  # One continuous endpoint: published numbers
  raw <- rbind(validate_yamaguchi_glaucoma(), validate_huang_table3(),
               validate_kang2026_table4())
  num <- number_verdict(raw)
  utils::write.csv(num, file.path(out_dir, "continuous.csv"),
                   row.names = FALSE)
  self_num <- number_self_check(raw)
  tab_num <- do.call(rbind, lapply(split(num, factor(num$study,
                                                     unique(num$study))),
                                   function(d) {
    data.frame(study = d$study[1], checks = nrow(d),
               PASS = count(d$verdict, "PASS"),
               EXPLAINED = count(d$verdict, "EXPLAINED"),
               FAIL = count(d$verdict, "FAIL"),
               INFO = count(d$verdict, "INFO"),
               max_abs_diff = max(d$abs_diff[d$verdict != "INFO"],
                                  -Inf))
  }))
  info <- num[num$verdict == "INFO", ]
  not_pass <- num[num$verdict %in% c("EXPLAINED", "FAIL"), ]

  verdicts <- c(res$verdict, res2$verdict, res3$verdict, num$verdict)
  fmt_diff <- function(x) ifelse(is.finite(x), sprintf("%.4f", x), "-")
  lines <- c(
    "# Validation of BayesianQDM against published results",
    "",
    sprintf("Generated by `inst/validation/run_validation.R` with BayesianQDM %s.",
            as.character(utils::packageVersion("BayesianQDM"))),
    "",
    "## Kang, Yamaguchi and Han (2026), Supplemental Tables 2 to 5 (one binary endpoint)",
    "",
    "| Table | Criterion | Checks | PASS | EXPLAINED | FAIL | Largest difference |",
    "|---|---|---|---|---|---|---|",
    sprintf("| %d | %s | %d | %d | %d | %d | %.4f |", tab$table, tab$criterion,
            tab$checks, tab$PASS, tab$EXPLAINED, tab$FAIL, tab$max_abs_diff),
    "",
    sprintf("Largest probability of a Miss decision: %.3g.", max(res$Miss)),
    "",
    paste("## Yamaguchi, Hu, Huang, Han, Homma and Wang (submitted),",
          "Table 3 and Supplemental Tables S1 and S2 (two binary endpoints)"),
    "",
    paste("| Table | Checks | PASS | EXPLAINED | FAIL | Not evaluable (NE) |",
          "Largest difference | Largest sensitivity |"),
    "|---|---|---|---|---|---|---|---|",
    sprintf("| %s | %d | %d | %d | %d | %d | %.4f | %.4f |", tab2$label,
            tab2$checks, tab2$PASS, tab2$EXPLAINED, tab2$FAIL, tab2$NE,
            tab2$max_abs_diff, tab2$max_sensitivity),
    "",
    paste("One check per scenario (largest difference of Go, Gray and NoGo).",
          "The scenarios that the article reports as not evaluable (NE)",
          "are PASS when the probabilities of BayesianQDM are NA",
          "(correlation outside the feasible range)."),
    sprintf("Largest probability of a Miss decision: %.3g.",
            max(res2$Miss, na.rm = TRUE)),
    "",
    paste("## Yamaguchi, Hu, Huang, Han, Homma and Wang (submitted),",
          "Table 4 and Supplemental Table S3 (two continuous endpoints)"),
    "",
    paste("| Table | Checks | PASS | EXPLAINED | FAIL | INFO |",
          "Largest difference (except INFO) | Largest difference / tolerance |"),
    "|---|---|---|---|---|---|---|---|",
    sprintf("| %s | %d | %d | %d | %d | %d | %.4f | %.2f |", tab3$label,
            tab3$checks, tab3$PASS, tab3$EXPLAINED, tab3$FAIL, tab3$INFO,
            tab3$max_abs_diff, tab3$max_ratio),
    "",
    sprintf(paste("One check per scenario (Go, NoGo and Gray), with %d",
                  "simulated studies and the moment-matching approximation.",
                  "For Supplemental Table S3 the verdict uses the predictive",
                  "scale of the article; the probabilities of BayesianQDM",
                  "(corrected scale) are in the columns package_Go,",
                  "package_NoGo and package_Gray of the details."), nsim2),
    "",
    if (nrow(info3) > 0L) {
      c("INFO rows (published; computed with the scale of the article;",
        "BayesianQDM), Go / NoGo / Gray:",
        "",
        "| Table | Panel | theta1 | theta2 | Published | Article scale | BayesianQDM |",
        "|---|---|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %s | %.4f / %.4f / %.4f | %.4f / %.4f / %.4f | %.4f / %.4f / %.4f |",
                info3$table, info3$panel, format(info3$theta1),
                format(info3$theta2), info3$pub_Go, info3$pub_NoGo,
                info3$pub_Gray, info3$Go, info3$NoGo, info3$Gray,
                info3$package_Go, info3$package_NoGo, info3$package_Gray),
        "")
    },
    "## Published numbers for one continuous endpoint",
    "",
    "| Source | Checks | PASS | EXPLAINED | FAIL | INFO | Largest difference (except INFO) |",
    "|---|---|---|---|---|---|---|",
    sprintf("| %s | %d | %d | %d | %d | %d | %s |", tab_num$study,
            tab_num$checks, tab_num$PASS, tab_num$EXPLAINED, tab_num$FAIL,
            tab_num$INFO, fmt_diff(tab_num$max_abs_diff)),
    "",
    if (nrow(not_pass) > 0L) {
      c("Numbers that are not PASS (except INFO):",
        "",
        "| Source | Quantity | Setting | Method | Published | Computed | Verdict |",
        "|---|---|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %s | %s | %.4f | %s |", not_pass$study,
                not_pass$quantity, not_pass$setting, not_pass$method,
                format(not_pass$published), not_pass$computed,
                not_pass$verdict),
        "")
    },
    if (nrow(info) > 0L) {
      c("INFO rows (published and computed values):",
        "",
        "| Source | Quantity | Setting | Published | Computed |",
        "|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %s | %.4f |", info$study, info$quantity,
                info$setting, format(info$published), info$computed),
        "")
    },
    sprintf(paste("Total: %d checks, %d PASS, %d EXPLAINED, %d FAIL,",
                  "%d INFO."), length(verdicts), count(verdicts, "PASS"),
            count(verdicts, "EXPLAINED"), count(verdicts, "FAIL"),
            count(verdicts, "INFO")),
    "",
    "Registered reasons:",
    paste("- `kang2026_mc`: the article computes the posterior predictive",
          "probability by Monte Carlo simulation. The difference is at most",
          "the probability of the outcomes whose exact predictive probability",
          "is within 0.01 of a threshold (plus 0.0005 for rounding), whose",
          "decisions Monte Carlo error can change."),
    paste("- `bivariate_mc`: the article and BayesianQDM estimate the",
          "posterior probabilities of the regions by Monte Carlo",
          "simulation. The difference is at most the probability of the",
          "outcomes whose posterior probability of the Go or NoGo region is",
          "within 0.01 of its threshold (plus 0.00005 for rounding), whose",
          "decisions Monte Carlo error can change."),
    paste("- `bivariate_cont_mc`: the article estimates the decision",
          "probabilities of two continuous endpoints by simulating studies",
          "(10000 assumed; the number is not stated) and the region",
          "probabilities by Monte Carlo simulation. The difference is at",
          "most 4 standard errors of both simulations plus the probability",
          "of the studies whose probability of the Go or NoGo region is",
          "within 0.01 of its threshold (plus 0.00005 for rounding)."),
    paste("- `c1_predictive_scale`: the article computes the posterior",
          "predictive probabilities of two continuous endpoints with the",
          "scale ((1 + kappa) / kappa) Lambda / ((nu - 1) m) for the mean of",
          "m future subjects, which BayesianQDM corrects to",
          "(1 / kappa + 1 / m) Lambda / (nu - 1). With the scale of the",
          "article the published values are reproduced within the",
          "tolerance of `bivariate_cont_mc`."),
    paste("- `s3_rho08`: Supplemental Table S3, panel (c) (rho = 0.8), is",
          "not reproduced from the stated thresholds with either scale",
          "(for example, the probability of Go under the null scenario is",
          "0.08 to 0.09 with the scale of the article and about 0.06 with",
          "the corrected scale, instead of 0.0498). Monte Carlo draws",
          "shared across the simulated studies, as in version 0.1.0, do",
          "not explain the difference. Reported for information."),
    paste("- `huang_mc`: the article summarizes Monte Carlo samples of the",
          "posterior distributions (its Section 6.3); the difference from",
          "the exact value is at most 0.01 beyond rounding."),
    paste("- `kang2026_table4`: the case study of Kang et al. (2026,",
          "Section 6) is not reproduced from the inputs stated in the",
          "article (vague priors, standard deviations 25, 30 and 20); the",
          "article used Markov chain Monte Carlo sampling. Reported for",
          "information until the inputs are clarified."),
    "",
    sprintf(paste("Self-check: shifting the published Go probabilities by 0.02",
                  "gives FAIL in %d of %d posterior rows (%s)."),
            self$n_fail, self$n_posterior,
            if (self$ok) "PASS" else "FAIL"),
    sprintf(paste("Self-check: shifting the published Go probabilities of",
                  "two binary endpoints by 0.15 gives FAIL in %d of %d",
                  "evaluable scenarios (%s)."),
            self2$n_fail, self2$n_evaluable,
            if (self2$ok) "PASS" else "FAIL"),
    sprintf(paste("Self-check: shifting the published Go probabilities of",
                  "two continuous endpoints by 0.2 gives FAIL in %d of %d",
                  "scenarios other than INFO (%s)."),
            self3$n_fail, self3$n_checked,
            if (self3$ok) "PASS" else "FAIL"),
    sprintf(paste("Self-check: shifting the published continuous numbers by",
                  "0.05 gives FAIL in %d of %d rows other than INFO (%s)."),
            self_num$n_fail, self_num$n_checked,
            if (self_num$ok) "PASS" else "FAIL"),
    "",
    paste("Details: `inst/validation/results/kang2026_binary.csv`,",
          "`inst/validation/results/yamaguchi_binary2.csv`,",
          "`inst/validation/results/yamaguchi_continuous2.csv` and",
          "`inst/validation/results/continuous.csv`.")
  )
  writeLines(lines, file.path(out_dir, "summary.md"))
  message(sprintf(paste("Validation: %d checks, %d PASS, %d EXPLAINED,",
                        "%d FAIL, %d INFO; self-checks %s"),
                  length(verdicts), count(verdicts, "PASS"),
                  count(verdicts, "EXPLAINED"), count(verdicts, "FAIL"),
                  count(verdicts, "INFO"),
                  if (self$ok && self2$ok && self3$ok && self_num$ok) {
                    "PASS"
                  } else {
                    "FAIL"
                  }))
})
