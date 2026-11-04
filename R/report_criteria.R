#' Decision criteria section of a report
#'
#' @param design A design.
#'
#' @return A list of blocks: the basis of the decision, the regions (two
#'   endpoints with the joint rule), the Go and NoGo criteria, and the
#'   decision rule.
#' @noRd
report_criteria <- function(design) {
  crit <- design$criterion
  num <- report_num
  vec <- function(v) if (length(v) == 1L) num(v) else sprintf("(%s)", num(v))
  pred <- crit$type == "predictive"
  two <- crit$n_endpoints == 2L
  joint <- two && crit$rule == "joint"
  future <- if (design$control == "fixed") {
    sprintf("$m_t = %d$ subjects", crit$m_t)
  } else {
    sprintf(paste("$m_t = %d$ subjects in the treatment group and",
                  "$m_c = %d$ subjects in the control group"), crit$m_t,
            crit$m_c)
  }
  values <- if (pred) {
    sprintf("the null value%s $\\theta_{NULL} = %s$", if (two) "s" else "",
            vec(crit$theta_null))
  } else {
    sprintf(paste("the target value%s $\\theta_{TV} = %s$ and the minimum",
                  "acceptable value%s $\\theta_{MAV} = %s$"),
            if (two) "s" else "", vec(crit$theta_tv), if (two) "s" else "",
            vec(crit$theta_mav))
  }
  regions_par <- paste("divide the plane of %s into %s regions (the columns",
                       "for endpoint 1 and the rows for endpoint 2):")
  less <- crit$direction == "less"
  favorable <- if (!any(less)) {
    character(0)
  } else if (!two) {
    "Smaller effects are favorable."
  } else if (all(less)) {
    "Smaller effects are favorable for both endpoints."
  } else {
    sprintf("Smaller effects are favorable for endpoint %d.", which(less))
  }
  # The sentence that introduces the table of the regions comes last
  basis <- if (!two) {
    c(if (pred) {
      sprintf(paste("The decision will be based on the posterior predictive",
                    "probability, given the data $D$ of the study, of the",
                    "treatment effect $\\tilde\\theta$ that will be observed",
                    "in a future study with %s, with %s."), future, values)
    } else {
      sprintf(paste("The decision will be based on the posterior probability",
                    "of the treatment effect $\\theta$ given the data $D$ of",
                    "the study, with %s."), values)
    }, favorable)
  } else if (joint) {
    if (pred) {
      c(sprintf(paste("The decision will be based on the posterior",
                      "predictive probabilities, given the data $D$ of the",
                      "study, of regions of the treatment effects",
                      "$(\\tilde\\theta_1, \\tilde\\theta_2)$ that will be",
                      "observed in a future study with %s."), future),
        favorable,
        sprintf(paste("The null values $\\theta_{NULL} = %s$", regions_par),
                vec(crit$theta_null), "$(\\tilde\\theta_1, \\tilde\\theta_2)$",
                "four"))
    } else {
      c(paste("The decision will be based on the posterior probabilities of",
              "regions of the treatment effects $(\\theta_1, \\theta_2)$",
              "given the data $D$ of the study."),
        favorable,
        sprintf(paste("The target values $\\theta_{TV} = %s$ and the minimum",
                      "acceptable values $\\theta_{MAV} = %s$", regions_par),
                vec(crit$theta_tv), vec(crit$theta_mav),
                "$(\\theta_1, \\theta_2)$", "nine"))
    }
  } else {
    c(if (pred) {
      sprintf(paste("The decision will be based on the posterior predictive",
                    "probabilities, given the data $D$ of the study, of each",
                    "treatment effect that will be observed in a future study",
                    "with %s, with %s."), future, values)
    } else {
      sprintf(paste("The decision will be based on the posterior",
                    "probabilities of each treatment effect given the data",
                    "$D$ of the study, with %s."), values)
    },
    favorable,
    sprintf("The criteria of the two endpoints are combined with \"%s\".",
            crit$rule))
  }
  basis <- paste(basis, collapse = " ")
  ev <- report_events(crit)
  join <- if (two && crit$rule == "or") " or " else " and "
  criterion <- function(e, gamma) {
    gamma <- rep_len(gamma, length(e))
    paste(sprintf("$%s \\ge %s$", e, vapply(gamma, num, character(1))),
          collapse = join)
  }
  rule <- paste("The decision will be Go if only the Go criterion is met,",
                "NoGo if only the NoGo criterion is met, and Gray if neither",
                "is met; if both are met, the decision is called Miss.")
  # For one endpoint, the Go and NoGo regions are disjoint (their
  # probabilities add up to at most 1, and to 1 for a predictive criterion)
  if (!two && crit$gamma_go + crit$gamma_nogo > 1) {
    rule <- paste(rule, "Because $\\gamma_{go} + \\gamma_{nogo} > 1$, the",
                  "two criteria cannot both be met.")
  } else if (!two && pred) {
    rule <- paste(rule, "Because the two probabilities add up to 1 and",
                  "$\\gamma_{go} + \\gamma_{nogo} \\le 1$, a Gray decision",
                  "is not possible.")
  }
  c(list(report_block("heading", text = "Decision criteria"),
         report_block("par", text = basis)),
    if (joint) list(report_regions(crit)),
    list(report_block("list", items = c(
      paste("Go criterion:", criterion(ev$go, crit$gamma_go)),
      paste("NoGo criterion:", criterion(ev$nogo, crit$gamma_nogo))
    )),
    report_block("par", text = rule)))
}
