rule_text <- paste("The decision will be Go if only the Go criterion is met,",
                   "NoGo if only the NoGo criterion is met, and Gray if",
                   "neither is met; if both are met, the decision is called",
                   "Miss.")

test_that("report_criteria() describes the posterior criterion of one endpoint", {
  expect_identical(
    report_criteria(ref_design("D1")),
    list(list(type = "heading", text = "Decision criteria"),
         list(type = "par",
              text = paste("The decision will be based on the posterior",
                           "probability of the treatment effect $\\theta$",
                           "given the data $D$ of the study, with the target",
                           "value $\\theta_{TV} = 0.2$ and the minimum",
                           "acceptable value $\\theta_{MAV} = 0.1$.")),
         list(type = "list",
              items = c("Go criterion: $\\Pr(\\theta > 0.2 \\mid D) \\ge 0.5$",
                        paste("NoGo criterion: $\\Pr(\\theta \\le 0.1 \\mid",
                              "D) \\ge 0.8$"))),
         list(type = "par",
              text = paste(rule_text, "Because $\\gamma_{go} + \\gamma_{nogo}",
                           "> 1$, the two criteria cannot both be met.")))
  )
  # gamma_go + gamma_nogo <= 1: no statement on Miss decisions
  expect_identical(report_criteria(ref_design("D6"))[[4L]][["text"]],
                   rule_text)
  # Smaller effects are favorable
  expect_match(report_criteria(ref_design_hr("H1"))[[2L]][["text"]],
               paste("with the target value $\\theta_{TV} = 0.7$ and the",
                     "minimum acceptable value $\\theta_{MAV} = 1$. Smaller",
                     "effects are favorable."), fixed = TRUE)
})

test_that("report_criteria() describes the predictive criterion of one endpoint", {
  res <- report_criteria(ref_design("D5"))
  expect_identical(
    res[[2L]][["text"]],
    paste("The decision will be based on the posterior predictive",
          "probability, given the data $D$ of the study, of the treatment",
          "effect $\\tilde\\theta$ that will be observed in a future study",
          "with $m_t = 30$ subjects in the treatment group and $m_c = 20$",
          "subjects in the control group, with the null value",
          "$\\theta_{NULL} = 0.05$.")
  )
  expect_identical(
    res[[3L]][["items"]],
    c("Go criterion: $\\Pr(\\tilde\\theta > 0.05 \\mid D) \\ge 0.9$",
      "NoGo criterion: $\\Pr(\\tilde\\theta \\le 0.05 \\mid D) \\ge 0.7$")
  )
  expect_match(report_criteria(ref_design("D9"))[[2L]][["text"]],
               "observed in a future study with $m_t = 10$ subjects, with",
               fixed = TRUE)
  unif <- prior_beta(1, 1)
  des <- design_binary(10, 10, unif, unif,
                       criterion_predictive(0.05, m_t = 30, m_c = 20,
                                            gamma_go = 0.6, gamma_nogo = 0.3))
  expect_identical(report_criteria(des)[[4L]][["text"]],
                   paste(rule_text, "Because the two probabilities add up to",
                         "1 and $\\gamma_{go} + \\gamma_{nogo} \\le 1$, a Gray",
                         "decision is not possible."))
})

test_that("report_criteria() describes the regions of two endpoints", {
  des <- ref_design_bin2("B6")
  res <- report_criteria(des)
  expect_length(res, 5L)
  expect_identical(
    res[[2L]][["text"]],
    paste("The decision will be based on the posterior probabilities of",
          "regions of the treatment effects $(\\theta_1, \\theta_2)$ given the",
          "data $D$ of the study. The target values",
          "$\\theta_{TV} = (0.15, 0.25)$ and the minimum acceptable values",
          "$\\theta_{MAV} = (0.1, 0.15)$ divide the plane of",
          "$(\\theta_1, \\theta_2)$ into nine regions (the columns for",
          "endpoint 1 and the rows for endpoint 2):")
  )
  expect_identical(res[[3L]], report_regions(des$criterion))
  expect_identical(
    res[[4L]][["items"]],
    c(paste("Go criterion: $\\Pr\\{(\\theta_1, \\theta_2) \\in R_{1} \\cup",
            "R_{2} \\mid D\\} \\ge 0.453$"),
      paste("NoGo criterion: $\\Pr\\{(\\theta_1, \\theta_2) \\in R_{9} \\mid",
            "D\\} \\ge 0.229$"))
  )
  expect_identical(res[[5L]][["text"]], rule_text)
  expect_identical(
    report_criteria(ref_design_bin2("B1"))[[2L]][["text"]],
    paste("The decision will be based on the posterior predictive",
          "probabilities, given the data $D$ of the study, of regions of the",
          "treatment effects $(\\tilde\\theta_1, \\tilde\\theta_2)$ that will",
          "be observed in a future study with $m_t = 4$ subjects in the",
          "treatment group and $m_c = 3$ subjects in the control group. The",
          "null values $\\theta_{NULL} = (0.1, 0)$ divide the plane of",
          "$(\\tilde\\theta_1, \\tilde\\theta_2)$ into four regions (the",
          "columns for endpoint 1 and the rows for endpoint 2):")
  )
})

test_that("report_criteria() describes the marginal rules of two endpoints", {
  res <- report_criteria(ref_design_bin2("B4"))
  expect_length(res, 4L)
  expect_identical(
    res[[2L]][["text"]],
    paste("The decision will be based on the posterior probabilities of each",
          "treatment effect given the data $D$ of the study, with the target",
          "values $\\theta_{TV} = (0.15, 0.15)$ and the minimum acceptable",
          "values $\\theta_{MAV} = (0.1, 0.1)$. The criteria of the two",
          "endpoints are combined with \"and\".")
  )
  expect_identical(
    res[[3L]][["items"]],
    c(paste("Go criterion: $\\Pr(\\theta_1 > 0.15 \\mid D) \\ge 0.5$ and",
            "$\\Pr(\\theta_2 > 0.15 \\mid D) \\ge 0.4$"),
      paste("NoGo criterion: $\\Pr(\\theta_1 \\le 0.1 \\mid D) \\ge 0.6$ and",
            "$\\Pr(\\theta_2 \\le 0.1 \\mid D) \\ge 0.5$"))
  )
  expect_identical(
    report_criteria(ref_design_bin2("B2"))[[2L]][["text"]],
    paste("The decision will be based on the posterior predictive",
          "probabilities, given the data $D$ of the study, of each treatment",
          "effect that will be observed in a future study with $m_t = 4$",
          "subjects in the treatment group and $m_c = 3$ subjects in the",
          "control group, with the null values $\\theta_{NULL} = (0.1, 0)$.",
          "The criteria of the two endpoints are combined with \"and\".")
  )
  vague <- prior_vague()
  des <- design_continuous2(12, 12, vague, vague,
                            crit_cont2(c(0.6, 0.5), c(0.5, 0.6), rule = "or"))
  expect_identical(
    report_criteria(des)[[3L]][["items"]][1L],
    paste("Go criterion: $\\Pr(\\theta_1 > 10 \\mid D) \\ge 0.6$ or",
          "$\\Pr(\\theta_2 > 2 \\mid D) \\ge 0.5$")
  )
})

test_that("report_criteria() says for which endpoints smaller effects are favorable", {
  vague <- prior_vague()
  one <- design_continuous2(
    12, 12, vague, vague,
    criterion_posterior(c(10, 1), c(5, 2), gamma_go = 0.5, gamma_nogo = 0.3,
                        direction = c("greater", "less"))
  )
  # Before the sentence that introduces the table of the regions
  expect_match(report_criteria(one)[[2L]][["text"]],
               paste("of the study. Smaller effects are favorable for",
                     "endpoint 2. The target values"), fixed = TRUE)
  both <- design_continuous2(
    12, 12, vague, vague,
    criterion_posterior(c(-10, 1), c(-5, 2), gamma_go = 0.5,
                        gamma_nogo = 0.3, direction = "less")
  )
  expect_match(report_criteria(both)[[2L]][["text"]],
               paste("of the study. Smaller effects are favorable for both",
                     "endpoints. The target values"), fixed = TRUE)
})
