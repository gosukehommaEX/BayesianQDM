test_that("report_events() writes the probabilities of one endpoint", {
  expect_identical(report_events(crit_post()),
                   list(go = "\\Pr(\\theta > 0.2 \\mid D)",
                        nogo = "\\Pr(\\theta \\le 0.1 \\mid D)"))
  expect_identical(report_events(crit_hr()),
                   list(go = "\\Pr(\\theta < 0.7 \\mid D)",
                        nogo = "\\Pr(\\theta \\ge 1 \\mid D)"))
  expect_identical(
    report_events(criterion_predictive(0.05, m_t = 30, m_c = 20,
                                       gamma_go = 0.9, gamma_nogo = 0.7)),
    list(go = "\\Pr(\\tilde\\theta > 0.05 \\mid D)",
         nogo = "\\Pr(\\tilde\\theta \\le 0.05 \\mid D)")
  )
})

test_that("report_events() writes the probabilities of regions", {
  expect_identical(
    report_events(ref_design_bin2("B6")$criterion),
    list(go = paste("\\Pr\\{(\\theta_1, \\theta_2) \\in R_{1} \\cup R_{2}",
                    "\\mid D\\}"),
         nogo = "\\Pr\\{(\\theta_1, \\theta_2) \\in R_{9} \\mid D\\}")
  )
  expect_identical(
    report_events(ref_design_bin2("B1")$criterion),
    list(go = paste("\\Pr\\{(\\tilde\\theta_1, \\tilde\\theta_2) \\in R_{1}",
                    "\\mid D\\}"),
         nogo = paste("\\Pr\\{(\\tilde\\theta_1, \\tilde\\theta_2) \\in",
                      "R_{4} \\mid D\\}"))
  )
})

test_that("report_events() writes the probabilities of each endpoint", {
  expect_identical(
    report_events(ref_design_bin2("B4")$criterion),
    list(go = c("\\Pr(\\theta_1 > 0.15 \\mid D)",
                "\\Pr(\\theta_2 > 0.15 \\mid D)"),
         nogo = c("\\Pr(\\theta_1 \\le 0.1 \\mid D)",
                  "\\Pr(\\theta_2 \\le 0.1 \\mid D)"))
  )
  crit <- criterion_predictive(c(0.1, 0), m_t = 4, m_c = 3, gamma_go = 0.6,
                               gamma_nogo = 0.5, rule = "and",
                               direction = c("greater", "less"))
  expect_identical(
    report_events(crit),
    list(go = c("\\Pr(\\tilde\\theta_1 > 0.1 \\mid D)",
                "\\Pr(\\tilde\\theta_2 < 0 \\mid D)"),
         nogo = c("\\Pr(\\tilde\\theta_1 \\le 0.1 \\mid D)",
                  "\\Pr(\\tilde\\theta_2 \\ge 0 \\mid D)"))
  )
})
