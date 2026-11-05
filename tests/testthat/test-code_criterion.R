test_that("code_criterion() writes the criteria without their defaults", {
  flat <- function(x) code_flat(code_criterion(x))
  expect_identical(flat(crit_post()),
                   paste("criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,",
                         "gamma_go = 0.5, gamma_nogo = 0.8)"))
  expect_identical(flat(ref_design("D5")$criterion),
                   paste("criterion_predictive(theta_null = 0.05, m_t = 30,",
                         "m_c = 20, gamma_go = 0.9, gamma_nogo = 0.7)"))
  # m_c equal to m_t is left out
  expect_identical(flat(ref_design("D9")$criterion),
                   paste("criterion_predictive(theta_null = 0, m_t = 10,",
                         "gamma_go = 0.6, gamma_nogo = 0.6)"))
  expect_identical(flat(ref_design_bin2("B3")$criterion),
                   paste("criterion_posterior(theta_tv = c(0.15, 0.15),",
                         "theta_mav = c(0.1, 0.1), gamma_go = 0.7,",
                         "gamma_nogo = 0.7, go_regions = 1:3,",
                         "nogo_regions = 7:9)"))
  expect_identical(flat(ref_design_bin2("B4")$criterion),
                   paste("criterion_posterior(theta_tv = c(0.15, 0.15),",
                         "theta_mav = c(0.1, 0.1), gamma_go = c(0.5, 0.4),",
                         "gamma_nogo = c(0.6, 0.5), rule = \"and\")"))
  expect_identical(flat(crit_hr()),
                   paste("criterion_posterior(theta_tv = 0.7, theta_mav = 1,",
                         "gamma_go = 0.6, gamma_nogo = 0.7, direction =",
                         "\"less\")"))
})

test_that("code_criterion() gives code that recreates the criteria", {
  crits <- list(
    crit_post(), crit_hr(), ref_design("D5")$criterion,
    ref_design_bin2("B2")$criterion, ref_design_bin2("B6")$criterion,
    criterion_posterior(c(0.2, -0.1), c(0.1, 0), 0.5, 0.5,
                        go_regions = c(1, 3), nogo_regions = c(6, 9),
                        direction = c("greater", "less")),
    criterion_posterior(c(0.2, 0.2), c(0.1, 0.1), 0.6, 0.6, rule = "or"),
    criterion_predictive(c(0.1, 0), m_t = 4, m_c = 3, gamma_go = 0.4,
                         gamma_nogo = 0.3, nogo_regions = c(3, 4),
                         direction = "less")
  )
  for (k in seq_along(crits)) {
    code <- code_flat(code_criterion(crits[[k]]))
    expect_identical(eval(parse(text = code)), crits[[k]], label = code)
  }
  expect_identical(
    code_flat(code_criterion(crits[[6]])),
    paste("criterion_posterior(theta_tv = c(0.2, -0.1), theta_mav = c(0.1, 0),",
          "gamma_go = 0.5, gamma_nogo = 0.5, go_regions = c(1, 3),",
          "nogo_regions = c(6, 9), direction = c(\"greater\", \"less\"))")
  )
  expect_match(code_flat(code_criterion(crits[[7]])),
               "gamma_go = 0.6, gamma_nogo = 0.6, rule = \"or\")", fixed = TRUE)
})
