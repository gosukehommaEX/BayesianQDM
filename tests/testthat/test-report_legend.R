test_that("report_legend() explains the columns of the scenarios", {
  expect_identical(report_legend(ref_design("D1")),
                   paste("$\\pi_t$ and $\\pi_c$ are the true response rates",
                         "of the treatment and control groups."))
  expect_identical(report_legend(ref_design("D8")),
                   "$\\pi_t$ is the true response rate of the treatment group.")
  expect_identical(report_legend(ref_design_cont("C1")),
                   paste("$\\mu_k$ and $\\sigma_k$ are the true mean and",
                         "standard deviation of group $k$."))
  expect_match(report_legend(ref_design_bin2("B4")),
               "$\\rho_k$ is the correlation between the responses",
               fixed = TRUE)
  expect_match(report_legend(ref_design_cont2("K1")),
               "$\\rho_k$ is the correlation between the two endpoints",
               fixed = TRUE)
  expect_identical(report_legend(ref_design_tte("T2")),
                   paste("$M_t$ and $M_c$ are the true median survival times",
                         "of the treatment and control groups."))
  expect_identical(report_legend(ref_design_tte("T1")),
                   paste("$M_t$ is the true median survival time of the",
                         "treatment group."))
})
