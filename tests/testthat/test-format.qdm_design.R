test_that("format.qdm_design() describes a design with concurrent control", {
  expect_identical(
    format(ref_design("D2")),
    c("Design: one binary endpoint",
      "  Treatment: n = 10, prior Beta(a = 1, b = 1)",
      "    External data: binary, n = 10, y = 6, a0 = 0.5",
      "  Control: concurrent, n = 8, prior Beta(a = 2, b = 3)",
      "    External data: binary, n = 20, y = 8, a0 = 0.5",
      "Posterior probability criterion",
      "  Go:   P(theta > 0.2 | data) >= 0.5",
      "  NoGo: P(theta <= 0.1 | data) >= 0.8")
  )
})

test_that("format.qdm_design() describes the other control forms", {
  f3 <- format(ref_design("D3"))
  expect_identical(f3[3:4],
                   c("  Control: hypothetical, prior Beta(a = 1, b = 1)",
                     "    Hypothetical control: data (binary, n = 12, y = 4)"))
  f4 <- format(ref_design("D4"))
  expect_identical(f4[3:4],
                   c("  Control: external data only, prior Beta(a = 1, b = 1)",
                     "    External data: binary, n = 20, y = 6, a0 = 0.5"))
  f5 <- format(ref_design("D5"))
  expect_identical(f5[4], paste("Posterior predictive probability criterion",
                                "(m_t = 30, m_c = 20)"))
})

test_that("format.qdm_design() rejects unknown endpoints", {
  des <- ref_design("D1")
  des$endpoint <- "unknown"
  expect_error(format(des), "Unknown design.")
})

test_that("format.qdm_design() describes continuous designs", {
  expect_identical(
    format(ref_design_cont("C2")),
    c("Design: one continuous endpoint",
      paste("  Treatment: n = 10, prior N-Inv-Chisq(mu0 = 5, kappa0 = 2,",
            "nu0 = 4, sigma0 = 20)"),
      "    External data: continuous, n = 6, mean = 20, sd = 30, a0 = 0.5",
      paste("  Control: concurrent, n = 8, prior N-Inv-Chisq(mu0 = -5,",
            "kappa0 = 1, nu0 = 3, sigma0 = 25)"),
      "    External data: continuous, n = 30, mean = -12, sd = 20, a0 = 0.4",
      "Posterior probability criterion",
      "  Go:   P(theta > 30 | data) >= 0.6",
      "  NoGo: P(theta <= 15 | data) >= 0.7",
      "Probabilities: numerical integration")
  )
  f3 <- format(ref_design_cont("C3", "MC", 500L))
  expect_identical(f3[2:4],
                   c("  Treatment: n = 12, prior Vague",
                     "  Control: hypothetical",
                     "    Hypothetical control: mean = (-10), r = 1.5"))
  expect_identical(f3[length(f3)],
                   "Probabilities: Monte Carlo simulation (nMC = 500)")
  f4 <- format(ref_design_cont("C4", "MM"))
  expect_identical(f4[3], "  Control: external data only, prior Vague")
  expect_identical(f4[length(f4)],
                   "Probabilities: moment-matching approximation")
})
