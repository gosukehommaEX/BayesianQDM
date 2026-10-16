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

test_that("format.qdm_design() describes designs with two binary endpoints", {
  f6 <- format(ref_design_bin2("B6"))
  expect_identical(f6[1:3],
                   c("Design: two binary endpoints",
                     paste("  Treatment: n = 8, prior Dirichlet(alpha =",
                           "(0.25, 0.25, 0.25, 0.25))"),
                     paste("    External data: two binary endpoints, n = 6,",
                           "counts (00, 01, 10, 11) = (2, 1, 1, 2), a0 = 0.7")))
  expect_identical(f6[length(f6) - 1:0],
                   c("Probabilities: Monte Carlo simulation (nMC = 20000)",
                     "Operating characteristics: enumeration of all outcomes"))
  f1 <- format(ref_design_bin2("B1", enumerate = FALSE))
  expect_identical(f1[length(f1) - 1:0],
                   c("Probabilities: exact computation",
                     "Operating characteristics: simulation"))
})

test_that("format.qdm_design() describes designs with two continuous endpoints", {
  f2 <- format(ref_design_cont2("K2"))
  expect_identical(f2[1:2],
                   c("Design: two continuous endpoints",
                     paste("  Treatment: n = 30, prior N-Inv-Wishart(mu0 =",
                           "(0, 1), kappa0 = 0.1, nu0 = 4, Lambda0 = (6400,",
                           "15; 15, 36))")))
  expect_identical(f2[length(f2)],
                   "Probabilities: moment-matching approximation")
  f3 <- format(ref_design_cont2("K3"))
  expect_identical(f3[3], paste("    External data: two continuous endpoints,",
                                "n = 8, mean = (2, 1), cov = (25, 3; 3, 9),",
                                "a0 = 0.5"))
  f4 <- format(ref_design_cont2("K4", "MC", 2000L))
  expect_identical(f4[3:4],
                   c("  Control: hypothetical",
                     "    Hypothetical control: mean = (-10, -1), r = 1.5"))
  expect_identical(f4[length(f4)],
                   "Probabilities: Monte Carlo simulation (nMC = 2000)")
})

test_that("format.qdm_design() describes a fixed control", {
  f8 <- format(ref_design("D8"))
  expect_identical(f8[3], "  Control: fixed value 0.2")
  f7 <- format(ref_design_cont("C7"))
  expect_identical(f7[3], "  Control: fixed value -10")
  expect_identical(f7[length(f7)],
                   "Probabilities: exact computation (fixed control)")
  fk <- format(ref_design_cont2("KF1"))
  expect_identical(fk[3], "  Control: fixed value (-10, -1)")
  fb <- format(ref_design_bin2("BF1"))
  expect_identical(fb[3], "  Control: fixed value (0.25, 0.50)")
})

test_that("format.qdm_design() describes time-to-event designs", {
  f1 <- format(ref_design_tte("T1"))
  expect_identical(
    f1,
    c("Design: one time-to-event endpoint (milestone survival at t = 3)",
      "  Treatment: n = 30, prior Gamma(shape = 0.1, rate = 1)",
      "  Control: fixed value 0.72",
      "Posterior probability criterion",
      "  Go:   P(theta > 0.15 | data) >= 0.6",
      "  NoGo: P(theta <= 0 | data) >= 0.8",
      "Probabilities: exact computation (fixed control)",
      "Conduct of the study (operating characteristics):",
      "  Enrollment: uniform from 0 to 1",
      "  Dropout: none",
      "  Analysis: at time 1e+06")
  )
  f4 <- format(ref_design_tte("T4", nMC = 500L))
  expect_identical(f4[length(f4) - 3L], "Conduct of the future study:")
  expect_identical(f4[7], "Probabilities: simulation of future studies (nMC = 500)")
  f6 <- format(ref_design_tte("T6"))
  expect_identical(f6[3], "    External data: time-to-event, events = 4, exposure = 30, a0 = 0.5")
  des <- ref_design_tte("T2")
  des$conduct <- NULL
  expect_false(any(grepl("Conduct", format(des))))
})

test_that("format.qdm_design() describes designs with the hazard ratio", {
  f2 <- format(ref_design_hr("H2"))
  expect_identical(
    f2[1:7],
    c("Design: one time-to-event endpoint (hazard ratio)",
      "  Treatment: n = 30, prior Gamma(shape = 0.1, rate = 1)",
      "  Control: fixed median survival time 6",
      "Posterior probability criterion",
      "  Go:   P(theta < 0.7 | data) >= 0.6",
      "  NoGo: P(theta >= 1 | data) >= 0.7",
      "Probabilities: exact computation (fixed control)")
  )
  expect_identical(format(ref_design_hr("H1"))[7],
                   "Probabilities: exact computation (F distribution)")
})
