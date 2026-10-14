test_that("format.qdm_decision() reports probabilities and criteria", {
  res <- qdm_decide(ref_design("D1"), data_binary(12, 9), data_binary(12, 4))
  expect_identical(
    format(res),
    c("Decision: Go",
      "  Go:   P(theta > 0.2 | data) = 0.8207 (threshold 0.5, met)",
      "  NoGo: P(theta <= 0.1 | data) = 0.07293 (threshold 0.8, not met)",
      "Posterior distributions",
      "  Treatment: Beta(a = 10, b = 4)",
      "  Control:   Beta(a = 5, b = 9)")
  )
})

test_that("format.qdm_decision() describes predictive criteria", {
  res <- qdm_decide(ref_design("D5"), data_binary(10, 4), data_binary(10, 4))
  out <- format(res, digits = 3)
  expect_identical(out[1], "Decision: Gray")
  expect_identical(out[2], paste("  Go:   P(future theta > 0.05 | data) =",
                                 "0.406 (threshold 0.9, not met)"))
})

test_that("format.qdm_decision() describes continuous decisions", {
  d <- ref_data_cont$C1
  res <- qdm_decide(ref_design_cont("C1"), d$t, d$c)
  expect_identical(
    format(res),
    c("Decision: Gray",
      "  Go:   P(theta > 30 | data) = 0.5 (threshold 0.6, not met)",
      "  NoGo: P(theta <= 15 | data) = 0.08784 (threshold 0.7, not met)",
      "Posterior distributions",
      "  Treatment: mean ~ t(df = 11, location = 18, scale = 7.217)",
      "  Control:   mean ~ t(df = 11, location = -12, scale = 7.217)")
  )
})

test_that("format.qdm_decision() describes two binary endpoints", {
  d_t <- data_binary2(c(2, 1, 3, 4))
  d_c <- data_binary2(c(4, 2, 3, 1))
  res <- qdm_decide(ref_design_bin2("B1"), d_t, d_c)
  out <- format(res, digits = 3)
  expect_identical(out[1:3],
                   c("Decision: Go",
                     "  Go:   P(regions {1} | data) = 0.422 (threshold 0.4, met)",
                     paste("  NoGo: P(regions {4} | data) = 0.131",
                           "(threshold 0.3, not met)")))
  expect_identical(out[4], paste("Region probabilities: R1 = 0.422,",
                                 "R2 = 0.222, R3 = 0.225, R4 = 0.131"))
  expect_identical(out[6], paste("  Treatment: Dirichlet(alpha = (2.25,",
                                 "1.25, 3.25, 4.25))"))
  res2 <- qdm_decide(ref_design_bin2("B2"), d_t, d_c)
  out2 <- format(res2, digits = 3)
  expect_identical(out2[2:4],
                   c("  Go:   met (rule \"and\")",
                     paste("    P(future theta1 > 0.1 | data) = 0.643",
                           "(threshold 0.6, met)"),
                     paste("    P(future theta2 > 0 | data) = 0.647",
                           "(threshold 0.5, met)")))
  expect_identical(out2[5], "  NoGo: not met (rule \"and\")")
})
