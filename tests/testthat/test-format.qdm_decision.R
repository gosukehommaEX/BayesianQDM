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
