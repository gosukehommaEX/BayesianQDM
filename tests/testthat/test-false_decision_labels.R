test_that("false_decision_labels() describes the false Go and NoGo", {
  labs <- false_decision_labels()
  expect_identical(
    labs$legend,
    c("False Go (Go criterion met when Go would be wrong)",
      "False NoGo (NoGo criterion met when NoGo would be wrong)")
  )
  expect_identical(labs$scenario,
                   c("Go would be wrong", "NoGo would be wrong"))
  # The labels do not use the argument names (no underscores)
  expect_false(any(grepl("_", unlist(labs), fixed = TRUE)))
})
