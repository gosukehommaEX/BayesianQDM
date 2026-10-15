test_that("print.qdm_decision() prints the decision and returns invisibly", {
  res <- qdm_decide(ref_design("D3"), data_binary(12, 9))
  expect_output(out <- withVisible(print(res)),
                "Decision: Go\n  Go:   P\\(theta > 0.2 \\| data\\)")
  expect_false(out$visible)
  expect_identical(out$value, res)
})
