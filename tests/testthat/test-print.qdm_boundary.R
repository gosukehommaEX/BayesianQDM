test_that("print.qdm_boundary() prints the tables", {
  b <- qdm_boundary(ref_design("D8"), n_t = c(10, 20))
  expect_output(out <- print(b), "Responders in the treatment group \\(y_t\\)")
  expect_identical(out, b)
  expect_output(print(b, n = 1), "1 more value of n_t")
  expect_output(print(b, digits = 2), "Observed difference")
  expect_error(print(b, n = 0), "'n' must be at least 1.")
  expect_error(print(b, digits = 1.5),
               "'digits' must contain whole numbers only.")
})

test_that("print.qdm_boundary() prints a subset as a data frame", {
  b <- qdm_boundary(ref_design("D8"))
  sub <- b[1:2, ]
  attr(sub, "axes") <- NULL
  expect_output(print(sub), "decision")
  expect_output(print(sub), "NoGo")
})
