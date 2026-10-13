test_that("print.qdm_hypothetical() prints the description invisibly", {
  h <- hypothetical(mean = 1, r = 1)
  expect_output(print(h), "Hypothetical control: mean = \\(1\\), r = 1")
  expect_invisible(print(h))
})
