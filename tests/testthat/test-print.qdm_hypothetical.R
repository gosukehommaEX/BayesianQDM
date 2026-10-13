test_that("print.qdm_hypothetical() prints the description invisibly", {
  h <- hypothetical(mean = 1, r = 1)
  expect_output(res <- withVisible(print(h)),
                "Hypothetical control: mean = \\(1\\), r = 1")
  expect_false(res$visible)
  expect_identical(res$value, h)
})
