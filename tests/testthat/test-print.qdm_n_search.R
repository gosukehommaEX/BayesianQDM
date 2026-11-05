test_that("print.qdm_n_search() prints the search and the probabilities", {
  res <- find_n_d1(n_t = seq(6, 20, by = 2))
  out <- capture.output(v <- withVisible(print(res)))
  expect_false(v$visible)
  expect_identical(v$value, res)
  expect_identical(out, c(
    format(res),
    "Probabilities for the candidate sample sizes",
    capture.output(print(res$curve, digits = 4L, row.names = FALSE)),
    "Decision probabilities with the chosen sample size",
    capture.output(print(res$oc, digits = 4L, row.names = FALSE))
  ))
  none <- suppressWarnings(find_n_d1(n_t = c(6, 8)))
  expect_identical(capture.output(print(none)), c(
    format(none),
    "Probabilities for the candidate sample sizes",
    capture.output(print(none$curve, digits = 4L, row.names = FALSE))
  ))
})
