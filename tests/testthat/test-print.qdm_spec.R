test_that("print.qdm_spec() prints the description and returns invisibly", {
  spec <- qdm_endpoint("binary") |> qdm_sample_size(12, 12)
  expect_output(res <- withVisible(print(spec)),
                "Design specification: one binary endpoint (incomplete)",
                fixed = TRUE)
  expect_false(res$visible)
  expect_identical(res$value, spec)
})

test_that("print.qdm_spec() passes the number of digits to format()", {
  spec <- qdm_endpoint("binary") |>
    qdm_sample_size(12, 12) |>
    qdm_prior(prior_beta(1, 1), prior_beta(1, 1)) |>
    qdm_criterion(crit_post())
  expect_identical(capture.output(print(spec, digits = 2L)),
                   format(spec, digits = 2L))
})
