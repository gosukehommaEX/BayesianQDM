test_that("spec_missing() lists the required arguments that are not set", {
  expect_identical(spec_missing(qdm_endpoint("binary")),
                   c("the sample sizes (qdm_sample_size())",
                     "the prior distributions (qdm_prior())",
                     "the decision criterion (qdm_criterion())"))
  spec <- qdm_endpoint("binary") |> qdm_sample_size(12, NULL)
  expect_identical(spec_missing(spec),
                   c("the sample sizes (qdm_sample_size())",
                     "the prior distributions (qdm_prior())",
                     "the decision criterion (qdm_criterion())"))
  spec <- spec |>
    qdm_sample_size(12, 12) |>
    qdm_criterion(crit_post())
  expect_identical(spec_missing(spec),
                   "the prior distributions (qdm_prior())")
})

test_that("spec_missing() is empty when the required arguments are set", {
  spec <- qdm_endpoint("binary") |>
    qdm_sample_size(12, 0) |>
    qdm_prior(prior_beta(1, 1)) |>
    qdm_criterion(crit_post())
  expect_identical(spec_missing(spec), character(0))
})
