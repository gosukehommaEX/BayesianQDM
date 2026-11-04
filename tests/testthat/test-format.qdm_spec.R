test_that("format.qdm_spec() lists the arguments set and those still needed", {
  expect_identical(
    format(qdm_endpoint("binary")),
    c("Design specification: one binary endpoint (incomplete)",
      "  Arguments set: none",
      paste("  Still needed: the sample sizes (qdm_sample_size()),",
            "the prior distributions (qdm_prior()), the decision criterion",
            "(qdm_criterion())"))
  )
  spec <- qdm_endpoint("tte") |>
    qdm_criterion(crit_tte()) |>
    qdm_sample_size(30, 0)
  expect_identical(
    format(spec),
    c("Design specification: one time-to-event endpoint (incomplete)",
      "  Arguments set: n_t, n_c, criterion",
      "  Still needed: the prior distributions (qdm_prior())")
  )
})

test_that("format.qdm_spec() names each endpoint", {
  labels <- c(binary = "one binary endpoint",
              continuous = "one continuous endpoint",
              binary2 = "two binary endpoints",
              continuous2 = "two continuous endpoints",
              tte = "one time-to-event endpoint")
  for (ep in names(labels)) {
    expect_identical(format(qdm_endpoint(ep))[1L],
                     sprintf("Design specification: %s (incomplete)",
                             labels[[ep]]))
  }
})

test_that("format.qdm_spec() describes the design of a complete specification", {
  spec <- qdm_endpoint("binary") |>
    qdm_sample_size(12, 12) |>
    qdm_prior(prior_beta(1, 1), prior_beta(1, 1)) |>
    qdm_criterion(crit_post())
  expect_identical(format(spec),
                   c("Design specification (complete):",
                     format(ref_design("D1"))))
  expect_identical(format(spec, digits = 2L),
                   c("Design specification (complete):",
                     format(ref_design("D1"), digits = 2L)))
})

test_that("format.qdm_spec() reports why a complete specification is not valid", {
  spec <- qdm_endpoint("binary") |>
    qdm_sample_size(0, 12) |>
    qdm_prior(prior_beta(1, 1), prior_beta(1, 1)) |>
    qdm_criterion(crit_post())
  expect_identical(format(spec),
                   c("Design specification: one binary endpoint (not valid)",
                     "  'n_t' must be at least 1."))
})
