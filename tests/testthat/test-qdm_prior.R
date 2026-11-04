test_that("qdm_prior() sets the prior distributions", {
  pt <- prior_beta(0.5, 0.5)
  pc <- prior_beta(1, 1)
  spec <- qdm_endpoint("binary") |> qdm_prior(prior_t = pt, prior_c = pc)
  expect_identical(spec$args[["prior_t"]], pt)
  expect_identical(spec$args[["prior_c"]], pc)
})

test_that("qdm_prior() without prior_c removes the control prior set before", {
  spec <- qdm_endpoint("binary") |>
    qdm_prior(prior_beta(1, 1), prior_beta(1, 1)) |>
    qdm_prior(prior_beta(0.5, 0.5))
  expect_identical(spec$args[["prior_t"]], prior_beta(0.5, 0.5))
  expect_false("prior_c" %in% names(spec$args))
})

test_that("qdm_prior() needs a design specification", {
  expect_error(qdm_prior(list(), prior_beta(1, 1)),
               "'x' must be a design specification", fixed = TRUE)
})
