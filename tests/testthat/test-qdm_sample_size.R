test_that("qdm_sample_size() sets n_t and n_c", {
  spec <- qdm_endpoint("binary") |> qdm_sample_size(n_t = 12, n_c = 10)
  expect_s3_class(spec, "qdm_spec", exact = TRUE)
  expect_identical(spec$args[["n_t"]], 12)
  expect_identical(spec$args[["n_c"]], 10)
})

test_that("qdm_sample_size() replaces the sample sizes set before", {
  spec <- qdm_endpoint("binary") |>
    qdm_sample_size(n_t = 12, n_c = 10) |>
    qdm_sample_size(n_t = 20, n_c = 0)
  expect_identical(spec$args[["n_t"]], 20)
  expect_identical(spec$args[["n_c"]], 0)
  expect_length(spec$args, 2L)
})

test_that("qdm_sample_size() needs a design specification", {
  expect_error(qdm_sample_size(ref_design("D1"), 12, 12),
               "'x' must be a design specification started by qdm_endpoint().",
               fixed = TRUE)
})
