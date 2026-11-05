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

test_that("qdm_sample_size() changes the sample sizes of a design", {
  des <- ref_design("D1")
  expect_identical(des |> qdm_sample_size(n_t = 20, n_c = 16) |> qdm_build(),
                   design_binary(20, 16, prior_beta(1, 1), prior_beta(1, 1),
                                 crit_post()))
  expect_identical(qdm_sample_size(des, 20, 16),
                   qdm_sample_size(design_spec(des), 20, 16))
  expect_error(qdm_sample_size(list(), 12, 12),
               paste("'x' must be a design specification started by",
                     "qdm_endpoint(), or a design"), fixed = TRUE)
})
