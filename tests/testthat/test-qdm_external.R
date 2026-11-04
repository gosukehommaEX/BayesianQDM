test_that("qdm_external() sets the external data of both groups", {
  et <- borrow(data_binary(10, 6), 0.5)
  ec <- borrow(data_binary(20, 8), 0.5)
  spec <- qdm_endpoint("binary") |>
    qdm_external(external_t = et, external_c = ec)
  expect_identical(spec$args[["external_t"]], et)
  expect_identical(spec$args[["external_c"]], ec)
})

test_that("qdm_external() removes the external data that it is not given", {
  et <- borrow(data_binary(10, 6), 0.5)
  ec <- borrow(data_binary(20, 8), 0.5)
  spec <- qdm_endpoint("binary") |>
    qdm_external(external_t = et, external_c = ec) |>
    qdm_external(external_c = ec)
  expect_false("external_t" %in% names(spec$args))
  expect_identical(spec$args[["external_c"]], ec)
})

test_that("qdm_external() needs a design specification", {
  expect_error(qdm_external(list(), borrow(data_binary(10, 6), 0.5)),
               "'x' must be a design specification", fixed = TRUE)
})
