test_that("qdm_control() sets a hypothetical or a fixed control", {
  hyp <- hypothetical(data_binary(12, 4))
  spec <- qdm_endpoint("binary") |> qdm_control(hypothetical_c = hyp)
  expect_identical(spec$args[["hypothetical_c"]], hyp)
  expect_false("fixed_c" %in% names(spec$args))
  spec <- spec |> qdm_control(fixed_c = 0.2)
  expect_identical(spec$args[["fixed_c"]], 0.2)
  expect_false("hypothetical_c" %in% names(spec$args))
})

test_that("qdm_control() needs a design specification", {
  expect_error(qdm_control(list(), fixed_c = 0.2),
               "'x' must be a design specification", fixed = TRUE)
})
