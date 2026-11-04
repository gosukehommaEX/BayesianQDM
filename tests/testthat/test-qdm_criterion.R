test_that("qdm_criterion() sets the decision criterion", {
  cr <- crit_post()
  spec <- qdm_endpoint("binary") |> qdm_criterion(cr)
  expect_identical(spec$args[["criterion"]], cr)
  cr2 <- crit_post(0.6, 0.7)
  spec <- spec |> qdm_criterion(cr2)
  expect_identical(spec$args[["criterion"]], cr2)
  expect_length(spec$args, 1L)
})

test_that("qdm_criterion() needs a design specification", {
  expect_error(qdm_criterion(crit_post(), crit_post()),
               "'x' must be a design specification", fixed = TRUE)
})
