test_that("qdm_computation() sets the computation of the probabilities", {
  spec <- qdm_endpoint("continuous") |>
    qdm_computation(method = "MC", nMC = 500L)
  expect_identical(spec$args[["method"]], "MC")
  expect_identical(spec$args[["nMC"]], 500L)
  spec2 <- qdm_endpoint("binary2") |>
    qdm_computation(nMC = 1000L, enumerate = FALSE)
  expect_identical(spec2$args[["nMC"]], 1000L)
  expect_identical(spec2$args[["enumerate"]], FALSE)
})

test_that("qdm_computation() replaces the arguments that it set before", {
  spec <- qdm_endpoint("continuous") |>
    qdm_computation(method = "MC", nMC = 500L) |>
    qdm_computation(method = "MM")
  expect_identical(names(spec$args), "method")
  expect_identical(spec$args[["method"]], "MM")
})

test_that("qdm_computation() rejects arguments that the endpoint does not use", {
  expect_error(qdm_endpoint("binary") |> qdm_computation(method = "MC"),
               paste("qdm_computation() cannot set 'method', which",
                     "design_binary() does not have."),
               fixed = TRUE)
  expect_error(qdm_endpoint("tte") |>
                 qdm_computation(method = "MC", enumerate = TRUE),
               paste("qdm_computation() cannot set 'method', 'enumerate',",
                     "which design_tte() does not have."),
               fixed = TRUE)
  expect_error(qdm_endpoint("continuous2") |>
                 qdm_computation(enumerate = TRUE),
               "which design_continuous2() does not have.", fixed = TRUE)
})

test_that("qdm_computation() changes the computation of a design", {
  expect_identical(
    ref_design_cont("C1") |> qdm_computation(method = "MM") |> qdm_build(),
    ref_design_cont("C1", method = "MM")
  )
})
