test_that("qdm_tte() sets the arguments of a time-to-event endpoint", {
  spec <- qdm_endpoint("tte") |>
    qdm_tte(milestone = 3, effect = "milestone", conduct = conduct_example,
            conduct_future = conduct_complete)
  expect_identical(spec$args[["milestone"]], 3)
  expect_identical(spec$args[["effect"]], "milestone")
  expect_identical(spec$args[["conduct"]], conduct_example)
  expect_identical(spec$args[["conduct_future"]], conduct_complete)
})

test_that("qdm_tte() replaces the arguments that it set before", {
  spec <- qdm_endpoint("tte") |>
    qdm_tte(milestone = 3, conduct = conduct_example) |>
    qdm_tte(effect = "hazard_ratio", conduct = conduct_complete)
  expect_setequal(names(spec$args), c("effect", "conduct"))
  expect_identical(spec$args[["effect"]], "hazard_ratio")
  expect_identical(spec$args[["conduct"]], conduct_complete)
})

test_that("qdm_tte() is used only with a time-to-event endpoint", {
  expect_error(qdm_endpoint("binary") |> qdm_tte(milestone = 3),
               "qdm_tte() is used only with qdm_endpoint(\"tte\").",
               fixed = TRUE)
  expect_error(qdm_tte(list(), milestone = 3),
               "'x' must be a design specification", fixed = TRUE)
})

test_that("qdm_tte() changes the time-to-event arguments of a design", {
  expect_identical(
    ref_design_tte("T2") |>
      qdm_tte(milestone = 4, conduct = conduct_example) |>
      qdm_build(),
    design_tte(20, 20, prior_gamma(0.1, 1), prior_gamma(0.1, 1), crit_tte(),
               milestone = 4, conduct = conduct_example)
  )
  expect_error(ref_design("D1") |> qdm_tte(milestone = 3),
               "qdm_tte() is used only with qdm_endpoint(\"tte\").",
               fixed = TRUE)
})
