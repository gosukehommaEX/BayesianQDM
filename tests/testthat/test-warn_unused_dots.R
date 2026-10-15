test_that("warn_unused_dots() is silent without arguments", {
  expect_silent(warn_unused_dots(list(), "plot.qdm_oc"))
})

test_that("warn_unused_dots() names the ignored arguments", {
  expect_warning(
    warn_unused_dots(list(point_size = 2), "plot.qdm_boundary"),
    "Ignoring unknown argument(s) of plot.qdm_boundary(): point_size.",
    fixed = TRUE
  )
  expect_warning(
    warn_unused_dots(list(linewidth = 2, 3), "plot.qdm_calibration"),
    paste("Ignoring unknown argument(s) of plot.qdm_calibration():",
          "linewidth, an unnamed argument."),
    fixed = TRUE
  )
  expect_warning(warn_unused_dots(list(3), "plot.qdm_oc"),
                 "plot.qdm_oc(): an unnamed argument.", fixed = TRUE)
})
