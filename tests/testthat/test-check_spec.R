test_that("check_spec() accepts design specifications only", {
  spec <- qdm_endpoint("binary")
  expect_identical(check_spec(spec), spec)
  expect_invisible(check_spec(spec))
  expect_error(check_spec(ref_design("D1")),
               "'x' must be a design specification started by qdm_endpoint().",
               fixed = TRUE)
})
