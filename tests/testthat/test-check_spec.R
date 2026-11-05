test_that("check_spec() accepts design specifications and designs", {
  spec <- qdm_endpoint("binary")
  expect_identical(check_spec(spec), spec)
  expect_invisible(check_spec(spec))
  # A design is converted to the specification that builds it
  expect_identical(check_spec(ref_design("D1")),
                   design_spec(ref_design("D1")))
  expect_error(check_spec(list()),
               paste("'x' must be a design specification started by",
                     "qdm_endpoint(), or a design created by a design_*()",
                     "function."), fixed = TRUE)
})
