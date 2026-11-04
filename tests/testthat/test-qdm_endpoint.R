test_that("qdm_endpoint() starts an empty specification for each endpoint", {
  for (ep in c("binary", "continuous", "binary2", "continuous2", "tte")) {
    spec <- qdm_endpoint(ep)
    expect_s3_class(spec, "qdm_spec", exact = TRUE)
    expect_identical(spec$endpoint, ep)
    expect_length(spec$args, 0L)
  }
})

test_that("qdm_endpoint() rejects an unknown endpoint", {
  msg <- paste0("'endpoint' must be one of \"binary\", \"continuous\", ",
                "\"binary2\", \"continuous2\", \"tte\".")
  expect_error(qdm_endpoint(), msg, fixed = TRUE)
  expect_error(qdm_endpoint("survival"), msg, fixed = TRUE)
  expect_error(qdm_endpoint(c("binary", "tte")), msg, fixed = TRUE)
  expect_error(qdm_endpoint(1), msg, fixed = TRUE)
  expect_error(qdm_endpoint(NA_character_), msg, fixed = TRUE)
})
