test_that("print.qdm_boundary() prints the counts and the first rows", {
  b <- qdm_boundary(ref_design("D8"))
  counts <- table(b$decision)
  header <- sprintf(paste("Decision boundary on a grid of 21 points",
                          "(Go: %d, NoGo: %d, Gray: %d, Miss: 0)"),
                    counts[["Go"]], counts[["NoGo"]], counts[["Gray"]])
  out <- capture.output(res <- withVisible(print(b)))
  expect_identical(out[1L], header)
  expect_match(out[2L], "n_t y_t")
  expect_length(out, 1L + 7L + 1L)
  expect_identical(out[9L], "... 15 more rows")
  expect_false(res$visible)
  expect_identical(res$value, b)
})

test_that("print.qdm_boundary() prints small grids in full", {
  b <- qdm_boundary(ref_design_tte("T1"), events_t = 6, exposure_t = 150)
  out <- capture.output(print(b))
  expect_identical(out[1L], paste("Decision boundary on a grid of 1 point",
                                  "(Go: 1, NoGo: 0, Gray: 0, Miss: 0)"))
  expect_length(out, 3L)
  expect_false(any(grepl("more rows", out)))
})

test_that("print.qdm_boundary() counts a single remaining row", {
  out <- capture.output(print(qdm_boundary(ref_design("D8"), n_t = 6)))
  expect_identical(out[length(out)], "... 1 more row")
})
