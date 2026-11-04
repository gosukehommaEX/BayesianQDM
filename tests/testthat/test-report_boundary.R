test_that("report_boundary() gives all the rows of the decision boundaries", {
  b <- qdm_boundary(ref_design("D1"))
  expect_identical(
    report_boundary(b),
    list(list(type = "heading", text = "Decision boundaries"),
         list(type = "par",
              text = paste("The observed data that lead to each decision",
                           "are as follows.")),
         list(type = "pre", lines = format(b, n = nrow(b))))
  )
  # 31 numbers of control responders, more than the 20 that format() shows
  # by default
  unif <- prior_beta(1, 1)
  b2 <- qdm_boundary(design_binary(30, 30, unif, unif, crit_post()))
  expect_true(any(grepl("more values", format(b2), fixed = TRUE)))
  expect_false(any(grepl("more values", report_boundary(b2)[[3L]][["lines"]],
                         fixed = TRUE)))
})

test_that("report_boundary() prints the data frame without the attributes", {
  b <- qdm_boundary(ref_design("D1"))
  attr(b, "axes") <- NULL
  expect_identical(report_boundary(b)[[3L]][["lines"]],
                   utils::capture.output(print(as.data.frame(b))))
})
