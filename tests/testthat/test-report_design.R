test_that("report_design() gives the four sections of a design in order", {
  des <- ref_design("D1")
  res <- report_design(des)
  expect_identical(res, c(report_study(des), report_model(des),
                          report_criteria(des), report_computation(des)))
  headings <- vapply(Filter(function(b) b[["type"]] == "heading", res),
                     function(b) b[["text"]], character(1))
  expect_identical(headings, c("Study design", "Statistical model",
                               "Decision criteria",
                               "Computation of the probabilities"))
})
