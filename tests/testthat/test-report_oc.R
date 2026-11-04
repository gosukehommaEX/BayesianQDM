test_that("report_oc() gives the method, the table and the legend", {
  oc <- qdm_oc(ref_design("D1"), scen_c)
  expect_identical(attr(oc, "method"), "exact enumeration of the outcomes")
  expect_identical(
    report_oc(oc, 3L),
    list(list(type = "heading", text = "Operating characteristics"),
         list(type = "par",
              text = paste("The probabilities of Go, NoGo, Gray and Miss",
                           "decisions under scenarios for the true parameters",
                           "(exact enumeration of the outcomes) are as",
                           "follows.")),
         report_scenarios(oc, 3L),
         list(type = "par", text = report_legend(ref_design("D1"))))
  )
})

test_that("report_oc() works without the method", {
  oc <- qdm_oc(ref_design("D1"), scen_c)
  attr(oc, "method") <- NULL
  expect_identical(report_oc(oc, 3L)[[2L]][["text"]],
                   paste("The probabilities of Go, NoGo, Gray and Miss",
                         "decisions under scenarios for the true parameters",
                         "are as follows."))
})
