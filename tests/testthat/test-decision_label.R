test_that("decision_label() maps the criteria met to decisions", {
  expect_identical(decision_label(c(TRUE, FALSE, FALSE, TRUE),
                                  c(FALSE, TRUE, FALSE, TRUE)),
                   c("Go", "NoGo", "Gray", "Miss"))
  expect_identical(decision_label(TRUE, FALSE), "Go")
})
