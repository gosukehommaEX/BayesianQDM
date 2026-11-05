test_that("code_conduct() writes the conduct without its defaults", {
  expect_identical(code_flat(code_conduct(conduct_complete)),
                   "conduct_tte(accrual_time = c(0, 1), cut_time = 1000000)")
  expect_identical(code_flat(code_conduct(conduct_example)),
                   paste("conduct_tte(accrual_time = c(0, 5, 10, 15),",
                         "accrual_rate = c(0.4, 2, 4.5), dropout_hazard =",
                         "0.00226, cut_time = 15.44)"))
  c3 <- conduct_tte(accrual_time = c(0, 2), dropout_time = c(0, 1),
                    dropout_hazard = c(0.01, 0.02), cut_events = 30)
  expect_identical(code_flat(code_conduct(c3)),
                   paste("conduct_tte(accrual_time = c(0, 2), dropout_time =",
                         "c(0, 1), dropout_hazard = c(0.01, 0.02),",
                         "cut_events = 30)"))
  for (x in list(conduct_complete, conduct_example, c3)) {
    expect_identical(eval(parse(text = code_flat(code_conduct(x)))), x)
  }
})
