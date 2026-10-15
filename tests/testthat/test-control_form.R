test_that("control_form() identifies the form of the control group", {
  ext <- borrow(data_binary(20, 6), 0.5)
  hyp <- hypothetical(data_binary(12, 4))
  expect_identical(control_form(5, NULL, NULL), "concurrent")
  expect_identical(control_form(5, ext, NULL), "concurrent")
  expect_identical(control_form(0, ext, NULL), "external")
  expect_identical(control_form(0, NULL, hyp), "hypothetical")
  expect_identical(control_form(0, NULL, NULL, 0.2), "fixed")
})

test_that("control_form() rejects inconsistent combinations", {
  ext <- borrow(data_binary(20, 6), 0.5)
  hyp <- hypothetical(data_binary(12, 4))
  expect_error(control_form(5, NULL, hyp), "can be used only when 'n_c' is 0")
  expect_error(control_form(5, NULL, NULL, 0.2),
               "'fixed_c' can be used only when 'n_c' is 0.")
  expect_error(control_form(0, NULL, NULL), "specify 'external_c'")
  expect_error(control_form(0, ext, hyp), "specify only one of")
  expect_error(control_form(0, ext, NULL, 0.2), "specify only one of")
  expect_error(control_form(0, NULL, hyp, 0.2),
               "only one of 'external_c', 'hypothetical_c' and 'fixed_c'.")
})
