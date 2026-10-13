test_that("print.qdm_oc() prints the method and the table", {
  oc <- qdm_oc(ref_design("D1"), scen_c[2:3, ])
  expect_output(res <- withVisible(print(oc)),
                paste0("Operating characteristics \\(exact enumeration of ",
                       "the outcomes\\)\n pi_t pi_c"))
  expect_false(res$visible)
  expect_identical(res$value, oc)
  expect_output(print(oc, digits = 2), "0.15")
})

test_that("print.qdm_oc() works after the attributes are dropped", {
  oc <- qdm_oc(ref_design("D1"), scen_c)
  attr(oc, "method") <- NULL
  out <- capture.output(print(oc))
  expect_match(out[1], "pi_t")
})
