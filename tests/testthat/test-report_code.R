test_that("report_code() gives the R code of the design", {
  des <- ref_design("D1")
  blocks <- report_code(des)
  expect_identical(blocks[[1]], report_block("heading", text = "R code"))
  expect_identical(blocks[[2]]$type, "par")
  expect_identical(blocks[[3]]$type, "pre")
  expect_identical(blocks[[3]]$lines, unclass(show_code(des)))
  # A specification is built, so the code ends with qdm_build()
  expect_identical(report_code(spec_d1()), blocks)
  expect_identical(report_code(qdm_oc(des, scen_c)), blocks)
  expect_identical(blocks[[3]]$lines[length(blocks[[3]]$lines)],
                   "  qdm_build()")
})
