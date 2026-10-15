test_that("print.qdm_design() prints the description and returns invisibly", {
  des <- ref_design("D3")
  expect_output(res <- withVisible(print(des)),
                "Design: one binary endpoint\n  Treatment: n = 12")
  expect_false(res$visible)
  expect_identical(res$value, des)
})
