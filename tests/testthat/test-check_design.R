test_that("check_design() accepts designs only", {
  des <- ref_design("D1")
  expect_identical(check_design(des), des)
  expect_invisible(check_design(des))
  expect_error(check_design(crit_post()),
               "'design' must be created by a design_\\*\\(\\) function")
})
