test_that("check_design() accepts designs only", {
  des <- ref_design("D1")
  expect_identical(check_design(des), des)
  expect_invisible(check_design(des))
  expect_error(check_design(crit_post()),
               "'design' must be created by a design_\\*\\(\\) function")
})

test_that("check_design() builds a design specification", {
  expect_identical(check_design(spec_d1()), ref_design("D1"))
  expect_error(check_design(qdm_endpoint("binary")),
               "The design specification is incomplete", fixed = TRUE)
  expect_error(check_design(crit_post()), "or by qdm_endpoint().",
               fixed = TRUE)
})
