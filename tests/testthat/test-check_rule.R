test_that("check_rule() returns the settings of the criteria", {
  one <- check_rule("and", 1L, 0.5, 0.6, 1L, 9L, FALSE, 9L)
  expect_identical(one, list(rule = "joint", gamma_go = 0.5,
                             gamma_nogo = 0.6, go = NULL, nogo = NULL))
  joint <- check_rule("joint", 2L, 0.5, 0.6, c(2, 1), 9, TRUE, 9L)
  expect_identical(joint$go, 1:2)
  expect_identical(joint$nogo, 9L)
  marg <- check_rule("or", 2L, 0.5, c(0.6, 0.7), 1L, 9L, FALSE, 9L)
  expect_identical(marg, list(rule = "or", gamma_go = c(0.5, 0.5),
                              gamma_nogo = c(0.6, 0.7), go = NULL,
                              nogo = NULL))
  expect_error(check_rule("or", 2L, 0.5, 0.5, 1L, 9L, TRUE, 9L),
               "'go_regions' and 'nogo_regions' are used only")
  expect_error(check_rule("and", 2L, 0.5, c(0.1, 0.2, 0.3), 1L, 9L, FALSE,
                          9L),
               "With rule = \"and\", 'gamma_go' and 'gamma_nogo' must")
})
