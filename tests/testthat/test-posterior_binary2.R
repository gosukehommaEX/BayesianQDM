test_that("posterior_binary2() reproduces Yamaguchi et al. (submitted), S1.1.2", {
  # Power priors of case study 3: Dir(1.65, 0.95, 0.95, 1.65) for the
  # treatment group (a0t = 0.7, counts (2, 1, 1, 2)) and
  # Dir(7.25, 1.25, 1.75, 0.75) for the control group (a0c = 0.5, counts
  # (14, 2, 3, 1)), with Dir(0.25, 0.25, 0.25, 0.25) initial priors
  des <- ref_design_bin2("B6")
  expect_equal(posterior_binary2(des$prior_t, c(0, 0, 0, 0), des$external_t),
               matrix(c(1.65, 0.95, 0.95, 1.65), 1L,
                      dimnames = list(NULL, c("00", "01", "10", "11"))))
  expect_equal(posterior_binary2(des$prior_c, c(0, 0, 0, 0),
                                 des$external_c)[1L, ],
               c("00" = 7.25, "01" = 1.25, "10" = 1.75, "11" = 0.75))
})

test_that("posterior_binary2() adds the counts row by row", {
  q <- dir_q()
  x <- rbind(c(1, 2, 3, 4), c(0, 0, 5, 0))
  expect_equal(unname(posterior_binary2(q, x)), x + 0.25)
  hyp <- hypothetical(data_binary2(c(2, 1, 1, 0)))
  expect_equal(unname(posterior_binary2(q, c(0, 0, 0, 0),
                                        hypothetical = hyp)),
               matrix(c(2.25, 1.25, 1.25, 0.25), 1L))
})
