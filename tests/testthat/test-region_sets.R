test_that("region_sets() gives the regions of the criteria", {
  cr <- criterion_posterior(c(0.2, 0.1), c(0, 0), 0.5, 0.5,
                            go_regions = c(1, 2), nogo_regions = 9)
  expect_identical(region_sets(cr), list(go = list(1:2), nogo = list(9L)))
  mp <- criterion_posterior(c(0.2, 0.1), c(0, 0), 0.5, 0.5, rule = "or")
  expect_identical(region_sets(mp),
                   list(go = list(1:3, c(1L, 4L, 7L)),
                        nogo = list(7:9, c(3L, 6L, 9L))))
  mq <- criterion_predictive(c(0, 0), m_t = 10, gamma_go = 0.5,
                             gamma_nogo = 0.5, rule = "and")
  expect_identical(region_sets(mq),
                   list(go = list(1:2, c(1L, 3L)),
                        nogo = list(3:4, c(2L, 4L))))
})
