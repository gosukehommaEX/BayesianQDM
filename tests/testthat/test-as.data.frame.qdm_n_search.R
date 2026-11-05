test_that("as.data.frame() gives the chosen sample sizes", {
  res <- find_n_d1(n_t = seq(6, 20, by = 2))
  expect_identical(as.data.frame(res),
                   data.frame(n_t = 14L, n_c = 14L, prob_go = res$prob_go,
                              target_go = 0.1, prob_nogo = res$prob_nogo,
                              target_nogo = 0.1, n_t_stable = 20L,
                              n_c_stable = 20L))
  expect_equal(as.data.frame(res)$prob_go, ref_find_n$prob_go[5L],
               tolerance = 1e-10)
  expect_identical(rownames(as.data.frame(res, row.names = "D1")), "D1")
  none <- suppressWarnings(find_n_d1(n_t = c(6, 8)))
  expect_true(is.na(as.data.frame(none)$n_t))
})
