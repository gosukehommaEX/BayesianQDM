test_that("oc_sample() dispatches on the endpoint of the design", {
  des <- ref_design("D3")
  sc <- data.frame(pi_t = 0.5)
  expect_identical(oc_sample(des, sc, 10L), oc_sample_binary(des, sc))
})
