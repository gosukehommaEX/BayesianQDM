test_that("oc_sample() dispatches on the endpoint of the design", {
  des <- ref_design("D3")
  sc <- data.frame(pi_t = 0.5)
  expect_identical(oc_sample(des, sc, 10L), oc_sample_binary(des, sc))
})

test_that("oc_sample() simulates continuous designs", {
  des <- ref_design_cont("C3", "MM")
  set.seed(1)
  a <- oc_sample(des, scen_cont_u, 30L)
  set.seed(1)
  expect_identical(a, oc_sample_continuous(des, scen_cont_u, 30L))
})

test_that("oc_sample() dispatches designs with two binary endpoints", {
  des <- ref_design_bin2("B1", enumerate = FALSE)
  set.seed(1)
  a <- oc_sample(des, scen_bin2, 20L)
  set.seed(1)
  expect_identical(a, oc_sample_binary2(des, scen_bin2, 20L))
})

test_that("oc_sample() dispatches designs with two continuous endpoints", {
  des <- ref_design_cont2("K4")
  set.seed(1)
  a <- oc_sample(des, scen_cont2("K4"), 10L)
  set.seed(1)
  expect_identical(a, oc_sample_continuous2(des, scen_cont2("K4"), 10L))
})
