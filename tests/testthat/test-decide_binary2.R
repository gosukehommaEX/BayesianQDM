test_that("decide_binary2() returns the region probabilities and posteriors", {
  # Python_reference/ref_binary2_qdm.py, case A1
  des <- ref_design_bin2("B1")
  res <- decide_binary2(des, data_binary2(c(2, 1, 3, 4)),
                        data_binary2(c(4, 2, 3, 1)))
  ref <- c(0.421763759173932, 0.221546787701068, 0.225275968676342,
           0.131413484448658)
  expect_equal(res$prob_regions, c(R1 = ref[1], R2 = ref[2], R3 = ref[3],
                                   R4 = ref[4]), tolerance = 1e-12)
  expect_equal(res$prob_go, ref[1], tolerance = 1e-12)
  expect_equal(res$prob_nogo, ref[4], tolerance = 1e-12)
  expect_equal(res$posterior_t,
               c("00" = 2.25, "01" = 1.25, "10" = 3.25, "11" = 4.25))
  # Marginal rule: one probability per endpoint
  res2 <- decide_binary2(ref_design_bin2("B2"), data_binary2(c(2, 1, 3, 4)),
                         data_binary2(c(4, 2, 3, 1)))
  expect_equal(res2$prob_go, c(ref[1] + ref[2], ref[1] + ref[3]),
               tolerance = 1e-12)
  expect_equal(res2$prob_nogo, c(ref[3] + ref[4], ref[2] + ref[4]),
               tolerance = 1e-12)
})

test_that("decide_binary2() uses the fixed control distribution", {
  res <- decide_binary2(ref_design_bin2("B7"), data_binary2(c(1, 0, 1, 1)),
                        NULL)
  expect_equal(res$posterior_c,
               c("00" = 2.25, "01" = 1.25, "10" = 1.25, "11" = 0.25))
  res8 <- decide_binary2(ref_design_bin2("B8"), data_binary2(c(1, 0, 1, 1)),
                         NULL)
  expect_equal(res8$posterior_c,
               c("00" = 7.25, "01" = 1.25, "10" = 1.75, "11" = 0.75))
})

test_that("decide_binary2() validates the data", {
  des <- ref_design_bin2("B1")
  expect_error(decide_binary2(des, data_binary(3, 1), NULL),
               "'data_t' must be created by data_binary2\\(\\).")
  expect_error(decide_binary2(des, data_binary2(c(1, 1, 1, 0)), NULL),
               "'data_c' must be created by data_binary2\\(\\) for a")
  expect_error(decide_binary2(ref_design_bin2("B7"),
                              data_binary2(c(1, 1, 1, 0)),
                              data_binary2(c(1, 1, 1, 0))),
               "'data_c' must be NULL for a design without a concurrent")
})
