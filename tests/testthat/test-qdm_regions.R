test_that("qdm_regions() follows Table 1 of Yamaguchi et al.", {
  r <- qdm_regions("posterior")
  expect_s3_class(r, "noquote")
  expect_identical(dim(r), c(3L, 3L))
  # R1 to R3 run down the first column (theta1 > TV1)
  expect_identical(unclass(r)[, 1], c("theta2 > TV2" = "R1",
                                      "MAV2 < theta2 <= TV2" = "R2",
                                      "theta2 <= MAV2" = "R3"))
  expect_identical(unclass(r)["theta2 <= MAV2", "theta1 <= MAV1"], "R9")
})

test_that("qdm_regions() follows Table 2 of Yamaguchi et al.", {
  r <- unclass(qdm_regions("predictive"))
  expect_identical(dim(r), c(2L, 2L))
  expect_identical(r["theta2 <= NULL2", "theta1 > NULL1"], "R2")
  expect_identical(r["theta2 > NULL2", "theta1 <= NULL1"], "R3")
  expect_error(qdm_regions("other"), "'arg' should be one of")
})

test_that("qdm_regions() shows the posterior regions by default", {
  expect_identical(qdm_regions(), qdm_regions("posterior"))
  expect_identical(sort(as.vector(unclass(qdm_regions()))), paste0("R", 1:9))
  expect_identical(sort(as.vector(unclass(qdm_regions("predictive")))),
                   paste0("R", 1:4))
})
