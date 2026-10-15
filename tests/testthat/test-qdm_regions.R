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

test_that("qdm_regions() reverses the inequalities for direction = \"less\"", {
  r <- unclass(qdm_regions("posterior", direction = c("greater", "less")))
  expect_identical(rownames(r), c("theta2 < TV2", "TV2 <= theta2 < MAV2",
                                  "theta2 >= MAV2"))
  expect_identical(colnames(r), c("theta1 > TV1", "MAV1 < theta1 <= TV1",
                                  "theta1 <= MAV1"))
  expect_identical(r["theta2 < TV2", "theta1 > TV1"], "R1")
  p <- unclass(qdm_regions("predictive", direction = "less"))
  expect_identical(dimnames(p),
                   list(c("theta2 < NULL2", "theta2 >= NULL2"),
                        c("theta1 < NULL1", "theta1 >= NULL1")))
  expect_identical(p["theta2 >= NULL2", "theta1 < NULL1"], "R2")
  expect_error(qdm_regions(direction = "up"),
               "'direction' must be \"greater\" or \"less\"", fixed = TRUE)
})
