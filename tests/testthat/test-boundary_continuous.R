test_that("boundary_continuous() matches the fixed-control references", {
  # Reference: Python_reference/ref_fixed_tte_qdm.py, designs C7 and C8
  b <- boundary_continuous(ref_design_cont("C7"), 12, c(30, 25, 10), 20,
                           NULL)
  expect_identical(b$grid$mean_t, c(10, 25, 30))
  expect_identical(names(b$grid), c("n_t", "mean_t", "sd_t"))
  expect_equal(b$p_go[2L], 0.797515237746, tolerance = 1e-10)
  expect_equal(b$p_nogo[2L], 0.00264736615439, tolerance = 1e-10)
  # With the vague prior, the mean has a t distribution with 11 degrees of
  # freedom, location ybar and scale s / sqrt(n); Go: mean + 10 > 30
  m <- c(10, 25, 30)
  expect_equal(b$p_go, stats::pt((m - 20) / (20 / sqrt(12)), 11))
  expect_equal(b$p_nogo, stats::pt((5 - m) / (20 / sqrt(12)), 11))
  b8 <- boundary_continuous(ref_design_cont("C8"), 15, 9, 12, NULL)
  expect_equal(b8$p_go, 0.845153264581, tolerance = 1e-10)
  expect_equal(b8$p_nogo, 1 - 0.845153264581, tolerance = 1e-10)
})

test_that("boundary_continuous() matches a concurrent-control reference", {
  # Reference: Python_reference/ref_continuous_qdm.py, design C2
  des <- ref_design_cont("C2")
  b <- boundary_continuous(des, 10, c(22, 40), 24, data_continuous(8, -8, 18))
  expect_equal(b$p_go[1L], 0.472868606954, tolerance = 1e-9)
  expect_equal(b$p_nogo[1L], 0.0337620028253, tolerance = 1e-9)
  expect_equal(b$grid$diff, c(30, 48))
  expect_identical(b$axes$y, "diff")
  expect_identical(b$axes$ylab, "Observed difference in means")
  expect_match(b$axes$caption, "^Control data: ")
})

test_that("boundary_continuous() chooses the horizontal axis", {
  des <- ref_design_cont("C3")
  b <- boundary_continuous(des, 12, c(0, 25), c(15, 20, 30), NULL)
  expect_identical(nrow(b$grid), 6L)
  expect_identical(b$grid$sd_t, rep(c(15, 20, 30), each = 2L))
  expect_identical(b$axes$x, "sd_t")
  expect_identical(b$axes$xlab, "Observed treatment standard deviation")
  expect_identical(b$axes$y, "mean_t")
  expect_null(b$axes$caption)
  dec <- decide_continuous(des, data_continuous(12, 25, 20), NULL)
  expect_equal(b$p_go[4L], dec$prob_go)
  expect_equal(b$p_nogo[4L], dec$prob_nogo)
  b2 <- boundary_continuous(des, c(20, 12), 25, 20, NULL)
  expect_identical(b2$grid$n_t, c(12L, 20L))
  expect_identical(b2$axes$x, "n_t")
  expect_identical(b2$axes$xlab, "Treatment group size")
  expect_equal(b2$p_go[1L], dec$prob_go)
  dec20 <- decide_continuous(des, data_continuous(20, 25, 20), NULL)
  expect_equal(b2$p_go[2L], dec20$prob_go)
})

test_that("boundary_continuous() validates its arguments", {
  des <- ref_design_cont("C3")
  expect_error(boundary_continuous(des, 12, 25, NULL, NULL),
               "'mean_t' and 'sd_t' must be given for a continuous endpoint.")
  expect_error(boundary_continuous(des, c(10, 12), 25, c(15, 20), NULL),
               "Only one of 'n_t' and 'sd_t' may have more than one value.")
  expect_error(boundary_continuous(des, 12, 25, 0, NULL),
               "'sd_t' must be greater than 0.")
  expect_error(boundary_continuous(des, 1, 25, 20, NULL),
               "'n_t' must be at least 2.")
  expect_error(boundary_continuous(ref_design_cont("C1"), 12, 25, 20,
                                   data_binary(12, 3)),
               "'data_c' must be created by data_continuous\\(\\).")
  expect_error(boundary_continuous(ref_design_cont("C1", "MM"), 4, 25, 20,
                                   data_continuous(12, 0, 25)),
               "\\(here 3 and 11\\)")
})
