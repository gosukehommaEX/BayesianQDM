test_that("boundary_continuous2() matches the fixed-control references", {
  # Reference: Python_reference/ref_fixed_tte_qdm.py, designs KF1 and KF2
  cov_kf1 <- matrix(c(6400, 190, 190, 36), 2L)
  b <- boundary_continuous2(ref_design_cont2("KF1"), 10, list(c(5, 0), 2),
                            cov_kf1, NULL)
  expect_identical(names(b$grid), c("n_t", "mean_t1", "mean_t2"))
  expect_identical(b$grid$mean_t1, c(0, 5))
  expect_equal(b$p_go[2L], 0.448613893688, tolerance = 1e-9)
  expect_equal(b$p_nogo[2L], 0.104726297391, tolerance = 1e-9)
  dec <- decide_continuous2(ref_design_cont2("KF1"),
                            data_continuous2(10, c(0, 2), cov_kf1), NULL)
  expect_equal(b$p_go[1L], dec$prob_go)
  expect_equal(b$p_nogo[1L], dec$prob_nogo)
  expect_identical(b$axes$x, "mean_t1")
  expect_identical(b$axes$y, "mean_t2")
  expect_null(b$axes$caption)
  b2 <- boundary_continuous2(ref_design_cont2("KF2"), 4, list(2, 1),
                             matrix(c(4, 1, 1, 2), 2L), NULL)
  expect_equal(b2$p_go, 0.545279559813, tolerance = 1e-9)
  expect_equal(b2$p_nogo, 0.13754578381, tolerance = 1e-9)
})

test_that("boundary_continuous2() matches a concurrent-control reference", {
  # Reference: Python_reference/ref_continuous2_qdm.py, design K1
  des <- ref_design_cont2("K1")
  b <- boundary_continuous2(des, 12, list(18, c(2.5, 0)), cov2(60, 5, 0.4),
                            ref_data_cont2$K1$c)
  expect_equal(b$p_go[2L], 0.560892392845125, tolerance = 1e-9)
  expect_equal(b$p_nogo[2L], 0.0606158356586099, tolerance = 1e-9)
  expect_equal(b$grid$diff1, c(28, 28))
  expect_equal(b$grid$diff2, c(1, 3.5))
  expect_identical(b$axes$x, "diff1")
  expect_identical(b$axes$y, "diff2")
  expect_match(b$axes$caption, "^Control data: ")
})

test_that("boundary_continuous2() returns two probabilities for a marginal rule", {
  # Reference: Python_reference/ref_continuous2_qdm.py, design K5 (rule
  # "and"; Go sets R1 to R3 and R1, R4, R7; NoGo sets R7 to R9 and R3, R6,
  # R9)
  ref <- c(0.320405184954998, 0.0837818501873882, 0.132197515292659,
           0.0286780743996835, 0.0109233687877224, 0.0206147670294312,
           0.150916740645318, 0.0694518911782248, 0.183030607524575)
  b <- boundary_continuous2(ref_design_cont2("K5"), 12, list(3, 1),
                            cov2(75, 5.5, 0.4), NULL)
  expect_length(b$p_go, 2L)
  expect_equal(b$p_go[[1L]], sum(ref[1:3]), tolerance = 1e-9)
  expect_equal(b$p_go[[2L]], sum(ref[c(1, 4, 7)]), tolerance = 1e-9)
  expect_equal(b$p_nogo[[1L]], sum(ref[7:9]), tolerance = 1e-9)
  expect_equal(b$p_nogo[[2L]], sum(ref[c(3, 6, 9)]), tolerance = 1e-9)
})

test_that("boundary_continuous2() validates its arguments", {
  des <- ref_design_cont2("K4")
  cv <- cov2(30, 3, 0.5)
  expect_error(boundary_continuous2(des, 15, c(1, 2), cv, NULL),
               "'mean_t' must be a list of two")
  expect_error(boundary_continuous2(des, 15, list(1, 2), NULL, NULL),
               "'cov_t' must be given.")
  expect_error(boundary_continuous2(des, c(10, 15), list(1, 2), cv, NULL),
               "'n_t' must have length 1.")
  expect_error(boundary_continuous2(des, 15, list(1, 2), diag(-1, 2), NULL),
               "'cov_t' must be positive definite.")
  expect_error(boundary_continuous2(ref_design_cont2("K1"), 12, list(1, 2),
                                    cv, ref_data_cont$C1$c),
               "'data_c' must be created by data_continuous2\\(\\).")
  expect_error(boundary_continuous2(ref_design_cont2("K1"), 5, list(1, 2),
                                    cv, ref_data_cont2$K1$c),
               "\\(here 3 and 10\\)")
})
