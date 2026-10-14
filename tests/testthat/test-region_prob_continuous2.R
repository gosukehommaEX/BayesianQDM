# References: Python_reference/ref_continuous2_qdm.py, section C (moment
# matching: Theorem 3 of the Moment_matching article with rectangle
# probabilities of the bivariate t distribution by adaptive quadrature;
# exact probabilities by Monte Carlo simulation with 2e7 draws). Regions
# R1 to R9 of the nine-region grid; the predictive designs K3 and K6 use
# R1, R3, R7 and R9.
ref_dist_cont2 <- list(
  K1_t = list(
    loc = matrix(c(18, 2.5), 1L), df = 10,
    scale = matrix(c(330, 11, 2.29166666666667), 1L)
  ),
  K1_c = list(
    loc = matrix(c(-10, -1), 1L), df = 10,
    scale = matrix(c(449.166666666667, 7.7, 3.3), 1L)
  ),
  K2_t = list(
    loc = matrix(c(4.98338870431894, 1.99667774086379), 1L), df = 33,
    scale = matrix(c(193.297585517318, 5.62116011161827, 1.08738514826748), 1L)
  ),
  K2_c = list(
    loc = matrix(c(-10, -1), 1L), df = 33,
    scale = matrix(c(170.668478808014, 4.28495922681969, 1.26975737440854), 1L)
  ),
  K3_t = list(
    loc = matrix(c(2.37257142857143, 1.00118633540373), 1L), df = 19,
    scale = matrix(c(1.67795549308551, 0.263590024497268, 0.60288915604117), 1L)
  ),
  K3_c = list(
    loc = matrix(c(-0.67910447761194, -1.32089552238806), 1L), df = 23,
    scale = matrix(c(1.30581977910577, 0.095846627989191,
                     0.600695357734365), 1L)
  ),
  K4_t = list(
    loc = matrix(c(4, 1.5), 1L), df = 13,
    scale = matrix(c(64.6153846153846, 3.23076923076923, 0.646153846153846), 1L)
  ),
  K4_c = list(
    loc = matrix(c(-10, -1), 1L), df = 13,
    scale = matrix(c(96.9230769230769, 4.84615384615385, 0.969230769230769), 1L)
  ),
  K5_t = list(
    loc = matrix(c(2.97520661157025, 1), 1L), df = 15,
    scale = matrix(c(376.17571659495, 10.0826446280992, 2.03168044077135), 1L)
  ),
  K5_c = list(
    loc = matrix(c(-10, -1), 1L), df = 8,
    scale = matrix(c(581.875, 12.46875, 2.96875), 1L)
  ),
  K6_t = list(
    loc = matrix(c(7, 2), 1L), df = 13,
    scale = matrix(c(9.87179487179487, -0.592307692307692,
                     0.394871794871795), 1L)
  ),
  K6_c = list(
    loc = matrix(c(0, 0), 1L), df = 13,
    scale = matrix(c(21.5384615384615, -1.29230769230769,
                     0.861538461538462), 1L)
  )
)

ref_mm_cont2 <- list(
  K1 = c(0.560892392845125, 0.0741296786423965, 0.0922152459136311,
         0.034539315250446, 0.00700495838687196, 0.0104135124279007,
         0.128560384036125, 0.0316286768388945, 0.0606158356586099),
  K2 = c(0.487770318626615, 0.0761832877719272, 0.0371748442271904,
         0.0670087840509609, 0.0173111723035729, 0.0105685362428432,
         0.182835577836786, 0.0664713667612531, 0.0546761121788512),
  K3 = c(0.626051840857435, 0, 0.0991345612327038, 0, 0, 0, 0.213908233621874,
         0, 0.060905364287987),
  K4 = c(0.475440903771989, 0.104478331943068, 0.0380887134173634,
         0.0739441049062476, 0.0374467301236908, 0.0200135306654595,
         0.0968179703715022, 0.0784964605938997, 0.0752732542067799),
  K5 = c(0.320405184954998, 0.0837818501873882, 0.132197515292659,
         0.0286780743996835, 0.0109233687877224, 0.0206147670294312,
         0.150916740645318, 0.0694518911782248, 0.183030607524575),
  K6 = c(0.478623343701189, -6.93889390390723e-18, 0.154954023476356, 0, 0, 0,
         0.322812062562545, 0, 0.0436105702599108)
)

test_that("region_prob_continuous2() computes moment-matching probabilities", {
  for (name in names(ref_mm_cont2)) {
    crit <- ref_design_cont2(name)$criterion
    dt <- ref_dist_cont2[[paste0(name, "_t")]]
    dc <- ref_dist_cont2[[paste0(name, "_c")]]
    if (crit$type == "posterior") {
      p <- region_prob_continuous2(crit, dt, dc, as.list(1:9), "MM")
      expect_equal(p[1L, ], ref_mm_cont2[[name]], tolerance = 1e-9,
                   info = name)
    } else {
      p <- region_prob_continuous2(crit, dt, dc, as.list(1:4), "MM")
      expect_equal(p[1L, ], ref_mm_cont2[[name]][c(1, 3, 7, 9)],
                   tolerance = 1e-9, info = name)
    }
  }
  # Unions of regions, and one-row distributions recycled
  crit <- ref_design_cont2("K1")$criterion
  dt <- ref_dist_cont2$K1_t
  dt2 <- list(loc = rbind(dt$loc, dt$loc + c(5, 1)),
              scale = rbind(dt$scale, dt$scale), df = dt$df)
  p <- region_prob_continuous2(crit, dt2, ref_dist_cont2$K1_c,
                               list(1:3, c(1L, 4L, 7L), 1:9), "MM")
  expect_identical(dim(p), c(2L, 3L))
  expect_equal(p[1L, ], c(sum(ref_mm_cont2$K1[1:3]),
                          sum(ref_mm_cont2$K1[c(1, 4, 7)]), 1),
               tolerance = 1e-9)
  # A larger treatment mean makes P(theta_1 > TV_1) larger
  expect_gt(p[2L, 1L], p[1L, 1L])
})

test_that("region_prob_continuous2() estimates probabilities by Monte Carlo", {
  ref_mc <- list(
    K1 = c(0.5604494, 0.0740703, 0.0923529, 0.034412, 0.007025, 0.01048045,
           0.12867035, 0.03174825, 0.06079135),
    K3 = c(0.62616815, 0, 0.09916835, 0, 0, 0, 0.2137463, 0, 0.0609172)
  )
  n_mc <- 200000L
  for (name in names(ref_mc)) {
    crit <- ref_design_cont2(name)$criterion
    dt <- ref_dist_cont2[[paste0(name, "_t")]]
    dc <- ref_dist_cont2[[paste0(name, "_c")]]
    n_reg <- if (crit$type == "posterior") 9L else 4L
    set.seed(11)
    p <- region_prob_continuous2(crit, dt, dc, as.list(seq_len(n_reg)), "MC",
                                 n_mc)[1L, ]
    ref <- ref_mc[[name]]
    if (n_reg == 4L) {
      ref <- ref[c(1, 3, 7, 9)]
    }
    se <- sqrt(ref * (1 - ref) * (1 / n_mc + 1 / 2e7))
    expect_true(all(abs(p - ref) <= 4.5 * se + 1e-12), info = name)
    expect_equal(sum(p), 1, tolerance = 1e-12)
  }
})
