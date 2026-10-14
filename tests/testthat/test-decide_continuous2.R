# References: Python_reference/ref_continuous2_qdm.py, sections A and C
# (moment matching; the regions R1 to R9 of the nine-region grid, of which
# the predictive designs K3 and K6 use R1, R3, R7 and R9).
ref_dec_cont2 <- list(
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

test_that("decide_continuous2() matches the references", {
  ref <- ref_dec_cont2
  decisions <- c(K1 = "Go", K2 = "Go", K3 = "Go", K4 = "Go", K5 = "Gray",
                 K6 = "Gray")
  go <- list(K1 = 1, K2 = 1, K3 = 1, K4 = 1:2, K5 = list(1:3, c(1, 4, 7)),
             K6 = 1)
  nogo <- list(K1 = 9, K2 = 9, K3 = 9, K4 = 9, K5 = list(7:9, c(3, 6, 9)),
               K6 = 9)
  for (name in names(ref)) {
    des <- ref_design_cont2(name)
    d <- ref_data_cont2[[name]]
    res <- decide_continuous2(des, d$t, d$c)
    cells <- if (des$criterion$type == "posterior") 1:9 else c(1, 3, 7, 9)
    expect_equal(unname(res$prob_regions), ref[[name]][cells],
                 tolerance = 1e-9, info = name)
    expect_identical(names(res$prob_regions),
                     paste0("R", seq_along(cells)))
    sum_cells <- function(s) {
      if (is.list(s)) {
        return(vapply(s, sum_cells, numeric(1)))
      }
      sum(ref[[name]][s])
    }
    expect_equal(res$prob_go, sum_cells(go[[name]]), tolerance = 1e-9,
                 info = name)
    expect_equal(res$prob_nogo, sum_cells(nogo[[name]]), tolerance = 1e-9,
                 info = name)
    expect_identical(names(res$posterior_t),
                     c("mu1", "mu2", "kappa", "nu", "L11", "L12", "L22"))
    dec <- qdm_decide(des, d$t, d$c)
    expect_identical(dec$decision, decisions[[name]], info = name)
  }
})

test_that("decide_continuous2() checks the data and the degrees of freedom", {
  des <- ref_design_cont2("K1")
  d <- ref_data_cont2$K1
  expect_error(decide_continuous2(des, data_continuous(12, 1, 2), d$c),
               "'data_t' must be created by data_continuous2\\(\\).")
  expect_error(decide_continuous2(des, d$t, NULL),
               "'data_c' must be created by data_continuous2\\(\\) for a")
  expect_error(decide_continuous2(ref_design_cont2("K4"), d$t, d$c),
               "'data_c' must be NULL for a design without a concurrent")
  # Vague prior with 6 subjects: 4 degrees of freedom
  small <- data_continuous2(6, c(1, 2), cov2(60, 5, 0.4))
  expect_error(decide_continuous2(des, small, d$c),
               "\\(here 4 and 10\\). Use method = \"MC\".")
  set.seed(3)
  res <- decide_continuous2(ref_design_cont2("K1", "MC", 50000L), d$t, d$c)
  expect_equal(sum(res$prob_regions), 1, tolerance = 1e-12)
  expect_equal(res$prob_go, res$prob_regions[["R1"]])
  expect_equal(res$prob_nogo, res$prob_regions[["R9"]])
  ref <- ref_dec_cont2$K1
  se <- sqrt(ref * (1 - ref) / 50000)
  expect_true(all(abs(res$prob_regions - ref) <= 4.5 * se + 0.002))
})
