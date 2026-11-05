# Value of the density curve at the grid point nearest to v (the thresholds
# are points of the grid)
density_at <- function(curve, v) {
  curve$y[which.min(abs(curve$x - v))]
}

# Area under a curve by the trapezoidal rule
trapezoid <- function(x, y) {
  sum(diff(x) * (y[-1L] + y[-length(y)]) / 2)
}

test_that("plot_decision_density() draws the density with the regions", {
  dec <- qdm_decide(ref_design("D1"), data_binary(12, 9), data_binary(12, 4))
  p <- plot(dec)
  b <- ggplot2::ggplot_build(p)
  # Layers: shaded regions, density, thresholds, their labels
  ribbon <- b$data[[1]]
  curve <- b$data[[2]]
  expect_setequal(unique(ribbon$fill), c("#658D1B", "#D91E49"))
  expect_equal(b$data[[3]]$xintercept, c(0.1, 0.2))
  go <- ribbon[ribbon$fill == "#658D1B", ]
  nogo <- ribbon[ribbon$fill == "#D91E49", ]
  expect_gte(min(go$x), 0.2)
  expect_lte(max(nogo$x), 0.1)
  expect_identical(plot_labels(p)$title, "Decision: Go")
  expect_identical(plot_labels(p)$subtitle, decision_subtitle(dec))
  expect_null(plot_labels(p)$caption)
})

test_that("plot_decision_density() agrees with independent references", {
  # References: Python_reference/ref_plot_decision.out (numerical
  # integration of the convolutions with scipy)
  cases <- list(
    D1 = list(dec = qdm_decide(ref_design("D1"), data_binary(12, 9),
                               data_binary(12, 4)),
              tv = 0.2, mav = 0.1, d_tv = 1.41734534210601,
              d_mav = 0.737202805086551, go = 0.820732075666572,
              nogo = 0.0729266458410321),
    D8 = list(dec = qdm_decide(ref_design("D8"), data_binary(20, 9)),
              tv = 0.2, mav = 0.05, d_tv = 3.35450803808753,
              d_mav = 0.568275766017905, go = 0.691441812114726,
              nogo = 0.0206296046344505),
    C1 = list(dec = qdm_decide(ref_design_cont("C1"), ref_data_cont$C1$t,
                               ref_data_cont$C1$c),
              tv = 30, mav = 15, d_tv = 0.0369316305634321,
              d_mav = 0.0136999196050039, go = 0.5,
              nogo = 0.0878402872646904),
    C7 = list(dec = qdm_decide(ref_design_cont("C7"), ref_data_cont$C7$t),
              tv = 30, mav = 15, d_tv = 0.0454718621061924,
              d_mav = 0.000808356477586253, go = 0.797515237746117,
              nogo = 0.00264736615439019),
    T2 = list(dec = qdm_decide(ref_design_tte("T2"), ref_data_tte$T2$t,
                               ref_data_tte$T2$c),
              tv = 0.15, mav = 0, d_tv = 2.24325856369969,
              d_mav = 0.982117879271843, go = 0.665455253873119,
              nogo = 0.0899531660406809),
    H1 = list(dec = qdm_decide(ref_design_hr("H1"), ref_data_hr$t,
                               ref_data_hr$c),
              tv = 0.7, mav = 1, d_tv = 1.20780573200716,
              d_mav = 0.250549025251693, go = 0.768463215969562,
              nogo = 0.0418231222288774)
  )
  for (nm in names(cases)) {
    a <- cases[[nm]]
    # The probabilities of the decision agree with the references
    expect_equal(unname(c(a$dec$prob_go, a$dec$prob_nogo)), c(a$go, a$nogo),
                 tolerance = 1e-6, label = nm)
    b <- ggplot2::ggplot_build(plot(a$dec))
    curve <- b$data[[2]]
    expect_equal(density_at(curve, a$tv), a$d_tv, tolerance = 1e-4,
                 label = paste(nm, "density at TV"))
    expect_equal(density_at(curve, a$mav), a$d_mav, tolerance = 1e-4,
                 label = paste(nm, "density at MAV"))
    ribbon <- b$data[[1]]
    go <- ribbon[ribbon$fill == "#658D1B", ]
    nogo <- ribbon[ribbon$fill == "#D91E49", ]
    expect_lt(abs(trapezoid(go$x, go$ymax) - a$go), 1e-3)
    expect_lt(abs(trapezoid(nogo$x, nogo$ymax) - a$nogo), 1e-3)
  }
})

test_that("plot_decision_density() shows the regions for smaller effects", {
  # Hazard ratio: Go below TV = 0.7 and NoGo above MAV = 1
  dec <- qdm_decide(ref_design_hr("H1"), ref_data_hr$t, ref_data_hr$c)
  ribbon <- ggplot2::ggplot_build(plot(dec))$data[[1]]
  expect_lte(max(ribbon$x[ribbon$fill == "#658D1B"]), 0.7)
  expect_gte(min(ribbon$x[ribbon$fill == "#D91E49"]), 1)
})

test_that("plot_decision_density() names a fixed control in the caption", {
  dec <- qdm_decide(ref_design("D8"), data_binary(20, 9))
  expect_identical(plot_labels(plot(dec))$caption,
                   "Fixed control value: 0.2")
  dec2 <- qdm_decide(ref_design_hr("H2"), ref_data_hr$t)
  expect_identical(plot_labels(plot(dec2))$caption,
                   "Fixed control median survival time: 6")
})
