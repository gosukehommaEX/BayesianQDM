# Exact references for qdm_find_n(): Python_reference/ref_find_n.out
# (rational arithmetic). Design D1 of helper-designs.R with n_t = n_c = n,
# the scenarios of cal_args$D1 (helper-calibration.R) and the targets 0.1.
ref_find_n <- data.frame(
  n = seq(6, 20, by = 2),
  prob_go = c(0.18844213248, 0.2186921927245824, 0.12643930365191258,
              0.14854199690661835, 0.0878492708944844, 0.10284600794598464,
              0.11640471498942191, 0.07267046857590724),
  prob_nogo = c(0.083356180992, 0.0659468285608192, 0.05214455192408728,
                0.041346697310887505, 0.03287510948470714,
                0.026210792900087765, 0.020934516331922774,
                0.016762286879802333)
)

find_n_d1 <- function(...) {
  a <- cal_args$D1
  qdm_find_n(ref_design("D1"), a$scenario_go, a$scenario_nogo,
             target_go = 0.1, target_nogo = 0.1, ...)
}
