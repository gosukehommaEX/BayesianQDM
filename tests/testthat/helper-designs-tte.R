# Designs shared by the tests of the time-to-event endpoint. They are the
# designs T1 to T6 of Python_reference/ref_fixed_tte_qdm.py, whose results
# (Gamma distribution functions, an independent numerical integration of
# the distribution of the difference of milestone survival probabilities,
# beta prime distributions of the total observed time with complete
# follow-up) are used as references.
crit_tte <- function(gamma_go = 0.6, gamma_nogo = 0.8, tv = 0.15, mav = 0) {
  criterion_posterior(theta_tv = tv, theta_mav = mav, gamma_go = gamma_go,
                      gamma_nogo = gamma_nogo)
}

# Complete follow-up: enrollment within one time unit, no dropout and an
# analysis long after the last event, so every subject has an event and the
# total observed time of n subjects is Gamma(n, hazard)
conduct_complete <- conduct_tte(accrual_time = c(0, 1), cut_time = 1e6)
# Enrollment with increasing rates over 15 time units, a small dropout hazard
# and an analysis at time 15.44
conduct_example <- conduct_tte(accrual_time = c(0, 5, 10, 15),
                               accrual_rate = c(0.4, 2, 4.5),
                               dropout_hazard = 0.00226, cut_time = 15.44)

ref_design_tte <- function(name, nMC = 10000L) {
  g <- prior_gamma(0.1, 1)
  switch(
    name,
    # Fixed control, posterior criterion
    T1 = design_tte(30, 0, g, criterion = crit_tte(), milestone = 3,
                    conduct = conduct_complete, fixed_c = 0.72),
    # Concurrent control, posterior criterion
    T2 = design_tte(20, 20, g, g, crit_tte(), milestone = 3,
                    conduct = conduct_example),
    # Fixed control, predictive criterion, complete follow-up of the future
    # study
    T3 = design_tte(30, 0, g,
                    criterion = criterion_predictive(0.05, m_t = 10,
                                                     gamma_go = 0.8,
                                                     gamma_nogo = 0.3),
                    milestone = 3, conduct = conduct_complete,
                    fixed_c = 0.55, nMC = nMC),
    # Concurrent control, predictive criterion
    T4 = design_tte(20, 20, g, g,
                    criterion_predictive(0.05, m_t = 10, m_c = 8,
                                         gamma_go = 0.8, gamma_nogo = 0.3),
                    milestone = 3, conduct = conduct_complete, nMC = nMC),
    # Hypothetical control
    T5 = design_tte(30, 0, g, g, crit_tte(), milestone = 3,
                    conduct = conduct_example,
                    hypothetical_c = hypothetical(data_tte(8, 60))),
    # External control only, external treatment data
    T6 = design_tte(30, 0, g, g, crit_tte(), milestone = 3,
                    conduct = conduct_example,
                    external_t = borrow(data_tte(4, 30), 0.5),
                    external_c = borrow(data_tte(20, 100), 0.5))
  )
}

# Observed data of the reference cases
ref_data_tte <- list(
  T1 = list(t = data_tte(6, 150), c = NULL),
  T2 = list(t = data_tte(5, 40), c = data_tte(9, 35)),
  T3 = list(t = data_tte(6, 70), c = NULL),
  T4 = list(t = data_tte(6, 70), c = data_tte(9, 60)),
  T5 = list(t = data_tte(6, 150), c = NULL),
  T6 = list(t = data_tte(6, 150), c = NULL)
)

# Reference probabilities of the Go and NoGo regions
ref_prob_tte <- list(
  T1 = c(go = 0.686259165773, nogo = 0.00106102273767),
  T2 = c(go = 0.665455253873, nogo = 0.0899531660407),
  T3 = c(go = 0.917785289262, nogo = 1 - 0.917785289262),
  T4 = c(go = 0.69356649484, nogo = 1 - 0.69356649484),
  T5 = c(go = 0.713482708146, nogo = 0.0137418340154),
  T6 = c(go = 0.920053158435, nogo = 0.00182044277956)
)

# Designs with the hazard ratio (direction = "less"), whose references are
# computed by Python_reference/ref_hazard_ratio.py (integration over the
# control hazard, exact operating characteristics with complete follow-up,
# beta distributions of the future total observed time).
crit_hr <- function(gamma_go = 0.6, gamma_nogo = 0.7, tv = 0.7, mav = 1) {
  criterion_posterior(theta_tv = tv, theta_mav = mav, gamma_go = gamma_go,
                      gamma_nogo = gamma_nogo, direction = "less")
}

ref_design_hr <- function(name, nMC = 10000L) {
  g <- prior_gamma(0.1, 1)
  pred <- function(m_c = 10L) {
    criterion_predictive(0.8, m_t = 10, m_c = m_c, gamma_go = 0.8,
                         gamma_nogo = 0.3, direction = "less")
  }
  switch(
    name,
    # Concurrent control, complete follow-up
    H1 = design_tte(30, 30, g, g, crit_hr(), effect = "hazard_ratio",
                    conduct = conduct_complete),
    # Fixed control median survival time 6
    H2 = design_tte(30, 0, g, criterion = crit_hr(), effect = "hazard_ratio",
                    conduct = conduct_complete, fixed_c = 6),
    # Hypothetical control
    H3 = design_tte(30, 0, g, g, crit_hr(), effect = "hazard_ratio",
                    conduct = conduct_example,
                    hypothetical_c = hypothetical(data_tte(15, 100))),
    # External control only, external treatment data
    H4 = design_tte(30, 0, g, g, crit_hr(), effect = "hazard_ratio",
                    conduct = conduct_example,
                    external_t = borrow(data_tte(6, 80), 0.5),
                    external_c = borrow(data_tte(30, 250), 0.5)),
    # Predictive criterion, fixed control, complete follow-up of the future
    # study
    P1 = design_tte(30, 0, g, criterion = pred(), effect = "hazard_ratio",
                    conduct = conduct_complete, fixed_c = 6, nMC = nMC),
    # Predictive criterion, concurrent control
    P2 = design_tte(30, 30, g, g, pred(8L), effect = "hazard_ratio",
                    conduct = conduct_complete, nMC = nMC)
  )
}

ref_data_hr <- list(t = data_tte(12, 200), c = data_tte(20, 180))

# Reference probabilities of the Go and NoGo regions
ref_prob_hr <- list(
  H1 = c(go = 0.7684632159696, nogo = 0.04182312222888),
  H2 = c(go = 0.8796664434856, nogo = 0.004229518766598),
  H3 = c(go = 0.9250389365495, nogo = 0.009394232547511),
  H4 = c(go = 0.7863715579954, nogo = 0.03986519327619),
  P1 = c(go = 0.833021537865, nogo = 1 - 0.833021537865),
  P2 = c(go = 0.751839324506, nogo = 1 - 0.751839324506)
)
