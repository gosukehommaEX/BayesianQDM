# Designs shared by the tests of the continuous endpoint. C1 to C6 are the
# designs of Python_reference/ref_continuous_qdm.py, whose results (closed
# forms, two numerical integrations and noncentral t distributions) are used
# as references; C7 and C8 (fixed control) are those of
# Python_reference/ref_fixed_tte_qdm.py.
crit_cont <- function(gamma_go = 0.6, gamma_nogo = 0.7, tv = 30, mav = 15) {
  criterion_posterior(theta_tv = tv, theta_mav = mav, gamma_go = gamma_go,
                      gamma_nogo = gamma_nogo)
}

ref_design_cont <- function(name, method = "NI", nMC = 10000L) {
  vague <- prior_vague()
  switch(
    name,
    # Concurrent control, vague priors (the setting of Kang et al., 2026,
    # Section 6, without external data)
    C1 = design_continuous(12, 12, vague, vague, crit_cont(),
                           method = method, nMC = nMC),
    # Normal-Inverse-chi-squared priors and external data for both groups
    C2 = design_continuous(10, 8, prior_nix(5, 2, 4, 20),
                           prior_nix(-5, 1, 3, 25), crit_cont(),
                           external_t = borrow(data_continuous(6, 20, 30), 0.5),
                           external_c = borrow(data_continuous(30, -12, 20),
                                               0.4),
                           method = method, nMC = nMC),
    # Hypothetical control
    C3 = design_continuous(12, 0, vague, criterion = crit_cont(),
                           hypothetical_c = hypothetical(mean = -10, r = 1.5),
                           method = method, nMC = nMC),
    # External control only, with the vague prior for the control group
    C4 = design_continuous(12, 0, prior_nix(0, 1, 2, 25), vague, crit_cont(),
                           external_c = borrow(data_continuous(30, -12, 20),
                                               0.4),
                           method = method, nMC = nMC),
    # Predictive criterion with a concurrent control
    C5 = design_continuous(12, 12, vague, vague,
                           criterion_predictive(10, m_t = 40, m_c = 30,
                                                gamma_go = 0.8,
                                                gamma_nogo = 0.7),
                           method = method, nMC = nMC),
    # Predictive criterion with a hypothetical control
    C6 = design_continuous(15, 0, vague,
                           criterion = criterion_predictive(5, m_t = 30,
                                                            m_c = 20,
                                                            gamma_go = 0.7,
                                                            gamma_nogo = 0.6),
                           hypothetical_c = hypothetical(mean = 0, r = 2),
                           method = method, nMC = nMC),
    # Fixed control, posterior criterion
    C7 = design_continuous(12, 0, vague, criterion = crit_cont(),
                           fixed_c = -10, method = method, nMC = nMC),
    # Fixed control, predictive criterion
    C8 = design_continuous(15, 0, vague,
                           criterion = criterion_predictive(5, m_t = 30,
                                                            gamma_go = 0.7,
                                                            gamma_nogo = 0.6),
                           fixed_c = 0, method = method, nMC = nMC)
  )
}

# Observed data of the reference cases (Python_reference/ref_continuous_qdm.py)
ref_data_cont <- list(
  C1 = list(t = data_continuous(12, 18, 25), c = data_continuous(12, -12, 25)),
  C2 = list(t = data_continuous(10, 22, 24), c = data_continuous(8, -8, 18)),
  C3 = list(t = data_continuous(12, 25, 20), c = NULL),
  C4 = list(t = data_continuous(12, 20, 22), c = NULL),
  C5 = list(t = data_continuous(12, 18, 25), c = data_continuous(12, -12, 25)),
  C6 = list(t = data_continuous(15, 9, 12), c = NULL),
  C7 = list(t = data_continuous(12, 25, 20), c = NULL),
  C8 = list(t = data_continuous(15, 9, 12), c = NULL)
)

scen_cont_c <- data.frame(mu_t = c(0, 20, 40), sd_t = 25, mu_c = 0, sd_c = 25)
scen_cont_u <- data.frame(mu_t = c(0, 20, 40), sd_t = 20)
