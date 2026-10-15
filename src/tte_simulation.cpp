// Simulation of studies with a time-to-event endpoint, for the operating
// characteristics (studies simulated under scenarios for the true hazards)
// and for posterior predictive probabilities (future studies simulated with
// hazards drawn from their posterior distributions).
//
// A study enrolls the subjects of both groups at independent times from a
// piecewise uniform distribution on the accrual intervals (with the
// probability of an interval proportional to its rate times its length).
// The event time of a subject is exponential with the hazard of its group,
// and the dropout time is piecewise exponential, independent of the event
// time. The analysis takes place at a fixed calendar time after the start of
// enrollment, or at the calendar time of the k-th event in the whole study
// (all subjects followed until their event or dropout if fewer than k events
// occur). A subject enrolled at e with event time T and dropout time C
// contributes the observed time min(T, C, cut - e) (0 if e > cut) and an
// event if T <= C and e + T <= cut. The summary of a group is its number of
// events and its total observed time.
//
// The random numbers of a subject are drawn in the order: a uniform number
// for the enrollment time, a standard exponential number for the event time
// and a standard exponential number for the dropout time, with the random
// number generator of R. The event time is the standard exponential number
// divided by the hazard, so the same random numbers give the event times
// under all scenarios.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

const double kInf = std::numeric_limits<double>::infinity();

struct Conduct {
  // Accrual: interval boundaries (0 = t_0 < ... < t_K) and cumulative
  // probabilities at t_1, ..., t_K (the last one is 1)
  std::vector<double> acc_time;
  std::vector<double> acc_cum;
  // Dropout: start times of the pieces (the first is 0, the last piece has no
  // end) and their hazards
  std::vector<double> drop_start;
  std::vector<double> drop_hazard;
  // Analysis at cut_time, or at the cut_events-th event if cut_events > 0
  double cut_time;
  int cut_events;
};

Conduct make_conduct(const Rcpp::NumericVector &acc_time,
                     const Rcpp::NumericVector &acc_cum,
                     const Rcpp::NumericVector &drop_start,
                     const Rcpp::NumericVector &drop_hazard, double cut_time,
                     int cut_events) {
  Conduct c;
  c.acc_time.assign(acc_time.begin(), acc_time.end());
  c.acc_cum.assign(acc_cum.begin(), acc_cum.end());
  c.drop_start.assign(drop_start.begin(), drop_start.end());
  c.drop_hazard.assign(drop_hazard.begin(), drop_hazard.end());
  c.cut_time = cut_time;
  c.cut_events = cut_events;
  return c;
}

// Enrollment time for the uniform number u (inverse distribution function of
// the piecewise uniform distribution). Intervals with probability 0 are
// never selected.
double enrollment_time(double u, const Conduct &c) {
  const int k_max = static_cast<int>(c.acc_cum.size());
  double prev = 0.0;
  int last_pos = 0;
  for (int k = 0; k < k_max; ++k) {
    const double cum = c.acc_cum[k];
    if (cum > prev) {
      last_pos = k;
      if (u < cum) {
        return c.acc_time[k] +
          (u - prev) / (cum - prev) * (c.acc_time[k + 1] - c.acc_time[k]);
      }
    }
    prev = cum;
  }
  // u at least the last cumulative probability (rounding): end of the last
  // interval with positive probability
  return c.acc_time[last_pos + 1];  // # nocov
}

// Dropout time for the standard exponential number x (inverse of the
// cumulative hazard of the piecewise exponential distribution); infinite if
// the cumulative hazard stays below x.
double dropout_time(double x, const Conduct &c) {
  const int j_max = static_cast<int>(c.drop_start.size());
  double cum = 0.0;
  for (int j = 0; j < j_max; ++j) {
    const double h = c.drop_hazard[j];
    const double end = j + 1 < j_max ? c.drop_start[j + 1] : kInf;
    if (h > 0.0) {
      if (end == kInf || x <= cum + h * (end - c.drop_start[j])) {
        return c.drop_start[j] + (x - cum) / h;
      }
      cum += h * (end - c.drop_start[j]);
    }
  }
  return kInf;
}

// Random numbers of n subjects: enrollment times, standard exponential
// numbers for the event times and dropout times.
void draw_subjects(int n, const Conduct &c, std::vector<double> &enroll,
                   std::vector<double> &base, std::vector<double> &drop) {
  for (int i = 0; i < n; ++i) {
    enroll[i] = enrollment_time(R::unif_rand(), c);
    base[i] = R::exp_rand();
    drop[i] = dropout_time(R::exp_rand(), c);
  }
}

// Numbers of events and total observed times of the two groups (the first
// n_t subjects are in the treatment group) for the hazards lam_t and lam_c.
// out: d_t, R_t, d_c, R_c. event_cal and work are work vectors of length n.
void summarize(int n_t, int n, const std::vector<double> &enroll,
               const std::vector<double> &base,
               const std::vector<double> &drop, double lam_t, double lam_c,
               const Conduct &c, std::vector<double> &event_cal,
               std::vector<double> &work, double out[4]) {
  // Calendar times of the events (infinite if dropout comes first)
  int n_event = 0;
  for (int i = 0; i < n; ++i) {
    const double t = base[i] / (i < n_t ? lam_t : lam_c);
    if (t <= drop[i]) {
      event_cal[i] = enroll[i] + t;
      ++n_event;
    } else {
      event_cal[i] = kInf;
    }
  }
  double cut = c.cut_time;
  if (c.cut_events > 0) {
    if (n_event >= c.cut_events) {
      std::copy(event_cal.begin(), event_cal.begin() + n, work.begin());
      std::nth_element(work.begin(), work.begin() + (c.cut_events - 1),
                       work.begin() + n);
      cut = work[c.cut_events - 1];
    } else {
      cut = kInf;
    }
  }
  out[0] = 0.0;
  out[1] = 0.0;
  out[2] = 0.0;
  out[3] = 0.0;
  for (int i = 0; i < n; ++i) {
    const double follow = cut - enroll[i];
    if (!(follow > 0.0)) {
      continue;
    }
    const double t = base[i] / (i < n_t ? lam_t : lam_c);
    const int g = i < n_t ? 0 : 2;
    out[g + 1] += std::min(std::min(t, drop[i]), follow);
    // A subject without an event (event_cal infinite) is not counted when
    // the analysis time is infinite (fewer events than cut_events)
    if (event_cal[i] < kInf && event_cal[i] <= cut) {
      out[g] += 1.0;
    }
  }
}

}  // namespace

// Summaries of nsim simulated studies under each scenario.
// n_t, n_c: numbers of subjects; hazard_t, hazard_c: hazards of the
// scenarios (hazard_c is not used if n_c = 0).
// acc_time, acc_cum, drop_start, drop_hazard, cut_time, cut_events: conduct
// of the study (see Conduct).
// Returns a matrix with nsim rows per scenario (scenario s in rows
// s * nsim + 1 to (s + 1) * nsim) and the columns d_t, R_t, d_c, R_c. The
// same random numbers are used for all scenarios.
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_tte_oc_data(int nsim, int n_t, int n_c,
                                    Rcpp::NumericVector hazard_t,
                                    Rcpp::NumericVector hazard_c,
                                    Rcpp::NumericVector acc_time,
                                    Rcpp::NumericVector acc_cum,
                                    Rcpp::NumericVector drop_start,
                                    Rcpp::NumericVector drop_hazard,
                                    double cut_time, int cut_events) {
  const Conduct c = make_conduct(acc_time, acc_cum, drop_start, drop_hazard,
                                 cut_time, cut_events);
  const int n = n_t + n_c;
  const int n_scen = static_cast<int>(hazard_t.size());
  std::vector<double> enroll(n);
  std::vector<double> base(n);
  std::vector<double> drop(n);
  std::vector<double> event_cal(n);
  std::vector<double> work(n);
  Rcpp::NumericMatrix out(nsim * n_scen, 4);
  double res[4];
  for (int i = 0; i < nsim; ++i) {
    draw_subjects(n, c, enroll, base, drop);
    for (int s = 0; s < n_scen; ++s) {
      summarize(n_t, n, enroll, base, drop, hazard_t[s],
                n_c > 0 ? hazard_c[s] : 1.0, c, event_cal, work, res);
      for (int k = 0; k < 4; ++k) {
        out(s * nsim + i, k) = res[k];
      }
    }
    if ((i + 1) % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }
  return out;
}

// Posterior predictive probabilities P(future effect > q) by simulation of
// nMC future studies for each posterior distribution.
// shape_t, rate_t: Gamma posterior distributions of the treatment hazard.
// shape_c, rate_c: those of the control hazard (same length as shape_t), or
//   empty vectors for a fixed control with milestone survival fixed_c (the
//   future study then has no control group).
// prior_t, prior_c: shapes and rates of the Gamma priors used to estimate
//   the milestone survival in the future study.
// m_t, m_c: numbers of subjects of the future study; conduct as in
//   cpp_tte_oc_data(); t0: milestone time.
// In each simulation, the hazards are drawn from their posterior
// distributions (treatment, then control), the future study is simulated,
// and the milestone survival of each group is estimated by
// exp(-t0 (a + d) / (b + R)) with the prior (a, b) and the future summary
// (d, R).
// [[Rcpp::export]]
Rcpp::NumericVector cpp_pred_tte(Rcpp::NumericVector shape_t,
                                 Rcpp::NumericVector rate_t,
                                 Rcpp::NumericVector shape_c,
                                 Rcpp::NumericVector rate_c, double fixed_c,
                                 Rcpp::NumericVector prior_t,
                                 Rcpp::NumericVector prior_c, int m_t,
                                 int m_c, Rcpp::NumericVector acc_time,
                                 Rcpp::NumericVector acc_cum,
                                 Rcpp::NumericVector drop_start,
                                 Rcpp::NumericVector drop_hazard,
                                 double cut_time, int cut_events, double t0,
                                 double q, int nMC) {
  const Conduct c = make_conduct(acc_time, acc_cum, drop_start, drop_hazard,
                                 cut_time, cut_events);
  const bool fixed = shape_c.size() == 0;
  const int n = m_t + (fixed ? 0 : m_c);
  const R_xlen_t n_dist = shape_t.size();
  std::vector<double> enroll(n);
  std::vector<double> base(n);
  std::vector<double> drop(n);
  std::vector<double> event_cal(n);
  std::vector<double> work(n);
  Rcpp::NumericVector out(n_dist);
  double res[4];
  for (R_xlen_t i = 0; i < n_dist; ++i) {
    long count = 0;
    for (int b = 0; b < nMC; ++b) {
      const double lam_t = R::rgamma(shape_t[i], 1.0 / rate_t[i]);
      const double lam_c =
        fixed ? 1.0 : R::rgamma(shape_c[i], 1.0 / rate_c[i]);
      draw_subjects(n, c, enroll, base, drop);
      summarize(m_t, n, enroll, base, drop, lam_t, lam_c, c, event_cal, work,
                res);
      const double est_t =
        std::exp(-t0 * (prior_t[0] + res[0]) / (prior_t[1] + res[1]));
      const double est_c = fixed ? fixed_c :
        std::exp(-t0 * (prior_c[0] + res[2]) / (prior_c[1] + res[3]));
      if (est_t - est_c > q) {
        ++count;
      }
    }
    out[i] = static_cast<double>(count) / nMC;
    Rcpp::checkUserInterrupt();
  }
  return out;
}
