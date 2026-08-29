#include <Rcpp.h>
#include "gengamma_orig.h" // based on code by Christopher Jackson et al. (flexsurv package)

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector simulate_renewal_orig(std::vector<double> time_vec, std::vector<double> modulant_vec,
                                    double shape, double k) {

  std::vector<double> event_times;
  event_times.reserve(1000);

  double t_diff = time_vec[1] - time_vec[0];
  double last_event = -0.05; //R::rexp(modulant_vec[0]);

  event_times.push_back(last_event);

  gengamma_orig::density gen_pdf;
  gengamma_orig::cdf gen_cdf;

  int n = time_vec.size();
  int ti = 0;

  double H_threshold = R::rexp(1);
  double H_level = 0;

  while (ti < n) {
    double dt = time_vec[ti] - last_event;

    if (dt < 0) {
      ti++;
      continue;
    }

    double current_scale = modulant_vec[ti];

    double pdf = gen_pdf(dt, current_scale, shape, k);
    double cdf = gen_cdf(dt, current_scale, shape, k);

    double surv = std::max(1.0 - cdf, 1e-15);
    double h_t = t_diff * pdf / surv;

    H_level += h_t;

    if (H_level >= H_threshold) {
      event_times.push_back(time_vec[ti]);
      last_event = time_vec[ti];

      H_level = 0;
      H_threshold = R::rexp(1);
    }

    ti++;
  }

  return wrap(event_times);
}
