#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector simulate_renewal(std::vector<double> time_vec, std::vector<double> modulant_vec,
                                    double mu, double sigma, double Q) {

  Environment flexsurv = Environment::namespace_env("flexsurv");

  Function dgengamma = flexsurv["dgengamma"];
  Function pgengamma = flexsurv["pgengamma"];

  std::vector<double> event_times;
  event_times.reserve(1000);

  double t_diff = time_vec[1] - time_vec[0];
  double last_event = 1.0;
  event_times.push_back(last_event);

  for (int ti = 1; ti < time_vec.size(); ti++) {

    double pdf = as<double>(dgengamma(time_vec[ti] - last_event, _["mu"] = modulant_vec, _["sigma"] = sigma, _["Q"] = Q));
    double cdf = as<double>(pgengamma(time_vec[ti] - last_event, _["mu"] = modulant_vec, _["sigma"] = sigma, _["Q"] = Q));

    double h_t = t_diff * (pdf / (1.0 - cdf));


    double accept_p = R::unif_rand();
    if (accept_p <= h_t) {
      event_times.push_back(time_vec[ti]);
      last_event = time_vec[ti];
    }
  }
  return(wrap(event_times));
}

