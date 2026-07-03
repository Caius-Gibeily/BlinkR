#ifndef GENGAMMA_ORIG_H
#define GENGAMMA_ORIG_H

#include <cmath>
#include <iostream>
#include <boost/math/distributions/gamma.hpp>
#include <algorithm>

#include "Rcpp.h"

namespace gengamma_orig {

inline bool bad(const double scale,
                const double shape,
                const double k) {
  if (scale <= 0 || shape <= 0 || k <= 0) {
    Rcpp::warning("Invalid generalized gamma parameters: scale, shape, and k must be positive.");
    return true;
  }
  return false;
}

class density {
public:
  typedef double result_type;

  inline double operator()(const double x,
                         const double scale,
                         const double shape,
                         const double k) const {

    if (bad(scale, shape, k)) return double(R_NaN);

    if (x <= 0.0) return 0.0;

    const double logdens = std::log(shape) - std::lgamma(k) +
      (shape * k - 1.0) * std::log(x) - (shape * k) * std::log(scale) -
      std::pow((x / scale), shape);

    return std::exp(logdens);
  }
};

class cdf {
public:
  typedef double result_type;

  inline double operator()(const double x,
                         const double scale,
                         const double shape,
                         const double k) const {

    if (bad(scale, shape, k)) return double(R_NaN);

    if (x <= 0.0) return 0.0;

    boost::math::gamma_distribution<double> pgamma(k, 1.0);

    const double y = std::log(x);
    const double w = (y - std::log(scale)) * shape;

    double exp_w = std::exp(w);
    if (std::isinf(exp_w)) return 1.0;

    return boost::math::cdf(pgamma, exp_w);
  }
};

}

#endif
