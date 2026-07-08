#ifndef GENGAMMA_H
#define GENGAMMA_H

#include <cmath>

#include <Rcpp.h>
#include <boost/math/distributions/gamma.hpp>

#include "distribution.h"

namespace {

namespace gengamma {

class density {
public:
  typedef double result_type;

  inline double operator()(const double x,
                         const double mu,
                         const double sigma,
                         const double Q) const {

    if (Q!=0) {
      const double y = std::log(x);
      const double w = (y - mu) / sigma;
      const double abs_q = std::abs(Q);
      const double qi = 1/(Q * Q);
      const double qw = Q * w;
      return -std::log(sigma*x) +
        std::log(abs_q) * (1 - 2 * qi) +
        qi * (qw - std::exp(qw)) - R::lgammafn(qi);
    } else {
      return R::dlnorm(x, mu, sigma, 1);
    }
  }
};

class cdf {
public:
  typedef double result_type;

  inline double operator()(const double q,
                         const double mu,
                         const double sigma,
                         const double big_q) const {
    /* check the arguments */

    if (big_q!=0) {
      const double y = std::log(q);
      const double w = (y - mu) / sigma;
      const double qq = 1/(big_q * big_q);
      const double expnu = std::exp(big_q * w) * qq;

      boost::math::gamma_distribution<double> pgamma(qq, 1.0);
      return boost::math::cdf(pgamma,expnu);

    } else {
      return R::plnorm(q, mu, sigma,1,0);
    }
  }


};
}
}

#endif
