/*
 Copyright (C) 2011 University of Southern California
 Authors: Andrew D. Smith, Song Qiang
 
 This file is part of rmap.
 
 rmap is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 rmap is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with rmap; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <cmath>
#include <limits>
#include <iomanip>
#include <numeric>
#include <iostream>
#include <sstream>

#include "distribution.hpp"

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>


using std::vector;
using std::pair;
using std::setw;
using std::max;
using std::min;
using std::cerr;
using std::endl;
using std::string;
using std::setprecision;


//////////////////////////////////////////////
//////       numeric compute            //////
//////////////////////////////////////////////

inline static double
sign(double x) {
  return (x >= 0) ? 1.0 : -1.0;
}

//static const double tolerance = 1e-10;

inline static double
invpsi(const double tolerance, const double x) {
  double L = 1.0, Y = std::exp(x);
  while (L > tolerance) {
    Y += L*sign(x - gsl_sf_psi(Y));
    L /= 2.0;
  }
  return Y;
}

static double
movement(const double curr, const double prev) {
  return std::abs(curr - prev)/std::max(std::fabs(curr), std::fabs(prev));
}

//////////////////////////////////////////////
//////       struct BetaBin             //////
//////////////////////////////////////////////

string
BetaBin::tostring() const {
  std::ostringstream os;
  os << setprecision(4) << alpha << " " << setprecision(4) << beta;
  return os.str();
}

double
BetaBin::operator()(const pair<double, double> &val) const
{
  const size_t x = static_cast<size_t>(val.first);
  const size_t n = static_cast<size_t>(x + val.second);
  return gsl_sf_lnchoose(n, x) +
  gsl_sf_lnbeta(alpha + x, beta + val.second) - lnbeta_helper;
}

void
BetaBin::fit(const vector<double> &vals_a, const vector<double> &vals_b,
             const vector<double> &p) {
  const double p_total = std::accumulate(p.begin(), p.end(), 0.0);
  const double alpha_rhs = inner_product(vals_a.begin(), vals_a.end(),
                                         p.begin(), 0.0)/p_total;
  const double beta_rhs = inner_product(vals_b.begin(), vals_b.end(),
                                        p.begin(), 0.0)/p_total;
  double prev_alpha = 0.0, prev_beta = 0.0;
  alpha = beta = 0.01;
  while (movement(alpha, prev_alpha) > tolerance &&
         movement(beta, prev_beta) > tolerance) {
    prev_alpha = alpha;
    prev_beta = beta;
    alpha = invpsi(tolerance, gsl_sf_psi(prev_alpha + prev_beta) + alpha_rhs);
    beta = invpsi(tolerance, gsl_sf_psi(prev_alpha + prev_beta) + beta_rhs);
  }
  lnbeta_helper = gsl_sf_lnbeta(alpha, beta);
}

/*
double
BetaBin::log_likelihood(const pair<double, double> &val) const
{
  const size_t x = static_cast<size_t>(val.first);
  const size_t n = static_cast<size_t>(x + val.second);
  return gsl_sf_lnchoose(n, x) +
  gsl_sf_lnbeta(alpha + x, beta + val.second) - lnbeta_helper;
}*/

//////////////////////////////////////////////
//////       struct NegBin              //////
//////////////////////////////////////////////

void
NegBin::set_helpers() {
  r_log_p_minus_lngamma_r = r*log(p) - gsl_sf_lngamma(r);
  r_log_p_minus_lngamma_r = log(1 - p);
}
