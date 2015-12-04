/*
 Copyright (C) 2015-2016 University of Southern California
 Authors: Andrew D. Smith, Xiaojing Ji
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
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
using std::isfinite;
using std::abs;

typedef vector< vector<double> > matrix;


//////////////////////////////////////////////
//////       numeric compute            //////
//////////////////////////////////////////////

inline double
log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  if (!isfinite(p) && !isfinite(q)) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}


inline double
log_sub_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  if (!isfinite(p) && !isfinite(q)) {return p;}
  return p + log(1.0 - exp(q - p));
}


inline static double
sign(double x) {
  return (x >= 0) ? 1.0 : -1.0;
}


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
  double y;
  if (val.second == -1) { // imputated CpG -- beta distribution
    const double v_smooth = min(max(val.first, 1e-2), 1.0 - 1e-2);
    y = (alpha - 1) * log(v_smooth) + (beta - 1) * log(1 - v_smooth)
        - lnbeta_helper;
  } else if (val.second == -2) {
    y = 1;
  } else {
    const size_t x = static_cast<size_t>(val.first);
    const size_t n = static_cast<size_t>(x + val.second);
    y = gsl_sf_lnchoose(n, x) + gsl_sf_lnbeta(alpha + x, beta + val.second)
        - lnbeta_helper;
  }
  return y;
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


//////////////////////////////////////////////
//////       struct NegBin              //////
//////////////////////////////////////////////

void
NegBin::set_helpers() {
  r_log_p_minus_lngamma_r = r*log(p) - gsl_sf_lngamma(r);
  r_log_p_minus_lngamma_r = log(1 - p);
}


//////////////////////////////////////////////
//////       struct CTHMM duration      //////
//////////////////////////////////////////////


void ExpTransEstimator::calc_internal_data(const vector<size_t> &t,
                                           const double v,
                                           vector<double> &ebt) const {
  // store ebt for the computation convenience
  for (size_t i = 0; i < t.size(); ++i) {
    ebt[i] = exp(- v * t[i]);
  }
}



double ExpTransEstimator::calc_llh(const matrix &r, const double u,
                                   const vector<double> &ebt) const {
  
  double val = 0;
  for (size_t i = 0; i < ebt.size(); ++i) {
    const double d = r[0][i]*log(1 - u + u * ebt[i])
                     + r[1][i]*log(u - u * ebt[i])
                     + r[2][i]*log(1 - u - (1 - u) * ebt[i])
                     + r[3][i]*log(u + (1 - u) * ebt[i]);
    
    val += d;
  }
  return val;
}


double ExpTransEstimator::llh_grad_a(const matrix &r,
                                     const vector<double> &ebt) const {
  double val = 0;
  for (size_t i = 0; i < ebt.size(); ++i) {
    val += r[0][i]*( (-1+ebt[i]) / (1-a+a*ebt[i]) )
           + r[1][i]*( (1) / (a) )
           + r[2][i]*( (1) / (1-a) )
           + r[3][i]*( (1-ebt[i]) / (a+(1-a)*ebt[i]) );
  }
  return val;
}


double ExpTransEstimator::llh_grad_b(const matrix &r,
                                     const vector<double> &ebt) const {
  double val = 0;
  for (size_t i = 0; i < ebt.size(); ++i) {
    val += r[0][i]*( (-a*b*ebt[i]) / (1-a+a*ebt[i]) )
           + (r[1][i] + r[2][i]) * ( (b*ebt[i]) / (1-ebt[i]) )
           + r[3][i]*( (-(1-a)*b*ebt[i]) / (a+(1-a)*ebt[i]) );
  }
  return val;
}


void ExpTransEstimator::GA_stepforward(const double grad_a,
                                       const double grad_b,
                                       const double &old_llh, double &new_llh,
                                       const matrix &r,
                                       const vector<size_t> &t)  {
  
  double try_step = 0.1 / abs(grad_a);
  
  double new_a = a + try_step * grad_a;
  double new_b = b + try_step * grad_b;
  
  double moving_llh = - std::numeric_limits<double>::max();
  
  vector<double> ebt (t.size(), 0);
  if (new_a > 0 && new_a < 1 && new_b > 0) {
    calc_internal_data(t, new_b, ebt);
    moving_llh = calc_llh(r, new_a, ebt);
  }
  
  size_t itr = 1;
  
  while (moving_llh <= old_llh && abs(moving_llh - old_llh) > tolerance
         && itr < max_iteration) {
    try_step = try_step / 2;
    new_a = a + try_step * grad_a;
    new_b = b + try_step * grad_b;
    
    if (new_a > 0 && new_a < 1 && new_b > 0) {
      calc_internal_data(t, new_b, ebt);
      moving_llh = calc_llh(r, new_a, ebt);
    }
    
    ++itr;
  }

  if (moving_llh > old_llh && new_a > 0 && new_a < 1 && new_b > 0) {
    a = new_a;
    b = new_b;
    new_llh = moving_llh;
  }
  else {
    new_llh = old_llh;
  }
}


void ExpTransEstimator::mle_GradAscent(const matrix &r,
                                       const vector<size_t> &t) {
  
  vector<double> ebt (t.size(), 0);
  calc_internal_data(t, b, ebt);
  double curr_llh = calc_llh(r, a, ebt);
  
  double grad_a, grad_b;
  grad_a = llh_grad_a(r, ebt);
  grad_b = llh_grad_b(r, ebt);
  
  double new_llh;
  GA_stepforward(grad_a, grad_b, curr_llh, new_llh, r, t);
  
  size_t itr = 1;
  
  while (abs(new_llh - curr_llh) > tolerance && itr < max_iteration) {
    curr_llh = new_llh;
    calc_internal_data(t, b, ebt);
    grad_a = llh_grad_a(r, ebt);
    grad_b = llh_grad_b(r, ebt);
    GA_stepforward(grad_a, grad_b, curr_llh, new_llh, r, t);
    ++itr;
  }
}


