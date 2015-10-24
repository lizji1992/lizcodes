/*
 Copyright (C) 2015-2016 University of Southern California
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

/*
double p_bb(const double t) const { return ; }
double p_bf(const double t) const { return a-a*ebt[i]; }
double p_fb(const double t) const { return 1-a-(1-a)*ebt[i]; }
double p_ff(const double t) const { return a+(1-a)*ebt[i]; }
*/
 

void ExpTransEstimator::calc_internal_data(const vector<size_t> &t,
                                           const double u, const double v,
                                           vector<double> &log_a,
                                           vector<double> &log_1_a,
                                           vector<double> &neg_bt) const {
  // store ebt for the computation convenience
  for (size_t i = 0; i < t.size(); ++i) {
    log_a[i] = log(u);
    log_1_a[i] = log(1-u);
    neg_bt[i] = - v * (t[i+1] - t[i]);
  }
  
}


double ExpTransEstimator::calc_log_llh(const matrix &r,
                                       const vector<double> &log_a,
                                       const vector<double> &log_1_a,
                                       const vector<double> &neg_bt) const {
  
  double lp_bb = log_sum_log(log_1_a[0], log_a[0] + neg_bt[0]);
  double lp_bf = log_sub_log(log_a[0], log_a[0] + neg_bt[0]);
  double lp_fb = log_sub_log(log_1_a[0], log_1_a[0] + neg_bt[0]);
  double lp_ff = log_sum_log(log_a[0], log_1_a[0] + neg_bt[0]);
  
  
  double val = r[0][0] + log(lp_bb);
  val = log_sum_log(val, r[1][0] + log(lp_bf));
  val = log_sum_log(val, r[2][0] + log(lp_fb));
  val = log_sum_log(val, r[3][0] + log(lp_ff));
  
  
  for (size_t i=0; i < neg_bt.size(); ++i) {
    lp_bb = log_sum_log(log_1_a[i], log_a[i] + neg_bt[i]);
    lp_bf = log_sub_log(log_a[i], log_a[i] + neg_bt[i]);
    lp_fb = log_sub_log(log_1_a[i], log_1_a[i] + neg_bt[i]);
    lp_ff = log_sum_log(log_a[i], log_1_a[i] + neg_bt[i]);
    
    val = log_sum_log(val, r[0][i] + log(lp_bb));
    val = log_sum_log(val, r[1][i] + log(lp_bf));
    val = log_sum_log(val, r[2][i] + log(lp_fb));
    val = log_sum_log(val, r[3][i] + log(lp_ff));
  }
  return val;
}



double ExpTransEstimator::llh_grad_a(const matrix &lr,
                                     const vector<double> &neg_bt) const {
  double val = 0;
  for (size_t i = 0; i < neg_bt.size(); ++i) {
    val += lr[0][i]*( (-1+exp(neg_bt[i])) / (1-a+a*exp(neg_bt[i])) )
           + lr[1][i]*( (1) / (a) )
           + lr[2][i]*( (1) / (1-a) )
           + lr[3][i]*( (1-exp(neg_bt[i])) / (a+(1-a)*exp(neg_bt[i])) );
  }
  return val;
}


double ExpTransEstimator::llh_grad_b(const matrix &lr,
                                     const vector<double> &neg_bt) const {
  double val = 0;
  for (size_t i = 0; i < neg_bt.size(); ++i) {
    val += lr[0][i]*( (-a*b*exp(neg_bt[i])) / (1-a+a*exp(neg_bt[i])) )
           + (lr[1][i] + lr[2][i])*( (b*exp(neg_bt[i])) / (1-exp(neg_bt[i])) )
           + lr[3][i]*( (-(1-a)*exp(b*neg_bt[i])) / (a+(1-a)*exp(neg_bt[i])) );
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
  
  vector<double> log_a (t.size(), 0);
  vector<double> log_1_a (t.size(), 0);
  vector<double> neg_bt (t.size(), 0);
  calc_internal_data(t, new_a, new_b, log_a, log_1_a, neg_bt);
  double moving_llh = calc_log_llh(r, log_a, log_1_a, neg_bt);
  
  std::cout << " start !!! " << " a: " << a << " b: "
            << b << " LLH: " << old_llh << endl;
  
  std::cout << " searching ... " << " a: " << new_a << " b: "
            << new_b << " LLH: " << moving_llh << endl;
  
  while (moving_llh <= old_llh && abs(moving_llh - old_llh) > tolerance) {
    try_step = try_step / 2;
    new_a = a + try_step * grad_a;
    new_b = b + try_step * grad_b;
    
    calc_internal_data(t, new_a, new_b, log_a, log_1_a, neg_bt);
    moving_llh = calc_log_llh(r, log_a, log_1_a, neg_bt);
    std::cout << " searching ... " << " a: " << new_a << " b: "
              << new_b << " LLH: " << moving_llh << endl;
  }
  
  std::cout << " final !!! " << " a: " << new_a << " b: "
            << new_b << " LLH: " << moving_llh << endl;
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
  
  cerr << "[ENTER MLE GRADASCENT]" << endl;
  
  matrix lr;
  lr.resize(4);
  for (size_t i = 0; i < 4; ++i) {
    lr[i].resize(r[i].size());
    for (size_t j = 0; j < r[i].size(); ++j) {
      lr[i][j] = exp(r[i][j]);
    }
  }
  
  vector<double> log_a (t.size(), 0);
  vector<double> log_1_a (t.size(), 0);
  vector<double> neg_bt (t.size(), 0);
  calc_internal_data(t, a, b, log_a, log_1_a, neg_bt);
  double curr_llh = calc_log_llh(r, log_a, log_1_a, neg_bt);
  
  double grad_a, grad_b;
  grad_a = llh_grad_a(lr, neg_bt);
  grad_b = llh_grad_b(lr, neg_bt);
  
  double new_llh;
  GA_stepforward(grad_a, grad_b, curr_llh, new_llh, r, t);
  
  while (abs(new_llh - curr_llh) > tolerance) {
    curr_llh = new_llh;
    calc_internal_data(t, a, b, log_a, log_1_a, neg_bt);
    grad_a = llh_grad_a(lr, neg_bt);
    grad_b = llh_grad_b(lr, neg_bt);
     GA_stepforward(grad_a, grad_b, curr_llh, new_llh, r, t);
  }
}


