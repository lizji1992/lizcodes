/*
 Copyright (C) 2020-2021 University of Southern California
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>


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
           + r[2][i]*( (1) / (a-1) )
           + r[3][i]*( (1-ebt[i]) / (a+(1-a)*ebt[i]) );
  }
  return val;
}


double ExpTransEstimator::llh_grad_b(const matrix &r,
                                     const vector<double> &ebt,
                                     const vector<size_t> &t) const {
  double val = 0;
  for (size_t i = 0; i < ebt.size(); ++i) {
    val += r[0][i]*( (-a*t[i]*ebt[i]) / (1-a+a*ebt[i]) )
           + (r[1][i] + r[2][i]) * ( (t[i]*ebt[i]) / (1-ebt[i]) )
           + r[3][i]*( (-(1-a)*t[i]*ebt[i]) / (a+(1-a)*ebt[i]) );
  }
  return val;
}


void ExpTransEstimator::GA_stepforward(const double grad_a,
                                       const double grad_b,
                                       const double &old_llh, double &new_llh,
                                       const matrix &r,
                                       const vector<size_t> &t)  {

  double try_step = step_size;

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
    try_step = try_step / 10;
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


void
ExpTransEstimator::GA_stepforward_BB(const double grad_a, const double grad_b,
                                     const double d_grad_a,
                                     const double d_grad_b,
                                     double &d_a, double &d_b,
                                     const double &old_llh, double &new_llh,
                                     const matrix &r,
                                     const vector<size_t> &t) {

  double try_step = fabs(d_grad_a * d_a + d_grad_b * d_b) /
  (d_grad_a * d_grad_a + d_grad_b * d_grad_b);

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
    try_step = try_step / 10;
    new_a = a + try_step * grad_a;
    new_b = b + try_step * grad_b;

    if (new_a > 0 && new_a < 1 && new_b > 0) {
      calc_internal_data(t, new_b, ebt);
      moving_llh = calc_llh(r, new_a, ebt);
    }
    ++itr;
  }

  if (moving_llh > old_llh && new_a > 0 && new_a < 1 && new_b > 0) {
    d_a = new_a - a;
    d_b = new_b - b;
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

  double a_old = a, b_old = b;
  double grad_a, grad_b;
  grad_a = llh_grad_a(r, ebt);
  grad_b = llh_grad_b(r, ebt, t);

  // first iteration
  double new_llh;
  GA_stepforward(grad_a, grad_b, curr_llh, new_llh, r, t);

  double d_a = a - a_old;
  double d_b = b - b_old;

  double new_grad_a, new_grad_b;
  double d_grad_a, d_grad_b;
  size_t itr = 1;

  while (abs(new_llh - curr_llh) > tolerance && itr < max_iteration) {
    curr_llh = new_llh;
    calc_internal_data(t, b, ebt);
    if (BB) {
      new_grad_a = llh_grad_a(r, ebt);
      new_grad_b = llh_grad_b(r, ebt, t);
      d_grad_a = new_grad_a - grad_a;
      d_grad_b = new_grad_b - grad_b;
      GA_stepforward_BB(new_grad_a, new_grad_b, d_grad_a, d_grad_b,
                        d_a, d_b, curr_llh, new_llh, r, t);
      grad_a = new_grad_a;
      grad_b = new_grad_b;
    } else {
      grad_a = llh_grad_a(r, ebt);
      grad_b = llh_grad_b(r, ebt, t);
      GA_stepforward(grad_a, grad_b, curr_llh, new_llh, r, t);
    }
    ++itr;
  }
  cerr << "Used " << itr << " iterations." << endl;
}

// CONJUGATE GRADIENT METHOD

double
f(const gsl_vector *x, void *params) {
  double *par = (double *) params;
  size_t n = par[0];
  double u = gsl_vector_get(x, 0);
  double v = gsl_vector_get(x, 1);

  double val = 0;
  for (size_t i = 0; i < n; ++i) {
    size_t s = i * 5;
    const double ebt = exp(- v * par[s+5]);
    const double d = par[s+1]*log(1 - u + u * ebt)
                     + par[s+2]*log(u - u * ebt)
                     + par[s+3]*log(1 - u - (1 - u) * ebt)
                     + par[s+4]*log(u + (1 - u) * ebt);
    val -= d;
  }
  return val;
}

void
df(const gsl_vector *x, void *params, gsl_vector *d) {
  double *par = (double *) params;
  size_t n = par[0];
  double u = gsl_vector_get(x, 0);
  double v = gsl_vector_get(x, 1);

  double da = 0;
  double db = 0;
  for (size_t i = 0; i < n; ++i) {
    const size_t s = i * 5;
    const double t = par[s+5];
    const double ebt = exp(- v * t);
    da -= par[s+1]*( (-1+ebt) / (1-u+u*ebt) )
          + par[s+2]*( (1) / (u) )
          + par[s+3]*( (1) / (u-1) )
          + par[s+4]*( (1-ebt) / (u+(1-u)*ebt) );

    db -= par[s+1]*( (-u*t*ebt) / (1-u+u*ebt) )
          + (par[s+2] + par[s+3]) * ( (t*ebt) / (1-ebt) )
          + par[s+4]*( (-(1-u)*t*ebt) / (u+(1-u)*ebt) );
  }
  //cerr << "da: " << da  << ", db: " << db << endl;
  gsl_vector_set(d, 0, da);
  gsl_vector_set(d, 1, db);
}

void
fdf(const gsl_vector *x, void *params, double *fp,
                       gsl_vector *d) {
  *fp = f(x, params);
  df(x, params, d);
}


void ExpTransEstimator::mle_CG(const matrix &r, const vector<size_t> &t) {

  // build up parameter list
  size_t n = t.size();
  double params[5*n+1];

  params[0] = n;  // the first parameter is the length of sequence
  for (size_t i = 0; i < n; ++i) {
    const size_t s = i * 5;
    params[s+1] = r[0][i];
    params[s+2] = r[1][i];
    params[s+3] = r[2][i];
    params[s+4] = r[3][i];
    params[s+5] = t[i];
  }

  gsl_multimin_function_fdf func;
  func.f = &f;
  func.df = &df;
  func.fdf = &fdf;
  func.n = 2;
  func.params = &params;

  // set starting point
  gsl_vector *x = gsl_vector_alloc(2);
  gsl_vector_set(x, 0, a);
  gsl_vector_set(x, 1, b);

  // allocate and set the minimizer and its type
  const gsl_multimin_fdfminimizer_type *type =
    gsl_multimin_fdfminimizer_conjugate_fr;
  gsl_multimin_fdfminimizer *minimizer =
    gsl_multimin_fdfminimizer_alloc (type, 2);

  // set tolerance and starting step size
  gsl_multimin_fdfminimizer_set(minimizer, &func, x, step_size, tolerance);

  size_t itr = 0;
  int status = 0;

  itr++;
  status = gsl_multimin_fdfminimizer_iterate(minimizer);
  // check for convergence
  status = gsl_multimin_test_gradient(minimizer->gradient, tolerance);
  double last_legal_a = a, last_legal_b = b;
  double new_a = gsl_vector_get(minimizer->x, 0);
  double new_b = gsl_vector_get(minimizer->x, 1);
  if (new_a > 0 && new_a < 1 && new_b > 0) {
    last_legal_a = new_a;
    last_legal_b = new_b;
  }

  while (status == GSL_CONTINUE && itr < max_iteration &&
         new_a > 0 && new_a < 1 && new_b > 0) {
    itr++;
    status = gsl_multimin_fdfminimizer_iterate(minimizer);

    if (status) break;
    status = gsl_multimin_test_gradient(minimizer->gradient, tolerance);

    new_a = gsl_vector_get(minimizer->x, 0);
    new_b = gsl_vector_get(minimizer->x, 1);
    if (new_a > 0 && new_a < 1 && new_b > 0) {
      last_legal_a = new_a;
      last_legal_b = new_b;
    }
    //cerr << "new a: " << new_a << " new b: " << new_b << endl;
  }
  cerr << "Used " << itr << " iterations." << endl;

  a = last_legal_a;
  b = last_legal_b;
  // free the memory
  gsl_multimin_fdfminimizer_free(minimizer);
  gsl_vector_free(x);
}
