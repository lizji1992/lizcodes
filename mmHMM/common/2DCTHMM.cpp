/*
  Copyright (C) 2008 Cold Spring Harbor Laboratory
  Authors: Andrew D. Smith

  This file is part of methpipe.

  methpipe is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  methpipe is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with methpipe; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "2DCTHMM.hpp"

#include <algorithm>
#include <iomanip>
#include <numeric>
#include <limits>
#include <cmath>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

// #pragma omp <rest of pragma>

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

typedef vector< vector<double> > matrix;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void
TwoVarHMM::set_parameters(const BetaBin _fg_emission,
                          const BetaBin _bg_emission,
                          const double _fg_rate, const double _bg_rate,
                          const double _p_sf, const double _p_sb,
                          const double _p_ft, const double _p_bt) {

  fg_emission = _fg_emission;
  bg_emission = _bg_emission;
  
  fg_rate = _fg_rate;
  bg_rate = _bg_rate;
 
  b = fg_rate + bg_rate;
  a = fg_rate/b;
  
  p_sf = _p_sf;
  p_sb = _p_sb;
  p_ft = _p_ft;
  p_bt = _p_bt;
  
}


inline double
TwoVarHMM::log_sum_log(const double p, const double q) const {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  if (!isfinite(p) && !isfinite(q)) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}


double
TwoVarHMM::log_sum_log_vec(const vector<double> &vals, size_t limit) const {
  const vector<double>::const_iterator x = 
    std::max_element(vals.begin(), vals.begin() + limit);
  const double max_val = *x;
  if (!isfinite(max_val)) {return vals[0];}
  const size_t max_idx = x - vals.begin();
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
#ifdef DEBUG
      assert(isfinite(sum));
#endif
    }
  }
  return max_val + log(sum);
}


double
TwoVarHMM::forward_algorithm(const vector<pair<double, double> > &meth,
                             const vector<size_t> &time,
                             const double lp_sf, const double lp_sb,
                             const double lp_ft, const double lp_bt,
                             matrix &ltp) {
  
  const size_t end = forward[0].size();
  
  forward[0][0] = bg_emission(meth[0]) + lp_sb; // background
  forward[1][0] = fg_emission(meth[0]) + lp_sf; // foreground
  
  for (size_t i = 1; i < end; ++i) {
    const size_t k = i - 1;
    const size_t dist = time[i] - time[k];
    
    const double ff = a + (1 - a) * exp(-(b * dist));
    const double bb = 1 - a + a * exp(-(b * dist));
    
    const double lp_ff = log(ff);
    const double lp_fb = log(1 - ff);
    const double lp_bf = log(1 - bb);
    const double lp_bb = log(bb);
    ltp[0][k] = lp_bb;
    ltp[1][k] = lp_bf;
    ltp[2][k] = lp_fb;
    ltp[3][k] = lp_ff;
    
    // background
    forward[0][i] = (bg_emission(meth[i]) +
                     log_sum_log(forward[0][k] + lp_bb, forward[1][k] + lp_fb));
    // foreground
    forward[1][i] = (fg_emission(meth[i]) +
                     log_sum_log(forward[0][k] + lp_bf, forward[1][k] + lp_ff));
  }
  
  return log_sum_log(forward[0][end - 1] + lp_bt, forward[1][end - 1] + lp_ft);
}



double
TwoVarHMM::backward_algorithm(const vector<pair<double, double> > &meth,
                              const vector<size_t> &time,
                              const double lp_sf, const double lp_sb,
                              const double lp_ft, const double lp_bt,
                              const matrix &ltp) {
  
  const size_t end = backward[0].size();
  
  backward[0][end - 1] = lp_bt; // background
  backward[1][end - 1] = lp_ft; // foreground
  
  
  for (size_t k = end - 1; k > 0; --k) {
    const size_t i = k - 1;
    
    const double lp_bb = ltp[0][i];
    const double lp_bf = ltp[1][i];
    const double lp_fb = ltp[2][i];
    const double lp_ff = ltp[3][i];
    
    const double bg_emi = bg_emission(meth[k]) + backward[0][k];
    const double fg_emi = fg_emission(meth[k]) + backward[1][k];
    
    // background
    backward[0][i] = log_sum_log(bg_emi + lp_bb, fg_emi + lp_bf);
    // foreground
    backward[1][i] = log_sum_log(bg_emi + lp_fb, fg_emi + lp_ff);
  }
  
  return log_sum_log(backward[0][0] + bg_emission(meth[0]) + lp_sb,
                     backward[1][0] + fg_emission(meth[0]) + lp_sf);
}


void
TwoVarHMM::update_trans_estimator(const vector<pair<double, double> > &meth,
                                  const double total, matrix &te,
                                  const matrix &ltp) const {
  
  for (size_t i = 0; i < te[0].size() - 1; ++i) {
    const size_t j = i + 1;
    te[0][i] = forward[0][i] + ltp[0][i] + bg_emission(meth[j])
               + backward[0][j] - total;
    te[1][i] = forward[0][i] + ltp[1][i] + fg_emission(meth[j])
               + backward[1][j] - total;
    te[2][i] = forward[1][i] + ltp[2][i] + bg_emission(meth[j])
               + backward[0][j] - total;
    te[3][i] = forward[1][i] + ltp[3][i] + fg_emission(meth[j])
               + backward[1][j] - total;
   // std::cout << te[0][i] << " " << te[1][i]
   //           << " " << te[2][i] << " " << te[3][i] << endl;
   }
}



void
TwoVarHMM::update_endprob(matrix &te) {
  size_t T = te[0].size();
  
  const double p_sum_bb = exp(log_sum_log_vec(te[0], T));
  const double p_sum_bf = exp(log_sum_log_vec(te[1], T));
  const double p_sum_fb = exp(log_sum_log_vec(te[2], T));
  const double p_sum_ff = exp(log_sum_log_vec(te[3], T));

  
  double fdenom = (p_sum_ff + p_sum_fb);
  double p_average_ff = max(p_sum_ff/fdenom - p_ft/2.0, MIN_PROB);
  double p_average_fb = max(p_sum_fb/fdenom - p_ft/2.0, MIN_PROB);
  
  
  double bdenom = (p_sum_bf + p_sum_bb);
  double p_average_bf = max(p_sum_bf/bdenom - p_bt/2.0, MIN_PROB);
  double p_average_bb = max(p_sum_bb/bdenom - p_bt/2.0, MIN_PROB);
  
  p_sb = (p_average_bb + p_average_fb)/2.0;
  p_sf = (p_average_bf + p_average_ff)/2.0;
  
}


void
TwoVarHMM::estimate_emissions(const vector<double> &meth_lp,
                              const vector<double> &unmeth_lp) {
  
  vector<double> fg_probs(meth_lp.size());
  vector<double> bg_probs(meth_lp.size());
  
  for (size_t i = 0; i < forward[0].size(); ++i) {
    const double bg = (forward[0][i] + backward[0][i]);
    const double fg = (forward[1][i] + backward[1][i]);
    const double denom = log_sum_log(fg, bg);

    fg_probs[i] = exp(fg - denom);
    bg_probs[i] = exp(bg - denom);
  }
  
  fg_emission.fit(meth_lp, unmeth_lp, fg_probs);
  bg_emission.fit(meth_lp, unmeth_lp, bg_probs);
}

  

double
TwoVarHMM::single_iteration(const vector<pair<double, double> > &meth,
                            const vector<size_t> &time,
                            const vector<double> &meth_lp,
                            const vector<double> &unmeth_lp) {
  
  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ft = log(p_ft);
  const double lp_bt = log(p_bt);
  
  
  // prepare internal data structures
  size_t data_size = meth.size();
  forward.resize(2);
  backward.resize(2);
  for (size_t i = 0; i < 2; ++i) {
    forward[i].resize(data_size);
    backward[i].resize(data_size);
  }
  
  ltp.resize(4);
  for (int i = 0; i < 4; ++i) {
    ltp[i].resize(data_size - 1);
  }
  
  
  // forward/backward algorithm
  const double forward_score =
      forward_algorithm(meth, time, lp_sf, lp_sb, lp_ft, lp_bt, ltp);
  
  const double backward_score =
      backward_algorithm(meth, time, lp_sf, lp_sb, lp_ft, lp_bt, ltp);
  
    
  if (DEBUG && (fabs(forward_score - backward_score) /
                max(forward_score, backward_score)) > 1e-10)
    cerr << "fabs(forward_score - backward_score)/"
         << "max(forward_score, backward_score) > 1e-10" << endl;
  
  double
  total_score = forward_score;

  // update emission
  estimate_emissions(meth_lp, unmeth_lp);
  
  
  // update transition parameters
  matrix te(4, vector<double> (meth.size(), 0));
  update_trans_estimator(meth, forward_score, te, ltp);
  
  ExpTransEstimator estimate_trans(a, b);
  estimate_trans.mle_GradAscent(te, time);
  a = estimate_trans.get_a();
  b = estimate_trans.get_b();
  fg_rate = a * b;
  bg_rate = b - fg_rate;
  update_endprob(te);
  
  return total_score;
}



double
TwoVarHMM::BaumWelchTraining(const vector<pair<double, double> > &meth,
                             const vector<size_t> &time) {
  
  cerr << "[ENTER BM-TRAINING]" << endl;
  
  vector<double> meth_lp(meth.size()), unmeth_lp(meth.size());
  for (size_t i = 0; i < meth.size(); ++i) {
    meth_lp[i] =
    log(min(max(meth[i].first/(meth[i].first + meth[i].second), 1e-2),
            1.0 - 1e-2));
    unmeth_lp[i] =
    log(1 - min(max(meth[i].first/(meth[i].first + meth[i].second), 1e-2),
                1.0 - 1e-2));
  }
  
  
  if (VERBOSE)
    cerr << setw(5)  << "ITR"
    << setw(18) << "F RATE"
    << setw(18) << "B RATE"
    << setw(18) << "F PARAMS"
    << setw(18) << "B PARAMS"
    << setw(14) << "DELTA"
    << endl;
  
  // observation likelihood
  double prev_score = -std::numeric_limits<double>::max();
  
  cerr << setw(5) << "0"
  << setw(18) << fg_rate
  << setw(18) << bg_rate
  << setw(18) << fg_emission.tostring()
  << setw(18) << bg_emission.tostring()
  << setw(14) << "NA"
  << setw(14) << "NA"
  << setw(14) << "NA"
  << endl;
  
  
  for(size_t i = 0; i < max_iterations; ++i) {
    
    double score = single_iteration(meth, time, meth_lp, unmeth_lp);
    
    if (VERBOSE) {
      cerr << setw(5) << i + 1
      << setw(18) << fg_rate
      << setw(18) << bg_rate
      << setw(18) << fg_emission.tostring()
      << setw(18) << bg_emission.tostring()
      << setw(14) << score
      << setw(14) << prev_score
      << setw(14) << (score - prev_score)/std::fabs(score)
      << endl;
    }
    
    if ((score - prev_score) < tolerance) {
      if (VERBOSE)
        cerr << "CONVERGED" << endl << endl;
      break;
    }
    prev_score = score;
    
  }
  
  return prev_score;
}


double
TwoVarHMM::PosteriorDecoding(const vector<pair<double, double> > &meth,
                             const vector<size_t> &time,
                             vector<int> &classes, vector<double> &llr_scores){
  

  size_t data_size = meth.size();
  forward.resize(2);
  backward.resize(2);
  for (size_t i = 0; i < 2; ++i) {
    forward[i].resize(data_size);
    backward[i].resize(data_size);
  }
  
  ltp.resize(4);
  for (int i = 0; i < 4; ++i) {
    ltp[i].resize(data_size - 1);
  }
  
  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ft = log(p_ft);
  const double lp_bt = log(p_bt);
  
  const double forward_score =
      forward_algorithm(meth, time, lp_sf, lp_sb, lp_ft, lp_bt, ltp);
  
  const double backward_score =
      backward_algorithm(meth, time, lp_sf, lp_sb, lp_ft, lp_bt, ltp);
  
  if (DEBUG && (fabs(forward_score - backward_score) /
                max(forward_score, backward_score)) > 1e-10)
    cerr << "fabs(forward_score - backward_score)/"
         << "max(forward_score, backward_score) > 1e-10" << endl;
    
  double  total_score = forward_score;
  
  classes.resize(data_size);
  llr_scores.resize(data_size);
  
  for (size_t i = 0; i < data_size; ++i) {
    
    double bscore = forward[0][i] + backward[0][i];
    double fscore = forward[1][i] + backward[1][i];
    double total_state_score = log_sum_log(bscore, fscore);
    
    if (bscore > fscore) {
      classes[i] = 0;
      llr_scores[i] = exp(bscore - total_state_score);
    } else {
      classes[i] = 1;
      llr_scores[i] = exp(fscore - total_state_score);
    }
  }
  
  return total_score;
}
