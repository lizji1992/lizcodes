/*
  Copyright (C) 2015-2016 University of Southern California
  Authors: Andrew D. Smith and Xiaojing Ji

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

#include "2DCTHMM.hpp"

#include <algorithm>
#include <iomanip>
#include <numeric>
#include <limits>
#include <cmath>

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

typedef vector< vector<double> > matrix;


////////////////////////////////////////////////////////////////////////

template <class T> void
select_vector_elements(const vector<T> &fullset, vector<T> &subset,
                       const vector<size_t> &index) {
  for (size_t i = 0; i < index.size(); ++i) {
    subset.push_back(fullset[index[i]]);
  }
}

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
  a = bg_rate/b;
  
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
    
    const double ff = a + (1 - a) * exp(-(b * time[k]));
    const double bb = 1 - a + a * exp(-(b * time[k]));
    
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
                                  matrix &r, const matrix &ltp) const {
  
  for (size_t i = 0; i < te[0].size() - 1; ++i) {
    const size_t j = i + 1;
    te[0][i] = forward[0][i] + ltp[0][i] + bg_emission(meth[j])
               + backward[0][j] - total;
    double denom = te[0][i];
    
    te[1][i] = forward[0][i] + ltp[1][i] + fg_emission(meth[j])
               + backward[1][j] - total;
    denom = log_sum_log(denom, te[1][i]);
    
    te[2][i] = forward[1][i] + ltp[2][i] + bg_emission(meth[j])
               + backward[0][j] - total;
    denom = log_sum_log(denom, te[2][i]);

    
    te[3][i] = forward[1][i] + ltp[3][i] + fg_emission(meth[j])
               + backward[1][j] - total;
    denom = log_sum_log(denom, te[3][i]);

    r[0][i] = exp(te[0][i] - denom);
    r[1][i] = exp(te[1][i] - denom);
    r[2][i] = exp(te[2][i] - denom);
    r[3][i] = exp(te[3][i] - denom);

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
TwoVarHMM::update_imputed_methylv(vector<pair<double, double> > &meth,
                                  const vector<double> &fg_probs,
                                  const vector<double> &bg_probs) const {
  for (size_t i = 0; i < fg_probs.size(); ++i) {
    if (meth[i].second < 0) {
      meth[i].first = (bg_probs[i] * meth[i].first) /
      (bg_probs[i] * meth[i].first + fg_probs[i] * (1 - meth[i].first));
    }
  }
}



void
TwoVarHMM::estimate_emissions(vector<pair<double, double> > &meth,
                              const vector<double> &meth_lp,
                              const vector<double> &unmeth_lp) {
  
  vector<double> fg_probs(meth_lp.size(), 0);
  vector<double> bg_probs(unmeth_lp.size(), 0);
  
  for (size_t i = 0; i < forward[0].size(); ++i) {
    const double bg = (forward[0][i] + backward[0][i]);
    const double fg = (forward[1][i] + backward[1][i]);
    const double denom = log_sum_log(fg, bg);
    bg_probs[i] = exp(bg - denom);
    fg_probs[i] = exp(fg - denom);
  }
  
  fg_emission.fit(meth_lp, unmeth_lp, fg_probs);
  bg_emission.fit(meth_lp, unmeth_lp, bg_probs);
  
}



double
TwoVarHMM::single_iteration(vector<pair<double, double> > &meth,
                            const vector<size_t> &time,
                            const vector<double> &meth_lp,
                            const vector<double> &unmeth_lp,
                            const size_t curr_itr) {
  
  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ft = log(p_ft);
  const double lp_bt = log(p_bt);
  
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
  estimate_emissions(meth, meth_lp, unmeth_lp);
  
  
  // update transition parameters
  matrix te(4, vector<double> (meth.size(), 0));
  matrix r(4, vector<double> (meth.size(), 0));
  update_trans_estimator(meth, forward_score, te, r, ltp);
  
  if (!NO_RATE_EST) {
    ExpTransEstimator estimate_trans(a, b, BB);
    estimate_trans.mle_GradAscent(r, time);
    double new_a = estimate_trans.get_a();
    double new_b = estimate_trans.get_b();
    double new_bg_rate = new_a * new_b;
    double new_fg_rate = new_b - new_bg_rate;
    double lim = max_iterations;
    if ((fabs(new_fg_rate-fg_rate)/std::min(new_fg_rate, fg_rate) < lim &&
        fabs(new_bg_rate-bg_rate)/std::min(new_bg_rate, bg_rate) < lim)) {
      a = new_a;
      b = new_b;
      fg_rate = new_fg_rate;
      bg_rate = new_bg_rate;
    }
  }

  update_endprob(te);
  
  return total_score;
}



double
TwoVarHMM::BaumWelchTraining(vector<pair<double, double> > &meth,
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
  
  size_t data_size = meth.size();
  forward.resize(2);
  backward.resize(2);
  for (size_t i = 0; i < 2; ++i) {
    forward[i].resize(data_size, 0);
    backward[i].resize(data_size, 0);
  }
  
  ltp.resize(4);
  for (int i = 0; i < 4; ++i) {
    ltp[i].resize(data_size - 1, 0);
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
    
    double score = single_iteration(meth, time, meth_lp, unmeth_lp,
                                    i+1);
    
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
    
    if ((score - prev_score) < tolerance && score >= prev_score) {
      if (VERBOSE)
        cerr << "CONVERGED" << endl << endl;
      break;
    }
    prev_score = score;
    
  }
  
  return prev_score;
}


double
TwoVarHMM::PosteriorDecoding(vector<pair<double, double> > &meth,
                             const vector<size_t> &time,
                             vector<int> &classes, vector<double> &llr_scores,
                             bool IMPUT){
  
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
  
  
  vector<double> fg_probs(data_size, 0);
  vector<double> bg_probs(data_size, 0);
  
  for (size_t i = 0; i < data_size; ++i) {
    
    const double bscore = forward[0][i] + backward[0][i];
    const double fscore = forward[1][i] + backward[1][i];
    const double denom = log_sum_log(bscore, fscore);
    bg_probs[i] = exp(bscore - denom);
    fg_probs[i] = exp(fscore - denom);
    
    llr_scores[i] = fg_probs[i];

    if (bscore > fscore) {
      classes[i] = 0;
      //if (meth[i].second >= 0) { // covered sites
      //  mean_bg_meth += meth[i].first/(meth[i].first + meth[i].second);
      //  n_bg_cpg++;
      //}
    } else {
      classes[i] = 1;
      //if (meth[i].second >= 0) { // covered sites
      //  mean_fg_meth += meth[i].first/(meth[i].first + meth[i].second);
      //  n_fg_cpg++;
      //}
    }
  }
  
  double mean_fg_meth = fg_emission.alpha / (fg_emission.alpha + fg_emission.beta);
  double mean_bg_meth = bg_emission.alpha / (bg_emission.alpha + bg_emission.beta);
 
  if (IMPUT) {
    for (size_t i = 0; i < data_size; ++i) {
      if (meth[i].second < 0) { // not-covered sites
        meth[i].first = mean_fg_meth*fg_probs[i] + mean_bg_meth*bg_probs[i];
      }
    }
  }
  return total_score;
}

