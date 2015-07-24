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

#include "NBVDHMM_BF.hpp"

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
TwoVarHMM::update_emission_matrix() {
  emission = vector<BetaBin> (num_states);
  for (size_t i = 0; i < fg_mode; ++i) {
    emission[i] = fg_emission;
  }
  for (size_t i = fg_mode; i < num_states; ++i) {
    emission[i] = bg_emission;
  }
}

void
TwoVarHMM::update_transition_matrix() {
  trans = matrix(num_states, vector<double>(num_states, 0));
  for (size_t i = 0; i < fg_mode; ++i) {
    trans[i][i] = 1 - fg_p;
    trans[i][i+1] = fg_p;
  }
  for (size_t j = 0; j < bg_mode; ++j) {
    trans[fg_mode+j][fg_mode+j] = 1 - bg_p;
    size_t bf = j == (bg_mode - 1) ? 0 : fg_mode+j+1;
    trans[fg_mode+j][bf] = bg_p;
  }
  
  start_trans = vector<double>(num_states, 0);
  start_trans[0] = 0.5;
  start_trans[fg_mode] = 0.5;
  
  end_trans = vector<double>(num_states, 1e-10);
}



void
TwoVarHMM::set_parameters(const BetaBin _fg_emission,
                          const BetaBin _bg_emission,
                          const size_t _fg_mode, const size_t _bg_mode,
                          const double _fg_p, const double _bg_p,
                          const double _p_sf, const double _p_sb,
                          const double _p_ft, const double _p_bt) {

  fg_emission = _fg_emission;
  bg_emission = _bg_emission;
  
  fg_mode = _fg_mode;
  bg_mode = _bg_mode;
  num_states = fg_mode + bg_mode;
  
  fg_p = _fg_p;
  bg_p = _bg_p;
  p_sf = _p_sf;
  p_sb = _p_sb;
  p_ft = _p_ft;
  p_bt = _p_bt;
  
  update_emission_matrix();
  update_transition_matrix();
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
				const size_t start, const size_t end,
        const vector<double> &lp_s, const vector<double> &lp_t,
        const vector< vector<double> > &lp_trans) {
  
  //cerr << "[ENTER FORWARD COMPUTING]" << endl;
  //cerr << start << ", " << end << endl;
  
  for (size_t i = 0; i < num_states; ++i) {
    forward[i][start] = emission[i](meth[start]) + lp_s[i];
  }

  for (size_t i = start + 1; i < end; ++i) {
    const size_t k = i - 1;
    // calculate forward score
    for (size_t s2 = 0; s2 < num_states; ++s2) {
      forward[s2][i] = forward[0][k] + lp_trans[0][s2] + emission[s2](meth[i]);
      for (size_t s1 = 1; s1 < num_states; ++s1) {
        forward[s2][i] = log_sum_log(forward[s2][i],
                                     forward[s1][k] + lp_trans[s1][s2] + emission[s2](meth[i]));
//        cerr << forward[s2][i] << " ";
      }
//      forward[s2][i] += emission[s2](meth[i]);
    }
  }
  double forward_score = forward[0][end - 1] + lp_t[0];
  for (size_t s = 1; s < num_states; ++s) {
    forward_score = log_sum_log(forward_score, forward[s][end-1] + lp_t[s]);
  }
  
  return forward_score;
}


double
TwoVarHMM::backward_algorithm(const vector<pair<double, double> > &meth,
				const size_t start, const size_t end,
				const vector<double> &lp_s, const vector<double> &lp_t,
				const vector< vector<double> > &lp_trans) {
  
  //cerr << "[ENTER BACKWARD COMPUTING]" << endl;
  //cerr << start << ", " << end << ", " << backward.size() << ", "
  //<< backward[0].size() << endl;
  for (size_t i = 0; i < num_states; ++i) {
    backward[i][end - 1] = lp_t[i];
  }
  
  for (size_t k = end - 1; k > start; --k) {
    size_t i = k - 1;
    // calculate backward score
    for (size_t s1 = 0; s1 < num_states; ++s1) {
      backward[s1][i] = backward[0][k] + lp_trans[s1][0] + emission[0](meth[k]);
      for (size_t s2 = 1; s2 < num_states; ++s2) {
        //<< "log "<< backward[s1][i] << "+" 
        //<< " log " << backward[s2][k]  << "+" << lp_trans[s1][s2] << "+" << emission[s2](meth[k]);
        backward[s1][i] = log_sum_log(backward[s1][i],
                                      backward[s2][k] + lp_trans[s1][s2] +
                                      emission[s2](meth[k]));
        //std::cout << "s1:" << s1 << " s2:" << i << backward[s1][i] << endl;
        //std::cout << "=" << backward[s1][i] << endl;
        //if (!isfinite(lp_trans[s1][s2]))
      }
    }
  }
 //cerr << endl;
  double backward_score = backward[0][start] + emission[0](meth[start])+lp_s[0];

  for (size_t s = 1; s < num_states; ++s) {
    backward_score = log_sum_log(backward_score, backward[s][start] +
                                 emission[s](meth[start])+lp_s[s]);
  }

  return backward_score;
}



void
TwoVarHMM::update_trans_estimator(const vector<pair<double, double> > &meth,
                                const size_t start, const size_t end,
                                const double total, vector<matrix> &te,
                                const vector<vector<double> > &lp_trans) const {
  
  //cerr << "[UPDATE TRANS ESTIMATOR]" << endl;
  for (size_t i = start + 1; i < end; ++i) {
    const size_t k = i - 1;
    for(size_t x = 0; x < num_states; ++x) {
      for (size_t y = 0; y < num_states; ++y) {
        te[x][y][k] = forward[x][k] + lp_trans[x][y] + emission[y](meth[i])
                      + backward[y][i] - total;
//        if (y != 0 && (y-x) > 1 && isfinite(te[x][y][k]))
  //        std::cout << "x:" << x << " y:" << y << " k:" << k << " : " << te[x][y][k] << endl;
      }
    }
  }
}



void
TwoVarHMM::update_transitions(const vector<matrix> &te) {
  //cerr << "[UPDATE TRANSITIONS]" << endl;
  size_t T = te[0][0].size();
  
  // Subtracting 1 from the limit of the summation because the final
  // term has no meaning since there is no transition to be counted
  // from the final observation (they all must go to terminal state)
  
  double fg_in = log_sum_log_vec(te[0][0], T - 1);
  double fg_out = log_sum_log_vec(te[0][1], T - 1);
  for(size_t i = 1; i < fg_mode; ++i) {
    fg_in = log_sum_log(fg_in, log_sum_log_vec(te[i][i], T - 1));
    fg_out = log_sum_log(fg_out, log_sum_log_vec(te[i][i+1], T - 1));
  }

  double p_fg_in = exp(fg_in);
  double p_fg_out = exp(fg_out);
  double fg_denom = p_fg_in + p_fg_out;
  fg_p = p_fg_out/fg_denom -
         std::accumulate(end_trans.begin(), end_trans.begin()+fg_mode, 0) / 2;
  
  
  double bg_in = log_sum_log_vec(te[num_states-1][num_states-1], T - 1);
  double bg_out = log_sum_log_vec(te[num_states-1][0], T - 1);
  for(size_t j = num_states-2; j >= fg_mode; ++j) {
    bg_in = log_sum_log(bg_in, log_sum_log_vec(te[j][j], T - 1));
    bg_out = log_sum_log(bg_out, log_sum_log_vec(te[j][j+1], T - 1));
  }
  
  double p_bg_in = exp(bg_in);
  double p_bg_out = exp(bg_out);
  double bg_denom = p_bg_in + p_bg_out;
  bg_p = p_bg_out/bg_denom -
         std::accumulate(end_trans.begin()+fg_mode, end_trans.end(), 0) / 2;
  
  // estimate start probabilities ????
  p_sf = (bg_p + 1 - fg_p)/2.0;
  p_sb = (1 - bg_p + fg_p)/2.0;
  //p_sf = exp(log_sum_log(te[0][0][0], te[0][1][0]));
  
}


void
TwoVarHMM::estimate_emissions(const vector<double> &meth_lp,
                              const vector<double> &unmeth_lp) {
  //cerr << "[ESTIMATE EMISSIONS]" << endl;
  vector<double> fg_probs(meth_lp.size());
  vector<double> bg_probs(meth_lp.size());
  
  for (size_t i = 0; i < forward[0].size(); ++i) {

    double fg = forward[0][i] + backward[0][i];
    double denom = fg;
    for(size_t k = 1; k < fg_mode; ++k) {
      fg = log_sum_log(fg, forward[k][i] + backward[k][i]);
      denom = log_sum_log(denom, fg);
    }
    
    double bg = forward[fg_mode][i] + backward[fg_mode][i];
    denom = log_sum_log(denom, bg);
    for(size_t k = fg_mode+1; k < num_states; ++k) {
      bg = log_sum_log(bg, forward[k][i] + backward[k][i]);
      denom = log_sum_log(denom, bg);
    }
    
    fg_probs[i] = exp(fg - denom);
    bg_probs[i] = exp(bg - denom);
    //cerr << forward[0][i] << " " << backward[0][i] << ", ";
    //cerr << fg_probs[i] << " " << bg_probs[i] << ", ";
  }
  // cerr << endl;
  
  fg_emission.fit(meth_lp, unmeth_lp, fg_probs);
  bg_emission.fit(meth_lp, unmeth_lp, bg_probs);
}

  

double
TwoVarHMM::single_iteration(const vector<pair<double, double> > &meth,
                            const vector<size_t> &reset_points,
                            const vector<double> &meth_lp,
                            const vector<double> &unmeth_lp) {
  
  //std::cout << "[ENTER SINGLE ITERATION]" << endl;
  
  double total_score = 0;
  
  // prepare forward/backward vectors
  size_t data_size = meth.size();
  //forward = matrix(num_states, vector<double>(data_size, 0));
  //backward = matrix(num_states, vector<double>(data_size, 0));
  forward.resize(num_states);
  backward.resize(num_states);
  for (size_t i = 0; i < num_states; ++i) {
    forward[i].resize(data_size);
    backward[i].resize(data_size);
  }
  
  // get log transition probability
  vector<double> lp_start_trans = vector<double>(num_states, 0);
  vector<double> lp_end_trans = vector<double>(num_states, 0);
  vector< vector<double> > lp_trans =
      vector< vector<double> >(num_states, vector<double>(num_states, 0));

  for (size_t i = 0; i < lp_start_trans.size(); ++i) {
    lp_start_trans[i] = log(start_trans[i]);
  }
  for (size_t i = 0; i < lp_end_trans.size(); ++i) {
    lp_end_trans[i] = log(end_trans[i]);
  }
  for (size_t i = 0; i < lp_trans.size(); ++i) {
    lp_trans.resize(trans[i].size());
    for (size_t j = 0; j < trans[i].size(); ++j) {
      lp_trans[i][j] = log(trans[i][j]);
     // std::cout << lp_trans[i][j] << " ";
    }
  //  std::cout << endl;
  }
 // std::cout << endl;
  
  // for estimating transitions
  vector<matrix> te(num_states,
                    matrix(num_states, vector<double> (meth.size(), 0)));
  
  // forward/backward algorithm
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double forward_score =
        forward_algorithm(meth, reset_points[i], reset_points[i + 1],
                          lp_start_trans, lp_end_trans, lp_trans);
    const double backward_score =
        backward_algorithm(meth, reset_points[i], reset_points[i + 1],
                           lp_start_trans, lp_end_trans, lp_trans);
  
    
    if (DEBUG && (fabs(forward_score - backward_score) /
                  max(forward_score, backward_score)) > 1e-10)
      cerr << "fabs(forward_score - backward_score)/"
           << "max(forward_score, backward_score) > 1e-10" << endl;
    
    update_trans_estimator(meth, reset_points[i], reset_points[i + 1],
                           forward_score, te, lp_trans);
    
    total_score += forward_score;
  }

  update_transitions(te);

  estimate_emissions(meth_lp, unmeth_lp);
  
  return total_score;
}


double
TwoVarHMM::compute_likelihood(const vector<pair<double, double> > &meth,
                              const vector<size_t> &reset_points,
                              const vector<double> &meth_lp,
                              const vector<double> &unmeth_lp) {
  
  double total_score = 0;
  
  // prepare forward/backward vectors
  size_t data_size = meth.size();
  forward.resize(num_states);
  backward.resize(num_states);
  for (size_t i = 0; i < num_states; ++i) {
    forward[i].resize(data_size);
    backward[i].resize(data_size);
  }
  
  // get log transition probability
  vector<double> lp_start_trans = vector<double>(num_states, 0);
  vector<double> lp_end_trans = vector<double>(num_states, 0);
  vector< vector<double> > lp_trans =
  vector< vector<double> >(num_states, vector<double>(num_states, 0));
  
  for (size_t i = 0; i < lp_start_trans.size(); ++i) {
    lp_start_trans[i] = log(start_trans[i]);
  }
  for (size_t i = 0; i < lp_end_trans.size(); ++i) {
    lp_end_trans[i] = log(end_trans[i]);
  }
  for (size_t i = 0; i < lp_trans.size(); ++i) {
    lp_trans.resize(trans[i].size());
    for (size_t j = 0; j < trans[i].size(); ++j) {
      lp_trans[i][j] = log(trans[i][j]);
      // std::cout << lp_trans[i][j] << " ";
    }
    //  std::cout << endl;
  }
  // std::cout << endl;

  // forward/backward algorithm
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double forward_score =
    forward_algorithm(meth, reset_points[i], reset_points[i + 1],
                      lp_start_trans, lp_end_trans, lp_trans);
    total_score += forward_score;
  }

  return total_score;
}



double
TwoVarHMM::BaumWelchTraining(const vector<pair<double, double> > &meth,
                             const vector<size_t> &reset_points) {
  
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
    << setw(18) << "F P"
    << setw(18) << "B P"
    << setw(18) << "F PARAMS"
    << setw(18) << "B PARAMS"
    << setw(10)  << "F mode"
    << setw(14) << "DELTA"
    << endl;
  
  
  double prev_score = -std::numeric_limits<double>::max(); // total likelihood

  cerr << setw(5) << "0"
  << setw(18) << fg_p
  << setw(18) << bg_p
  << setw(18) << fg_emission.tostring()
  << setw(18) << bg_emission.tostring()
  << setw(5) << fg_mode
  << setw(14) << "NA"
  << setw(14) << "NA"
  << setw(14) << "NA"
  << endl;
  
  size_t i = 0;
  double score = single_iteration(meth, reset_points, meth_lp, unmeth_lp);
  if (VERBOSE) {
    cerr << setw(5) << i + 1
    << setw(18) << fg_p
    << setw(18) << bg_p
    << setw(18) << fg_emission.tostring()
    << setw(18) << bg_emission.tostring()
    << setw(5) << fg_mode
    << setw(14) << score
    << setw(14) << prev_score
    << setw(14) << (score - prev_score)/std::fabs(score)
    << endl;
  }
  
  prev_score = score;
  
  BetaBin fg_emission_better = fg_emission;
  BetaBin bg_emission_better = bg_emission;
  size_t fg_mode_better = fg_mode;
  double fg_p_better = fg_p;
  double bg_p_better = bg_p;
  double p_sf_better = p_sf;
  double p_sb_better = p_sb;
  double p_ft_better = p_ft;
  double p_bt_better = p_bt;
  
  for(i = 1; i < max_iterations; ++i) {
    // start from estimated transition/emission, alter fg_mode
    //cerr << "NEW ITERATION " << i+1 << endl;
    size_t start_k = 1;
    if (fg_mode > mode_search_k) start_k = fg_mode - mode_search_k;
    for(size_t k = start_k; k <= (fg_mode+mode_search_k); ++k) {
      // start trainning a branch hmm
      TwoVarHMM hmm_branch(mode_search_k, MIN_PROB, tolerance,
                           max_iterations, VERBOSE);
      hmm_branch.set_parameters(fg_emission, bg_emission, k, bg_mode,
                                fg_p, bg_p, p_sf, p_sb, p_ft, p_bt);
      
      hmm_branch.single_iteration(meth, reset_points, meth_lp, unmeth_lp);
      double new_score = hmm_branch.compute_likelihood(meth, reset_points,
                                                       meth_lp, unmeth_lp);
      
      //cerr << "BRANCH: " << k << " NEW SCORE: " << new_score << endl;
      
      if (new_score > score) {
        score = new_score;
        fg_emission_better = hmm_branch.fg_emission;
        bg_emission_better = hmm_branch.bg_emission;
        fg_mode_better = hmm_branch.fg_mode;
        fg_p_better = hmm_branch.fg_p;
        bg_p_better = hmm_branch.bg_p;
        p_sf_better = hmm_branch.p_sf;
        p_sb_better = hmm_branch.p_sb;
        p_ft_better = hmm_branch.p_ft;
        p_bt_better = hmm_branch.p_bt;
      }
    }
    
    //cerr << "SET PARAMETERS FOR ITR " << i+1 << endl;
    set_parameters(fg_emission_better, bg_emission_better,
                   fg_mode_better, bg_mode, fg_p_better, bg_p_better,
                   p_sf_better, p_sb_better, p_ft_better, p_bt_better);
    
    
    if (VERBOSE) {
      cerr << setw(5) << i + 1
      << setw(18) << fg_p
      << setw(18) << bg_p
      << setw(18) << fg_emission.tostring()
      << setw(18) << bg_emission.tostring()
      << setw(5) << fg_mode
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
   // cerr << "END ITERATION " << i+1 << endl;

  }
  
  return score;
}


double
TwoVarHMM::PosteriorDecoding(const vector<pair<double, double> > &meth,
                             const vector<size_t> &reset_points,
                             vector<int> &classes,
                             vector<double> &llr_scores) {
  
  //cerr << "[ENTER POSTERIOR DECODING]" << endl;
  double total_score = 0;
  
  //cerr << fg_p << " ," << bg_p << " ," << p_sf << " ," << p_sb
  //     << " ," << p_ft << " ," << p_bt << endl;
  // prepare forward/backward vectors
  size_t data_size = meth.size();
  forward.resize(num_states);
  backward.resize(num_states);
  for (size_t i = 0; i < num_states; ++i) {
    forward[i].resize(data_size);
    backward[i].resize(data_size);
  }
  
  // get log transition probability
  vector<double> lp_start_trans = vector<double>(num_states, 0);
  vector<double> lp_end_trans = vector<double>(num_states, 0);
  vector< vector<double> > lp_trans =
  vector< vector<double> >(num_states, vector<double>(num_states, 0));
  
  for (size_t i = 0; i < lp_start_trans.size(); ++i) {
    lp_start_trans[i] = log(start_trans[i]);
  }
  for (size_t i = 0; i < lp_end_trans.size(); ++i) {
    lp_end_trans[i] = log(end_trans[i]);
  }
  for (size_t i = 0; i < lp_trans.size(); ++i) {
    lp_trans.resize(trans[i].size());
    for (size_t j = 0; j < trans[i].size(); ++j) {
      lp_trans[i][j] = log(trans[i][j]);
    }
  }

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double forward_score =
    forward_algorithm(meth, reset_points[i], reset_points[i + 1],
                      lp_start_trans, lp_end_trans, lp_trans);
    const double backward_score =
    backward_algorithm(meth, reset_points[i], reset_points[i + 1],
                       lp_start_trans, lp_end_trans, lp_trans);
    
    if (DEBUG && (fabs(forward_score - backward_score) /
                  max(forward_score, backward_score)) > 1e-10)
      cerr << "fabs(forward_score - backward_score)/"
      << "max(forward_score, backward_score) > 1e-10" << endl;

    total_score += forward_score;
  }
  
  classes.resize(data_size);

  llr_scores.resize(data_size);
  for (size_t i = 0; i < data_size; ++i) {
    
    int best_state = 0;
    double best_state_score = forward[0][i] + backward[0][i];
    double total_state_score = best_state_score;
    
    for (int s = 1; s < num_states; ++s) {
      int state_score = forward[s][i] + backward[s][i];
      if (state_score >= best_state_score) {
        best_state = s;
        best_state_score = state_score;
      }
      total_state_score = log_sum_log(total_state_score, state_score);
    }

    classes[i] = best_state;
    llr_scores[i] = exp(best_state_score - total_state_score);
  }

  return total_score;
}

