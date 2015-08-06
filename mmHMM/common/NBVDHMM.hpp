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

#ifndef TWO_STATE_HMM_HPP
#define TWO_STATE_HMM_HPP

#include "smithlab_utils.hpp"
#include "distribution.hpp"
#include <memory>

using std::vector;
using std::pair;

typedef vector< vector<double> > matrix;


class TwoVarHMM {
public:
  
  TwoVarHMM(const double mp, const double tol, const size_t max_itr,
            const bool v, bool d = false) :
    MIN_PROB(mp), tolerance(tol), max_iterations(max_itr),
    VERBOSE(v), DEBUG(d) {}
  
  
  void
  set_parameters(const BetaBin _fg_emission, const BetaBin _bg_emission,
                 const size_t _fg_mode, const size_t _bg_mode,
                 const double _fg_p, const double _bg_p,
                 const double _p_sf, const double _p_sb,
                 const double _p_ft, const double _p_bt);
  
  void
  update_emission_matrix();
  void
  update_transition_matrix();
  
  
  double
  BaumWelchTraining(const vector<pair<double, double> > &meth,
                    const vector<size_t> &reset_points);
  
  
  double
  PosteriorDecoding(const vector<pair<double, double> > &meth,
                    const vector<size_t> &reset_points,
                    vector<int> &classes, vector<double> &llr_scores);
  
  double
  ViterbiDecoding(const vector<pair<double, double> > &meth,
                  const vector<size_t> &reset_points,
                  vector<int> &classes) const;
    
  BetaBin
  get_fg_emission() const {return fg_emission;}

  BetaBin
  get_bg_emission() const {return bg_emission;}
  
  size_t
  get_fg_mode() const {return fg_mode;}
  
  double
  get_fg_p() const {return fg_p;}
  
  double
  get_bg_p() const {return bg_p;}
  
  double
  get_p_sf() const {return p_sf;}
  
  double
  get_p_sb() const {return p_sb;}
  
  double
  get_p_ft() const {return p_ft;}
  
  double
  get_p_bt() const {return p_bt;}
  
private:
  
  double
  single_iteration(const vector<pair<double, double> > &meth,
                   const vector<size_t> &reset_points,
                   const vector<double> &meth_lp,
                   const vector<double> &unmeth_lp);
  double
  compute_likelihood(const vector<pair<double, double> > &meth,
                     const vector<size_t> &reset_points,
                     const vector<double> &meth_lp,
                     const vector<double> &unmeth_lp);

  double
  forward_algorithm(const vector<pair<double, double> > &meth,
                    const size_t start, const size_t end,
                    const vector<double> &lp_s, const vector<double> &lp_t,
                    const vector< vector<double> > &lp_trans);
  double 
  backward_algorithm(const vector<pair<double, double> > &meth,
                     const size_t start, const size_t end,
                     const vector<double> &lp_s, const vector<double> &lp_t,
                     const vector< vector<double> > &lp_trans);
  
  void
  update_trans_estimator(const vector<pair<double, double> > &meth,
                         const size_t start, const size_t end,
                         const double total, vector<matrix> &et,
                         const vector<vector<double> > &lp_trans) const;
  
  void
  update_transitions(const vector<matrix> &te);
  
  void
  estimate_emissions(const vector<double> &meth_lp,
                     const vector<double> &unmeth_lp);
  
  
  double
  log_sum_log(const double p, const double q) const;
  
  double
  log_sum_log_vec(const std::vector<double> &vals, size_t limit) const;
  

  //  HMM internal data
  
  vector< vector<double> > forward;
  vector< vector<double> > backward;
  
  ////////  model structure  ////////
  size_t num_states;
  size_t fg_mode;
  size_t bg_mode;
  
  BetaBin fg_emission, bg_emission;
  vector<BetaBin> emission;
  
  double fg_p, bg_p; // failure probability (transition probability)
  double p_sf, p_sb;
  double p_ft, p_bt;
  
  matrix trans;
  vector<double> start_trans, end_trans;
  
  //  parameters
  double MIN_PROB;
  double tolerance;
  size_t max_iterations;
  bool VERBOSE;
  bool DEBUG;
  
  mutable size_t emission_correction_count;
};


vector<vector<double> >
initialize_trans(const size_t fg_mode, const double fg_p,
                 const size_t bg_mode, const double bg_p);
#endif
