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

#ifndef TWO_STATE_CTHMM_HPP
#define TWO_STATE_CTHMM_HPP

#include "smithlab_utils.hpp"
#include "distribution.hpp"
#include <memory>

using std::vector;
using std::pair;

typedef vector< vector<double> > matrix;


class TwoVarHMM {
public:
  
  TwoVarHMM(const double tol, const double minprob, const size_t max_itr,
            const bool v, bool d = false) :
    tolerance(tol), MIN_PROB(minprob), max_iterations(max_itr),
    VERBOSE(v), DEBUG(d) {}
  
  
  void
  set_parameters(const BetaBin _fg_emission, const BetaBin _bg_emission,
                 const double _fg_rate, const double _bg_rate,
                 const double _p_sf, const double _p_sb,
                 const double _p_ft, const double _p_bt);
  
  
  double
  BaumWelchTraining(const vector<pair<double, double> > &meth,
                    const vector<size_t> &time);
  

  double
  PosteriorDecoding(const vector<pair<double, double> > &meth,
                    const vector<size_t> &time,
                    vector<int> &classes, vector<double> &llr_scores);
 
  
private:
  
  double
  single_iteration(const vector<pair<double, double> > &meth,
                   const vector<size_t> &time,
                   const vector<double> &meth_lp,
                   const vector<double> &unmeth_lp);

  double
  forward_algorithm(const vector<pair<double, double> > &meth,
                    const vector<size_t> &time,
                    const double lp_sf, const double lp_sb,
                    const double lp_ft, const double lp_bt,
                    matrix &ltp);
  double
  backward_algorithm(const vector<pair<double, double> > &meth,
                     const vector<size_t> &time,
                     const double lp_sf, const double lp_sb,
                     const double lp_ft, const double lp_bt,
                     const matrix &ltp);
  
  
  void
  update_trans_estimator(const vector<pair<double, double> > &meth,
                         const double total, matrix &te,
                         const matrix &ltp) const;
  
  void
  update_endprob(matrix &te);
  
  
  void
  estimate_emissions(const vector<double> &meth_lp,
                     const vector<double> &unmeth_lp);
  
  
  double
  log_sum_log(const double p, const double q) const;
  
  double
  log_sum_log_vec(const std::vector<double> &vals, size_t limit) const;
  

  //  HMM internal data
  
  matrix forward;
  matrix backward;
  matrix ltp; // store the calculated transition probabilities across CpGs
  
  ////////  model structure  ////////
  BetaBin fg_emission, bg_emission;

  double fg_rate, bg_rate;
  double a, b; // a=fg_rate/(fg_rate+bg_rate), b=fg_rate+bg_rate
  double p_sf, p_sb;
  double p_ft, p_bt;
  
  
  //  parameters
  double tolerance;
  double MIN_PROB;
  size_t max_iterations;
  bool VERBOSE;
  bool DEBUG;
  
  // internal data
  vector<double> lr;
};


#endif