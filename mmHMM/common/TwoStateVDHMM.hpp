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

#ifndef TWO_STATE_VDHMM_HPP
#define TWO_STATE_VDHMM_HPP

#include <memory>
#include <utility>
#include <string>
#include <vector>

#include "smithlab_utils.hpp"
#include "Distro.hpp"
#include "BetaBin.hpp"


class TwoStateVDHMM {
public:
  
    TwoStateVDHMM(const double mp, const double tol,
                  const size_t max_itr, const bool v,
                  const size_t _MAX_LEN) : 
        MAX_LEN(_MAX_LEN), MIN_PROB(mp), tolerance(tol),
        max_iterations(max_itr), VERBOSE(v) {}

    void
    init(const std::vector<std::pair<double, double> > &_observations,
         const std::vector<size_t> &_reset_points);
    
    void
    initialize_parameters(const betabin &_fg_emission,
                          const betabin &_bg_emission,
                          const Distro &_fg_duration,
                          const Distro &_bg_duration);
    
    double
    BaumWelchTraining();

    void
    get_posterior_scores(std::vector<double> & scores,
                         std::vector<bool> & classes);

private:

    //////////// methods ////////////
    double
    single_iteration();
    double 
    forward_algorithm(const size_t start, const size_t end);
    double 
    backward_algorithm(const size_t start, const size_t end);

    double
    fg_segment_log_likelihood(const size_t start, const size_t end);

    double
    bg_segment_log_likelihood(const size_t start, const size_t end);
    
    void
    estimate_state_posterior(const size_t start, const size_t end);

    void estimate_parameters();
    void update_observation_likelihood();

    double
    log_sum_log_vec(const std::vector<double> &vals, size_t limit) const;
    double
    log_sum_log(const double p, const double q) const;
  


    ////////   data   ////////
    std::vector<std::pair<double, double> > observations;
    std::vector<size_t> reset_points;
    std::vector<double> meth_lp, unmeth_lp;
    std::vector<double> fg_log_likelihood, bg_log_likelihood;
    

   //  HMM internal data 
    size_t MAX_LEN;
    betabin fg_emission,  bg_emission;
    Distro fg_duration,  bg_duration;
    
    double lp_sf, lp_sb, lp_ft, lp_bt;
    
    std::vector<std::pair<double, double> > forward;
    std::vector<std::pair<double, double> > backward;

    std::vector<double> fg_posteriors;
    std::vector<double> bg_posteriors;
    
    static const size_t FG_TO_BG_TRANSITION = 1;
    static const size_t BG_TO_FG_TRANSITION = 2;

    // parameters
    double MIN_PROB;
    double tolerance;
    size_t max_iterations;
    bool VERBOSE;
    bool DEBUG;
};

#endif
