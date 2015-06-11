/*    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Jenny Qu
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/****************************************

 ***************************************/
#ifndef POST_APPROX_HPP
#define POST_APPROX_HPP

#include <string>
#include <vector>
#include <iostream>
#include <iomanip> // std::setw
#include <gsl/gsl_sf_gamma.h>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTree.hpp"
#include "HistMethEvo.hpp"
#include "BetaBin.hpp"


void
approx_posterior(const std::vector<std::vector<std::pair<double, double> > > &meth,
                 const std::vector<size_t> &reset_points, 
                 const HistoryME &history,  const size_t MAXITER,
                 std::vector<std::vector<double> > &node_posterior);

void 
candidate_hme(const size_t NHME,
              const std::vector<std::vector<double> > &node_posterior,
              std::vector<std::vector<std::vector<bool> > > &hmecand,
              std::vector<std::vector<double> > &hmecand_ps);

void 
trans_posterior(const HistoryME &history,
                const size_t NHME,
                const std::vector<std::vector<std::vector<bool> > > &hmecand,
                const std::vector<std::vector<double> > &hmecand_ps,
                const size_t pos,
                std::vector<std::vector<double> > &constrained_trans_post );

void
post_scores_to_coeff_app(const bool VERBOSE,
                         const std::vector<std::vector<std::pair<double, double> > > &values,
                         const std::vector<size_t> &reset_points,
                         const HistoryME &history,
                         const size_t NHME,
                         const std::vector<std::vector<std::vector<bool> > > &hmecand,
                         const std::vector<std::vector<double> > &hmecand_ps,
                         const std::vector<bool> &unknown,
                         std::vector<double> &coeff);


#endif
