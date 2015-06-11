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

#include <cmath>

#include <iostream>
#include <vector>
#include <numeric>
#include <string>

#include <gsl/gsl_cdf.h>

#include "vdhmm-utils.hpp"
#include "Distro.hpp"

using std::string;
using std::cerr;
using std::endl;

// Maximization
STATE_LABELS 
max3(const double gain, const double same, const double loss)
{
    STATE_LABELS l = SAME;
    if (gain > same && gain > loss) l = GAIN;
    if (loss > gain && loss > same) l = LOSS;
    return l;
}

//////////////////////////////////////////////
//////       log_sum_log                //////
//////////////////////////////////////////////

double
log_sum_log_vec(const std::vector<double> &vals, const size_t limit) 
{
    const std::vector<double>::const_iterator x = 
        std::max_element(vals.begin(), vals.begin() + limit);
    const double max_val = *x;
    const size_t max_idx = x - vals.begin();
    double sum = 1.0;
    for (size_t i = 0; i < limit; ++i) 
    {
        if (i != max_idx) 
        {
            sum += exp(vals[i] - max_val);
        }
    }
    return max_val + log(sum);
}

double 
lp_segment_longer_than(
    const size_t len,
    const Distro& duration_distro)
{
    if (duration_distro.tostring().find("nbd") != string::npos)
    {
        const vector<double> params = duration_distro.get_params();
        const double n = 1 / params[1];
        const double p = n / (params[0] + n);
    
        return log(gsl_cdf_negative_binomial_Q(len, p, n));
    }
    else if (duration_distro.tostring().find("pois") != string::npos)
    {
        const vector<double> params = duration_distro.get_params();
        const double mu = params.front();
        
        return log(gsl_cdf_poisson_Q(len, mu));
    }
    else if (duration_distro.tostring().find("geo") != string::npos)
    {
        const vector<double> params = duration_distro.get_params();
        const double p = params.front();
        
        return log(gsl_cdf_geometric_Q(len, p));
    }
    else 
    {
        cerr << "BayBinomHMR::first_segment_len_log_prob:"
             <<" unsupported duration distribution" << endl;
        abort();
    }
}



