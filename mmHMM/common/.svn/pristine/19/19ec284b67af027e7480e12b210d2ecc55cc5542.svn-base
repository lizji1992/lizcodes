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

#ifndef VDHMM_UTILS_HPP
#define VDHMM_UTILS_HPP

#include <cmath>

#include <vector>
#include <numeric>

#include "Distro.hpp"

enum STATE_LABELS {GAIN, SAME, LOSS};

struct Triplet
{
    double gain, same, loss;
};

// Maximization
STATE_LABELS 
max3(const double gain, const double same, const double loss);

inline double
log_sum_log(const double p, const double q)
{
    if (p == 0) {return q;}
    else if (q == 0) {return p;}
    const double larger = (p > q) ? p : q;
    const double smaller = (p > q) ? q : p;
    return larger + log(1.0 + exp(smaller - larger));
}

inline double
log_sum_log(const double p, const double q, const double r) 
{
    return log_sum_log(log_sum_log(p, q), r);
}

double
log_sum_log_vec(const std::vector<double> &vals, const size_t limit);

double 
lp_segment_longer_than(
    const size_t len,
    const Distro& duration_distro);
#endif
