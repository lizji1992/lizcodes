/*
 *    Part of SubunitFinder
 *
 *    Copyright (C) 2008 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Liz Ji
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

#include "UnitCluster.hpp"
#include <cmath>
#include <algorithm>
#include <vector>
using std::vector;


Panel::Panel(const vector<Subunit> &subunits, const pair<size_t, size_t> &b,
             const size_t &window_size) {
  size_t s = subunits[b.first].start;
  size_t e = subunits[b.second].end;
  xscale = static_cast<float> (e - s) / static_cast<float> (window_size);
  signal = vector<float> (window_size, 0);
  
  size_t idx_lhs = 0, idx_rhs = 0;
  float score = 0, max_score = 0;
  for (size_t i = b.first; i <= b.second; ++i) {
    idx_rhs = round(static_cast<float> (subunits[i].end - s + 1) /
                         xscale);
    idx_rhs = idx_rhs > 0 ? idx_rhs - 1: idx_rhs;
    score = subunits[i].score_sum();
    size_t lfill = idx_lhs < window_size? idx_lhs : window_size - 1;
    size_t rfill = idx_rhs < window_size? idx_rhs : window_size - 1;
    std::fill(signal.begin()+lfill, signal.begin()+rfill+1, score);
    idx_lhs = idx_rhs;
    if (score > max_score) max_score = score;
  }
  // Y normalization
  for (size_t i = 0; i < signal.size(); ++i) {
    signal[i] = signal[i] / max_score;
  }
  yscale = max_score;
}