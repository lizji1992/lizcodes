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

#include "ProcSubunit.hpp"
#include <iostream>
#include <vector>
using std::vector;


bool
Subunit::contains(const GenomicRegion& other) const {
  return chr == other.get_chrom() && start <= other.get_start()
         && other.get_end() <= end;
}

bool
Subunit::front(const GenomicRegion& other) const {
  return chr == other.get_chrom() && end < other.get_start();
}

bool
Subunit::behind(const GenomicRegion& other) const {
  return chr == other.get_chrom() && other.get_end() < start;
}

bool
Subunit::empty() const {
  bool empty = true;
  for (size_t i = 0; i < source.size() && empty; ++i) {
    if (source[i]) empty = false;
  }
  return empty;
}

bool
Subunit::zero_density() const {
  bool zero_density = true;
  for (size_t i = 0; i < cpgbound.size() && zero_density; ++i) {
    if (cpgbound[i].first != 0) zero_density = false;
  }
  return zero_density;
}

float
Subunit::score_sum() const {
  return std::accumulate(score.begin(), score.end(), 0);
}

///////////////////////////////////////////////////////////////////////////
// FUNCTIONS TO OPERATE BIT VECTORS

vector<bool>
bv_union(const vector<bool> &s1, const vector<bool> &s2) {
  vector<bool> uni(s1.size());
  for(size_t i = 0; i < s1.size(); ++i) {
    uni[i] = s1[i] | s2[i];
  }
  return uni;
}

vector<bool>
bv_diff(const vector<bool> &s1, const vector<bool> &s2) {
  // remove the intersection of s1 & s2 from s1
  vector<bool> diff(s1.size());
  for(size_t i = 0; i < s1.size(); ++i) {
    diff[i] = s1[i] & ~(s1[i] & s2[i]);
  }
  return diff;
}

vector<bool>
bv_exist(const vector<bool> &s) {
  vector<bool> exist;
  for(size_t i = 0; i < s.size(); ++i) {
    if (s[i]) exist.push_back(i);
  }
  return exist;
}

vector<size_t>
operator+(const vector<size_t> &v1, const vector<size_t> &v2) {
  vector<size_t> v(v1.size());
  for(size_t i = 0; i < v1.size(); ++i) { v[i] = v1[i] + v2[i]; }
  return v;
}

template< typename T >
void
set_pairvec_first(vector<pair<T, T> > &pv, const vector<T> &v) {
  for(size_t i = 0; i < v.size(); ++i) { pv[i].first = v[i]; }
}

template< typename T >
void
set_pairvec_second(vector<pair<T, T> > &pv, const vector<T> &v) {
  for(size_t i = 0; i < v.size(); ++i) { pv[i].second = v[i]; }
}