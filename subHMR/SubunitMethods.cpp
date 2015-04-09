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



///////////////////////////////////////////////////////////////////////////
// FUNCTIONS TO OPERATE BIT VECTORS
void
bv_print(const vector<int> &s) {
  for(size_t i = 0; i < s.size(); ++i) {
    std::cout << s[i] << " , ";
  }
  std::cout << std::endl;
}

vector<int>
bv_union(const vector<int> &s1, const vector<int> &s2) {
  vector<int> uni(s1.size());
  for(size_t i = 0; i < s1.size(); ++i) {
    uni[i] = s1[i] || s2[i];
  }
  return uni;
}

vector<int>
bv_diff(const vector<int> &s1, const vector<int> &s2) {
  // remove the intersection of s1 & s2 from s1
  vector<int> diff(s1.size());
  for(size_t i = 0; i < s1.size(); ++i) {
    diff[i] = s1[i] - (s1[i] && s2[i]);
  }
  return diff;
}

vector<int>
bv_exist(const vector<int> &s) {
  vector<int> exist;
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