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

#ifndef PROCESS_SUBUNIT_HPP
#define PROCESS_SUBUNIT_HPP

#include <string>
#include <vector>
#include "GenomicRegion.hpp"

using std::string;
using std::vector;
using std::pair;


struct Endpoint {
  Endpoint(string c, const size_t s, const bool isf) :
  chr(c), start(s), count(1), is_first(isf) {}
  bool operator<(const Endpoint &other) const {
    return (chr < other.chr ||
            (chr == other.chr &&
             (start < other.start ||
              (start == other.start &&
               is_first < other.is_first))));
  }
  bool operator==(const Endpoint &other) const {
    return (chr == other.chr &&
            start == other.start && is_first == other.is_first);
  }
  string chr;
  size_t start;
  size_t count;
  bool is_first;
  vector<int> source;
};


class Subunit {
public:
  Subunit(const string c, const size_t b, const size_t e, vector<int> n,
          size_t d, size_t ct) :
          chr(c), start(b), end(e), source(n), ncpg(d), count(ct), strand('+'),
          istart(0), iend(0), icover(0), freq(0), score(0) {}
  
  bool contains(const GenomicRegion &other) const;
  bool front(const GenomicRegion &other) const;
  bool behind(const GenomicRegion &other) const;
  bool operator>(const Subunit &other) const {
    return (score > other.score);
  }
  
  string chr;
  size_t start;
  size_t end;
  vector<int> source;
  vector<pair<size_t, size_t> > cpgbound;
  size_t ncpg;
  size_t count;
  char strand;
  size_t istart;
  size_t iend;
  size_t icover; // island cover
  float freq;
  float score;
};

// bit vector operation
vector<int>
bv_union(const vector<int> &s1, const vector<int> &s2);

void
bv_print(const vector<int> &s);

vector<int>
bv_diff(const vector<int> &s1, const vector<int> &s2);

vector<int>
bv_exist(const vector<int> &s);

// vector operations
vector<size_t>
operator+(const vector<size_t> &v1, const vector<size_t> &v2);

template< typename T >
void
set_pairvec_first(vector<pair<T, T> > &pv, const vector<T> &v);
template< typename T >
void
set_pairvec_second(vector<pair<T, T> > &pv, const vector<T> &v);

#endif
