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

#ifndef SUBUNIT_METHODS_HPP
#define SUBUNIT_METHODS_HPP

#include <string>
#include <vector>
#include <tr1/unordered_map>
#include "GenomicRegion.hpp"
#include "ProcSubunit.hpp"


using std::string;
using std::vector;
using namespace std::tr1;
using std::pair;

/*
struct mate {
  mate(const size_t i, const double s) : idx(i), sim(s) {}
  bool operator<(const mate &other) const {
    return (sim < other.sim);
    
    size_t idx;
    double sim;
};
class StableMatching { // Gale–Shapley algorithm
public:
  StableMatching(const vector<vector <double> > sim);
  
  bool contains(const GenomicRegion &other) const;
  
  
  unordered_map<size_t, size_t> x_id, y_id;
  vector<vector<mate> > x_targets, y_suitors;
  vector<bool> x_engaged, x_unengaged;
  vector<pair<size_t, size_t> > matching;
};
*/
#endif
