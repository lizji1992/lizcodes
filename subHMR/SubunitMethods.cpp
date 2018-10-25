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

#include "SubunitMethods.hpp"

#include <string>
#include <vector>
#include "SubunitMethods.hpp"

using std::unordered_map;
using std::string;
using std::vector;
using std::pair;

using std::vector;
/*
StableMatching
StableMatching::StableMatching(const vector<double> sim) {
  size_t num = xid.size();
  unordered_map<int, string> x_id, y_id;
  vector< vector<double>> sxy, syx;
  vector<bool> x_engaged;
  vector<vector<bool> > x_proposed;
  vector<vector<bool> > y_suitors;
  
  vector<pair<int, int> > matching;
   static void
   sim2ranked_mates(const vector<vector<double> > &sim_matrix,
   vector<vector<size_t> > &mates) {
   mates.resize(sim_matrix.size());
   for (size_t i = 0; i < sim_matrix.size(); ++i) {
   vector<mate> rmates(sim_matrix[i].size() - 1);
   for (size_t j=0; j < sim_matrix[i].size(); ++j) {
   if (j != i) rmates.push_back(mate(j, sim_matrix[i][j]) );
   }
   sort(rmates.begin(), rmates.end());
   mates[i].resize(rmates.size());
   for (size_t j=0; j < rmates.size(); ++j) {mates[i][j] = rmates[j].idx;}
   }
   }
}
*/

/*
 static void
 calc_similarity(vector<Subunit> &smallunits,
 vector< vector<double> > &sim_matrix) {
 const size_t m = sim_matrix[0].size();
 for (size_t i = 0; i < smallunits.size(); ++i) {
 for (size_t j = 0; j < m; ++j) {
 for (size_t k = j; k < m; ++k) {
 sim_matrix[j][k] += smallunits[i].source[j] * smallunits[i].source[k];
 }
 }
 }
 // normalization
 for (size_t j = 0; j < m; ++j) {
 for (size_t k = j+1; k < m; ++k) {
 sim_matrix[k][j] = sim_matrix[j][k] / sim_matrix[k][k];
 sim_matrix[j][k] = sim_matrix[j][k] / sim_matrix[j][j];
 }
 sim_matrix[j][j] = 1;
 }
 }*/
