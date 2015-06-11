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

#ifndef UNITCLUSTER_HPP
#define UNITCLUSTER_HPP

#include <string>
#include <vector>
#include "ProcSubunit.hpp"


using std::string;
using std::vector;
using std::pair;

class Panel{
public:
  Panel(const vector<float> s, const float x, const float y) :
        signal(s), xscale(x), yscale(y) {}
  Panel(const vector<Subunit> &subunits, const pair<size_t, size_t> &b,
        const size_t &window_size);

  vector<float> signal;
  float xscale;
  float yscale;
};

#endif
