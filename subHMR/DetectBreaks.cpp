/*    BinarizeCpG:
 *
 *    Copyright (C) University of Southern California
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

#include <cstring>
#include <string> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>
#include <cmath>
#include <algorithm>
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::pair;
using std::endl;
using std::cerr;
using std::sort;
using std::max;
using std::min;

///////////////////////////////////////////////////////////////////////////////
// DATA STRUCTURES

struct interval {
  interval(const string c, const size_t s, const size_t e, const string n) :
    chr(c), start(s), end(e), state(n), size(e-s), break_type("na") {}

  string chr;
  size_t start;
  size_t end;
  string state;
  size_t size;
  string break_type;
};

static interval
interval_from_str(string buffer) {
  std::istringstream os(buffer);
  string chr, state;
  size_t start, end;
  os >> chr >> start >> end >> state;
  return(interval(chr, start, end, state));
}

static void
load_units(const string &interval_file, vector<interval> &units) {
  std::fstream in(interval_file.c_str());
  string buffer;
  while (getline(in, buffer)) {
    units.push_back(interval_from_str(buffer));
  }
}

static void
separate_regions(const vector<interval> &units, vector<size_t> &reset_points) {
  for (size_t i = 0; i < units.size(); ++i) {
    const bool adjacent = (i > 0 && units[i].chr == units[i-1].chr) ? 
      units[i].start == units[i-1].end : false;
    if (!adjacent)
      reset_points.push_back(i);
  }
  reset_points.push_back(units.size());
}

static void
mark_break(vector<interval> &units, size_t mostleft, size_t mostright,
    size_t degree) {
  for (size_t i = mostleft+1; i < mostright; i++) {
    size_t search_start = max(mostleft, i-degree);
    size_t search_end = min(mostright, i+degree);
    bool over_domain = false;
    for (size_t j = search_start; j < i && !over_domain; j++) {
      for (size_t k = i+1; k <= search_end && !over_domain; k++) {
        if (units[j].state == units[k].state &&
          units[i].size < (min(units[j].size, units[k].size)) / 2) {
          units[i].break_type = "break";
          if (units[j].break_type != "break")
            units[j].break_type = "broken";
          if (units[k].break_type != "break")
            units[k].break_type = "broken";
        }
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////


static void
write_units(const vector<interval> &units, const string outfile) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  for (size_t i=0; i < units.size(); ++i) {
    out << units[i].chr << '\t' << units[i].start << '\t'
        << units[i].end << '\t' << units[i].state << '\t'
        << units[i].break_type << endl;
  }
}


////////////////////////////////////////////////////////////////////////////////

int
main(int argc, const char **argv) {

  try {

    /* FILES */
    string outfile;

    size_t degree = 4;
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "DetectBreaks",
                           "<interval_files>");
    opt_parse.add_opt("degree", 'd', "the degree of checking",
                      false , degree);
    opt_parse.add_opt("output", 'o', "output file", true, outfile);
    
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() < 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const vector<string> interval_file(leftover_args);
    /**********************************************************************/
    
    std::cout << "Load units ... " << endl;
    vector<interval> units;
    load_units(interval_file[0], units);
    std::cout << "Load units: over. \n" << endl;
    
    std::cout << "Set boundaries ... " << endl;
    vector<size_t> reset_points;
    separate_regions(units, reset_points);
    std::cout << "Set boundaries: over. \n" << endl;
    
    std::cout << "Mark breaks ... " << endl;
    for (size_t i = 0; i < reset_points.size()-1; ++i) {
      size_t mostleft = max(static_cast<size_t> (0), reset_points[i]);
      size_t mostright = min(units.size()-1, reset_points[i+1]-1);
      mark_break(units, mostleft, mostright, degree);
    }
    
    std::cout << "Mark breaks: over ... " << endl;
    
    std::cout << "Write units ..." << endl;
    write_units(units, outfile);
    std::cout << "Write units: over " << endl;
  }
  catch (SMITHLABException &e) {
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
