/*
 *    Copyright (C) 2015 University of Southern California
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
#include <limits>
#include <cmath>
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "ProcSubunit.hpp"

using std::string;
using std::vector;
using std::pair;
using std::endl;
using std::cerr;
using std::sort;
using std::unordered_map;


///////////////////////////////////////////////////////////////////////////
// DATA STRUCTURES

struct cpg {
  cpg(const string c, const size_t s) : chr(c), start(s) {}
  
  bool operator<(const cpg &other) const {
    return (chr < other.chr ||
            (chr == other.chr &&
             (start < other.start)));
  }
  bool operator==(const cpg &other) const {
    return (chr == other.chr && start == other.start);
  }
  string chr;
  size_t start;
};

struct interval {
  interval(const string c, const size_t s, const size_t e, const string t) :
    chr(c), start(s), end(e), state(t), num_cpg(0) {}
  
  bool contains_cpg(const cpg &target) const{
    return(chr == target.chr && start <= target.start && end >= target.start);
  }
  bool operator<(const interval &other) const {
    return (chr < other.chr ||
            (chr == other.chr && start < other.start) ||
             (chr == other.chr && start == other.start && end < other.end));
  }
  bool operator==(const interval &other) const {
    return (chr == other.chr && start == other.start && end == other.end);
  }
  
  string chr;
  size_t start;
  size_t end;
  string state;
  size_t num_cpg;
};

static cpg
cpg_from_buffer(string buffer) {
  std::istringstream os(buffer);
  string chr;
  size_t start;
  os >> chr >> start;
  return cpg(chr, start);
}

static interval
interval_from_buffer(string buffer) {
  std::istringstream os(buffer);
  string chr, state;
  size_t start, end;
  os >> chr >> start >> end >> state;
  return interval(chr, start, end, state);
}


static void
load_cpgs(const string &cpgfile, vector<cpg> &cpgs) {
    std::fstream in(cpgfile.c_str());
    if (!in)
      throw BEDFileException("cannot open input file " + cpgfile);
    string buffer;
    while (getline(in, buffer)) {
      cpgs.push_back(cpg_from_buffer(buffer));
    }
}


static void
load_intervals(const string &interval_file, vector<interval> &intervals) {
  std::fstream in(interval_file.c_str());
  if (!in)
    throw BEDFileException("cannot open input file " + interval_file);
  string buffer;
  while (getline(in, buffer)) {
    intervals.push_back(interval_from_buffer(buffer));
  }
}


static void
count_cpgs(vector<interval> &intervals, const vector<cpg> &cpgs) {
  size_t cpg_pointer = 0;
  for (size_t i=0; i < intervals.size(); ++i) {
    while (cpg_pointer < cpgs.size() &&
           !intervals[i].contains_cpg(cpgs[cpg_pointer])) {
      cpg_pointer++;
    }
    while (cpg_pointer < cpgs.size() &&
           intervals[i].contains_cpg(cpgs[cpg_pointer])) {
      intervals[i].num_cpg++;
      cpg_pointer++;
    }
  }
}
  

static void
intervals_to_file(const string &outfile, const vector<interval> &intervals) {
  std::ofstream of;
  of.open(outfile.c_str());
  
  for (size_t i = 0; i < intervals.size(); ++i) {
     of << intervals[i].chr << '\t'
        << intervals[i].start << '\t' << intervals[i].end << '\t'
        << intervals[i].state << '\t' << intervals[i].num_cpg << endl;
  }
}


////////////////////////////////////////////////////////////////////////////////

int
main(int argc, const char **argv) {

  try {

    /* FILES */
    string cpgfile;
    string outfile;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "Assign CpGs to intervals",
                           "<interval-file>");
    opt_parse.add_opt("cpg", 'c', "Cpg index file", true , cpgfile);
    opt_parse.add_opt("output", 'o', "Out put file", true , outfile);

    
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
    const string interval_file = leftover_args[0];
    /**********************************************************************/

    vector<cpg> cpgs;
    load_cpgs(cpgfile, cpgs);
    std::cout << "Got " << cpgs.size() << " CpGs." << endl;
    //sort(cpgs.begin(), cpgs.end());
    //std::cout << "Sort: over." << endl;

    vector<interval> intervals;
    load_intervals(interval_file, intervals);
    std::cout << "Got " << intervals.size() << " intervals." << endl;
    //sort(intervals.begin(), intervals.end());
    //std::cout << "Sort: over." << endl;
    
    count_cpgs(intervals, cpgs);
    std::cout << "Count CpGs: over." << endl;
    intervals_to_file(outfile, intervals);
    std::cout << "Write: over." << endl;

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
