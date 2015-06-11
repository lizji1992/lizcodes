/*    subUnitFiner: 
 *
 *    Copyright (C) 2014 Andrew D. Smith
 *    Authors: Andrew D. Smith and Liz Ji
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
#include <tr1/unordered_map>
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "ProcSubunit.hpp"
#include "UnitCluster.hpp"

using std::string;
using std::vector;
using std::pair;
using namespace std::tr1;
using std::endl;
using std::cerr;
using std::sort;



///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// DATA STRUCTURES

static void
load_units(const size_t &num_files, const string &interval_file,
           vector<Endpoint> &endpoints, unordered_map <size_t, string> &fidmap) {
    // assign sample id (id_idx)
    const int id_idx = fidmap.size();
    const string source = basename(interval_file);
    fidmap[id_idx] = source;
    // construct source vector;
    vector<bool> sources(num_files, false);
    sources[id_idx] = 1;
    //load endpoints of units and cpgs
    std::fstream in(interval_file.c_str());
    if (!in)
      throw BEDFileException("cannot open input file " + interval_file);
    string buffer;
    while (getline(in, buffer)) {
      GenomicRegion interval = GenomicRegion(buffer);
      Endpoint lhs = Endpoint(interval, true, sources);
      Endpoint rhs = Endpoint(interval, false, sources);
      endpoints.push_back(lhs);
      endpoints.push_back(rhs);
    }
}

static void
collapse_endpoints(vector<Endpoint> &ep) {
  size_t j = 0;
  for (size_t i = 1; i < ep.size(); ++i) {
    if (ep[i] == ep[j]) {
      ep[j].count++;
      ep[j].source = bv_union(ep[i].source, ep[j].source);
    }
    else ep[++j] = ep[i];
  }
  ep.erase(ep.begin() + j + 1, ep.end());
}

static void
part_intervals(vector<Endpoint> &endpoints, vector<Subunit> &subunits,
               vector<pair<size_t, size_t> > &island_bound) {

  collapse_endpoints(endpoints);

  vector<Subunit> newunits;
  vector<bool> inter(endpoints[0].source.size(), 0);

  size_t count = 0;
  size_t firstone = 0;
  size_t lastend = 0;
  size_t cover = 0;
  size_t j = 0;
  for (size_t i = 0; i < endpoints.size() - 1; ++i) {
    if (endpoints[i].is_first) {
      if (subunits.size() > 0 && (endpoints[i].start > lastend ||
          endpoints[i].chr != subunits.back().chr)) {
        for (j = firstone; j < subunits.size(); ++j) {
          subunits[j].istart = subunits[firstone].start;
          subunits[j].iend = subunits.back().end;
          subunits[j].icover = cover;
        }
        island_bound.push_back(pair<size_t, size_t> (firstone,
                                                     subunits.size() - 1));
        firstone = subunits.size();
      }
      count += endpoints[i].count;
      if (cover < count)  cover = count;
      inter = bv_union(inter, endpoints[i].source);
    } else {
      count -= endpoints[i].count;
      inter = bv_diff(inter, endpoints[i].source);
    }
    if (count >= 1) {
      Subunit addone = Subunit(endpoints[i].chr, endpoints[i].start,
                               endpoints[i+1].start, inter, count);
      addone.cpgbound = vector<pair<size_t, size_t> >(inter.size(),
                        pair<size_t, size_t> (0, 0));
      subunits.push_back(addone);
      lastend = endpoints[i+1].start;
    }
  } // assign final island
  for (size_t j = firstone; j < subunits.size(); ++j) {
    subunits[j].istart = subunits[firstone].start;
    subunits[j].iend = subunits.back().end;
    subunits[j].icover = cover;
  }
  island_bound.push_back(pair<size_t, size_t> (firstone, subunits.size() - 1));
}


////////////////////////////////////////////////////////////////////////////////
// SCORING FUNCTIONS
static void
score_a_subunit(Subunit &subunit) {
  for (size_t i = 0; i < subunit.source.size(); ++i) { // i sample
    if (subunit.source[i]) subunit.score[i] = 1;
  }
}


////////////////////////////////////////////////////////////////////////////////

static string
join_signal(const vector<float> &v, string sep) {
  string s;
  for(size_t i = 0; i < v.size(); ++i ) {
      s = s + toa(v[i]) + sep;
  }
  s.erase(s.end() - sep.size(), s.end());
  return s;
}


static void
islands_to_file(const string &outfile, const vector<Panel> &islands) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  for (size_t i = 0; i < islands.size(); ++i) {
     out << join_signal(islands[i].signal, ",") << '\t'
         << islands[i].xscale << '\t' << islands[i].yscale << endl;
  }
}



int
main(int argc, const char **argv) {

  try {

    /* FILES */
    string ppdir;
    string outfile;

    bool VERBOSE = false;
    size_t window_size = 100;
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "subUnitFinder", 
                           "<interval-files>");
    opt_parse.add_opt("postprob", 'p', "dir of post-probs files",
                    true , ppdir);
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("window", 'w', "window_size to scale the unit in",
                      false , window_size);
    // opt_parse.add_opt("verbose", 'v', "print more run info",
    //                false , VERBOSE);
    
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
    const vector<string> interval_files(leftover_args);
    /**********************************************************************/
    vector<Endpoint> endpoints;
    unordered_map <size_t, string> fidmap;
    size_t num_files = interval_files.size();

    for (size_t i = 0; i < num_files; ++i) {
      if (VERBOSE)
        cerr << interval_files[i] << endl;
      load_units(num_files, interval_files[i], endpoints, fidmap);
    }
    std::cout << "Got " << endpoints.size() << " end points." << endl;

    sort(endpoints.begin(), endpoints.end());
    std::cout << "Sort: over." << endl;

    vector<Subunit> subunits;
    vector<pair<size_t, size_t> > island_bound;
    part_intervals(endpoints, subunits, island_bound);
    std::cout << "Got " << subunits.size() << " subunits." << endl;
/*
    vector<vector<GenomicRegion> > cpgs;
    for (size_t i = 0; i < num_files; ++i) {
      load_cpgs(i, fidmap[i], ppdir, subunits, cpgs);
    }
    std::cout << "Load cpgs: over." << endl;
*/
    for (size_t i = 0; i < subunits.size(); ++i) {
      score_a_subunit(subunits[i]);
    }
    std::cout << "Calculate fragment scores: over." << endl;
  
    vector<Panel> islands;
    for (size_t i = 0; i < island_bound.size(); ++i) {
      Panel new_island = Panel(subunits, island_bound[i], window_size);
      islands.push_back(new_island);
    }
    std::cout << "Scale islands: over." << endl;
    std::cout << "Got " << islands.size() << " scaled islands." << endl;


    islands_to_file(outfile, islands);
    std::cout << "Write islands: over." << endl;
 
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
