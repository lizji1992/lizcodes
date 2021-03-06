/*    subUnitFiner: 
 *
 *    Copyright (C) 2015-2016 Andrew D. Smith
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
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "ProcSubunit.hpp"

using std::unordered_map;
using std::string;
using std::vector;
using std::pair;
using std::endl;
using std::cerr;
using std::sort;



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
part_intervals(const size_t &min_cover, vector<Endpoint> &endpoints,
               vector<Subunit> &smallunits) {

  collapse_endpoints(endpoints);

  vector<Subunit> newunits;
  vector<bool> inter(endpoints[0].source.size(), 0);

  size_t count = 0;
  size_t firstone = 0;
  size_t lastend = 0;
  size_t cover = 0;
  for (size_t i = 0; i < endpoints.size() - 1; ++i) {
    if (endpoints[i].is_first) {
      if (smallunits.size() > 0 && (endpoints[i].start > lastend ||
          endpoints[i].chr != smallunits.back().chr)) {
        for (size_t j = firstone; j < smallunits.size(); j++) {
          smallunits[j].istart = smallunits[firstone].start;
          smallunits[j].iend = smallunits.back().end;
          smallunits[j].icover = cover;
        }
        firstone = smallunits.size();
      }
      count += endpoints[i].count;
      if (cover < count)  cover = count;
      inter = bv_union(inter, endpoints[i].source);
    } else {
      count -= endpoints[i].count;
      inter = bv_diff(inter, endpoints[i].source);
    }
    if (count >= min_cover) {
      Subunit addone = Subunit(endpoints[i].chr, endpoints[i].start,
                               endpoints[i+1].start, inter, count);
      addone.cpgbound = vector<pair<size_t, size_t> >(inter.size(),
                        pair<size_t, size_t> (0, 0));
      smallunits.push_back(addone);
      lastend = endpoints[i+1].start;
    }
  } // assign final island
  for (size_t j = firstone; j < smallunits.size(); j++) {
    smallunits[j].istart = smallunits[firstone].start;
    smallunits[j].iend = smallunits.back().end;
    smallunits[j].icover = cover;
  }
}

static GenomicRegion
construct_GenomicRegion(string buffer) {
  std::istringstream os(buffer);
  string chr;
  size_t start, end;
  float score;
  os >> chr >> start >> end >> score;
  return GenomicRegion(chr, start, end, "X", score, '+');
}


static void
load_cpgs(const size_t &idx, const string &sample_name, const string &ppdir,
          vector<Subunit> &smallunits, vector<vector<GenomicRegion> > &cpgs) {
  const string pp_file = path_join(ppdir, sample_name + ".hmrpp");
  std::fstream in(pp_file.c_str());
  if (!in)
    throw BEDFileException("cannot open input file " + pp_file);
  vector<GenomicRegion> newcpgs;
  string buffer;
  bool hit_next = false;
  getline(in, buffer);
  GenomicRegion cpg = construct_GenomicRegion(buffer);
  for (size_t i = 0; i < smallunits.size() && in.good(); ++i) {
    while (smallunits[i].behind(cpg) && in.good()) {
      getline(in, buffer);
      cpg = construct_GenomicRegion(buffer);
    }
    if(smallunits[i].contains(cpg)) { // hit the first cpg
      if (!hit_next) // not overlapped with the last fragment
        smallunits[i].cpgbound[idx].first = newcpgs.size() + 1;
      while (smallunits[i].contains(cpg) && in.good()) {
        if (!hit_next) { // not overlapped with the last fragment
          smallunits[i].cpgbound[idx].second = newcpgs.size() + 1;
          // check if hit the next fragment
          if (i < smallunits.size()-1 && smallunits[i+1].contains(cpg)) {
            smallunits[i+1].cpgbound[idx].first = newcpgs.size() + 1;
            smallunits[i+1].cpgbound[idx].second = newcpgs.size() + 1;
            hit_next = true;
          }
          newcpgs.push_back(cpg);
        } else {
          hit_next = false;
        }
        getline(in, buffer);
        cpg = construct_GenomicRegion(buffer);
      } // hit the last cpg
    }
  }
  cpgs.push_back(newcpgs);
}

////////////////////////////////////////////////////////////////////////////////
// QUALITY CONTROL

static void
rm_junk(const float &min_del_freq, const size_t &min_size,
        vector<Subunit> &smallunits) {
  size_t j = 0;
  size_t num = smallunits[0].source.size();
  for (size_t i = 0; i < smallunits.size(); ++i) {
    float freq = static_cast<float> (smallunits[i].count) /
    static_cast<float> (num);
    const size_t unitsize = smallunits[i].end - smallunits[i].start;
    if (freq >= min_del_freq && unitsize >= min_size)
      smallunits[j++] = smallunits[i];
  }
  smallunits.erase(smallunits.begin() + j, smallunits.end());
}


static void
score_a_subunit(Subunit &subunit, const vector<vector<GenomicRegion> > &cpgs) {
  for (size_t i = 0; i < subunit.source.size(); ++i) { // i sample
    float p = 0, np = 0;
    if (subunit.cpgbound[i].first != 0) { // cover hmr cpg(s)
      for (size_t k = subunit.cpgbound[i].first - 1;
           k < subunit.cpgbound[i].second; ++k) {
        p += cpgs[i][k].get_score();
      }
      np = (subunit.cpgbound[i].second - subunit.cpgbound[i].first + 1) - p;
      subunit.score[i] = subunit.source[i] ? p : np;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

static string
join_source(const vector<bool> &s, string sep) {
  string js;
  for(size_t i = 0; i < s.size(); ++i ) {
    js = js + toa(s[i]) + sep;
  }
  js.erase(js.end()-1, js.end());
  return js;
}



static void
subunits_to_file(const string &outfile, const vector<Subunit> &subunits) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  for (size_t i = 0; i < subunits.size(); ++i) {
     out << subunits[i].chr << '\t'
         << subunits[i].start << '\t' << subunits[i].end << '\t'
         << subunits[i].istart << '\t' << subunits[i].iend << '\t'
         << join_source(subunits[i].source, ",") << '\t'
         << subunits[i].num_cpg() << '\t' << subunits[i].score_sum() << endl;
    }
}

static void
write_list(const string &listfile, unordered_map <size_t, string> &fidmap) {
  std::ofstream of;
  if (!listfile.empty()) of.open(listfile.c_str());
  std::ostream out(listfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  for (size_t i = 0; i < fidmap.size(); ++i) {
    out << i << '\t' << fidmap[i] << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////

int
main(int argc, const char **argv) {

  try {

    /* FILES */
    string ppdir;
    string outfile, listfile;

    bool VERBOSE = false;

    size_t min_cover = 1;
    float min_del_freq = 0.5;
    size_t min_size = 100;
    //float size_factor = 0.5;
    float degcutoff = 1;
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "subUnitFinder", 
                           "<interval-files>");
    opt_parse.add_opt("postprob", 'p', "dir of post-probs files",
                    true , ppdir);
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("list", 'l', "sample, id list",
                      false , listfile);
    opt_parse.add_opt("cover", 'c', "report intervals covered this many times",
                      false , min_cover);
    opt_parse.add_opt("freq", 'f', "the cutoff to drop small intervals",
                      false , min_del_freq);
    opt_parse.add_opt("min_size", 's', "min size", false, min_size);
    //opt_parse.add_opt("size", 'S', "the size factor", false, size_factor);
    opt_parse.add_opt("degcutoff", 'd', "merging degeneration", false, degcutoff);

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

    vector<Subunit> smallunits;
    part_intervals(min_cover, endpoints, smallunits);
    std::cout << "Got " << smallunits.size() << " smallunits." << endl;
    
    rm_junk(min_del_freq, min_size, smallunits);
    std::cout << "Got " << smallunits.size()
    << " smallunits left after throwing away too small ones." << endl;

    vector<vector<GenomicRegion> > cpgs;
    for (size_t i = 0; i < num_files; ++i) {
      load_cpgs(i, fidmap[i], ppdir, smallunits, cpgs);
    }
    std::cout << "Load cpgs: over." << endl;
    
    for (size_t i = 0; i < smallunits.size(); ++i) {
      score_a_subunit(smallunits[i], cpgs);
    }
    
    subunits_to_file(outfile, smallunits);
    write_list(listfile, fidmap);
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
