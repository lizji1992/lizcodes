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

struct cpg {
  cpg(const string c, const size_t p) : chr(c), pos(p) {}
  
  size_t distance_to_next(cpg &other) const {
    size_t dist = 0;
    if (other.chr == chr) {
      dist = other.pos - pos;
    }
    else {
      dist = std::numeric_limits<size_t>::max();
    }
    return dist;
  }
  
  bool in_interval(const Subunit &interval) const{
    return(chr == interval.chr && pos >= interval.start && pos <= interval.end);
  }
  bool in_interval(const GenomicRegion &interval) const{
    return(chr == interval.get_chrom() &&
            pos >= interval.get_start() && pos <= interval.get_end());
  }
  
  string chr;
  size_t pos;
  vector<bool> member;
};

static cpg
cpg_from_str(string buffer) {
  std::istringstream os(buffer);
  string chr;
  size_t pos;
  os >> chr >> pos;
  return(cpg(chr, pos));
}

static void
load_cpgs(const string cpgfile, const string indexfile,
          const size_t desert_size, const size_t fill_num,
          const size_t num_files, vector<cpg> &cpgs, vector<string> &chrs) {
  std::fstream in(cpgfile.c_str());
  if (!in)
    throw BEDFileException("cannot open input file " + cpgfile);
  
  std::ofstream of;
  of.open(indexfile.c_str());
  
  string buffer;
  cpg lastcpg = cpg("chr1", 10468); // the first cpg
  string lastchr = "chrStart";
  while (getline(in, buffer)) {
    cpg newcpg = cpg_from_str(buffer);
    newcpg.member = vector<bool> (num_files, false);
    if (newcpg.chr != lastchr) {
      chrs.push_back(newcpg.chr);
      lastchr = newcpg.chr;
    }
    if (newcpg.chr == lastcpg.chr &&
        lastcpg.distance_to_next(newcpg) > desert_size) {
      size_t startpos = lastcpg.pos;
      size_t bin = lastcpg.distance_to_next(newcpg) / fill_num;
      for (size_t i=0; i < fill_num - 1; ++i) {
        cpg fillcpg = cpg(lastcpg.chr, startpos + bin);
        fillcpg.member = vector<bool> (num_files, false);
        cpgs.push_back(fillcpg);
        of << endl;
      }
    }
    cpgs.push_back(newcpg);
    lastcpg = newcpg;
    of << newcpg.chr << '\t' << newcpg.pos << endl;
  }
}


static void
load_units(const size_t &num_files, const string &interval_file,
           vector<cpg> &cpgs,
           const vector<string> chrs,
           vector <string> &sample_name) {
  // assign sample id (id_idx)
  const int id_idx = sample_name.size();
  const string source = basename(interval_file);
  sample_name.push_back(source);
  // construct source vector;
  vector<bool> sources(num_files, false);
  sources[id_idx] = 1;
  //load endpoints of units and cpgs
  std::fstream in(interval_file.c_str());
  if (!in)
    throw BEDFileException("cannot open input file " + interval_file);

  size_t i = 0;
  string buffer;
  pair<string, bool> lastchr(chrs[0], false); // lastchr in chrs?
  if (find(chrs.begin(), chrs.end(), lastchr.first) != chrs.end()) {
    lastchr.second = true;
  }
  while (getline(in, buffer)) {
    GenomicRegion interval = GenomicRegion(buffer);
    pair<string, bool> newchr(interval.get_chrom(), false);
    if (newchr.first == lastchr.first ) {
      newchr.second = lastchr.second;
    }
    else{
      if (find(chrs.begin(), chrs.end(), newchr.first) != chrs.end()) {
        newchr.second = true;
      }
    }
    if (newchr.second) {
      while (i < cpgs.size() && !cpgs[i].in_interval(interval)) {
        ++i;
      }
      while (i < cpgs.size() && cpgs[i].in_interval(interval)) {
        cpgs[i].member[id_idx] = true;
        ++i;
      }
    }
  }
}



////////////////////////////////////////////////////////////////////////////////

static string
join_source(const vector<bool> &s, const string sep) {
  string js = "";
  if (s.size() == 1) {
    js = toa(s[0]);
  }
  else if (s.size() > 1) {
    for(size_t i = 0; i < s.size() - 1; ++i ) {
      js = js + toa(s[i]) + sep;
    }
    js = js + toa(s[s.size()-1]);
  }
  return js;
}

static string
join_strings(const vector<string> &s, const string sep) {
  string js = "";
  if (s.size() == 1) {
    js = s[0];
  }
  else if (s.size() > 1) {
    for(size_t i = 0; i < s.size() - 1; ++i ) {
      js = js + s[i] + sep;
    }
    js = js + s[s.size()-1];
  }
  return js;
}

static void
write_signal(const vector<cpg> &cpgs, const vector <string> &sample_name,
             const string outdir) {
  string lastchr = cpgs[0].chr;
  string file = path_join(outdir, "subhmr_"+lastchr+"_binary.txt");
  std::ofstream of;
  of.open(file.c_str());
  of << "subhmr" << '\t' << lastchr << endl;
  of << join_strings(sample_name, "\t") << endl;

  for (size_t i=0; i < cpgs.size(); ++i) {
    if (cpgs[i].chr != lastchr) {
      of.close();
        // write a new chromosome
      lastchr = cpgs[i].chr;
      file = path_join(outdir, "subhmr_"+lastchr+"_binary.txt");
      of.open(file.c_str());
      of << "subhmr" << '\t' << lastchr << endl;
      of << join_strings(sample_name, "\t") << endl;
    }
    of << join_source(cpgs[i].member, "\t") << endl;
  }
}


////////////////////////////////////////////////////////////////////////////////

int
main(int argc, const char **argv) {

  try {

    /* FILES */
    string cpgfile, indexfile, outdir;


    float ppcutoff = 0.95;
    size_t desert_size = 1000;
    size_t fill_num = 20;
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "Binarize Cpgs",
                           "<interval_files>");
    opt_parse.add_opt("cpg", 'C', "CpG index file", true, cpgfile);
    opt_parse.add_opt("cutoff", 'c', "the cutoff to binarize value",
                      false , ppcutoff);
    opt_parse.add_opt("desert", 'd', "the size of desert",
                      false , desert_size);
    opt_parse.add_opt("fill", 'f', "the number of CpGs to fill in desert",
                      false , fill_num);
    opt_parse.add_opt("output", 'o', "output dir",
                      true, outdir);
    opt_parse.add_opt("map", 'm', "the indexes on genome of reported entries",
                      true, indexfile);
    
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
    
    std::cout << "Load CpGs ... " << endl;
    vector<string> chrs;
    vector<cpg> cpgs;
    size_t num_files = interval_files.size();
    load_cpgs(cpgfile, indexfile, desert_size, fill_num, num_files, cpgs, chrs);
    std::cout << "Load cpgs: over. \n" << endl;
    
    
    //unordered_map <string, vector<Endpoint> > endpoints;
    std::cout << "Load intervals ... " << endl;
    vector <string> sample_name;
    for (size_t i = 0; i < num_files; ++i) {
      load_units(num_files, interval_files[i], cpgs, chrs, sample_name);
    }
    std::cout << "Load intervals: over " << endl;
    write_signal(cpgs, sample_name, outdir);
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
