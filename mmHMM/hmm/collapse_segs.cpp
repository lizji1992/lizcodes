/*
 *  Continuous-time Variable-duration HMM
 *
 * Copyright (C) 2015-2016 University of Southern California
 *                         Andrew D Smith
 * Author: Andrew D. Smith, Song Qiang
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include <numeric>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <unistd.h>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"


using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::max;
using std::min;
using std::pair;
using std::setw;



static void
load_domains(const string &segments_file, vector<GenomicRegion> &domains) {
  
  std::ifstream in(segments_file.c_str());
  GenomicRegion read_domain;
  while (in >> read_domain) {
    string state = read_domain.get_name();
    state.erase(0, 1);
    read_domain.set_name(state);
    domains.push_back(read_domain);
  }
}


static void
build_hmr_domains(const vector<GenomicRegion> &domains,
                  vector<GenomicRegion> &hmrs, const string bg_class) {

  size_t prev_end = 0;
  size_t hmr_count = 0;
  size_t hmr_sum_cpgs = 0;
  bool in_HMR = false;
  
  for (size_t i = 0; i < domains.size(); ++i) {
    if (domains[i].get_name() == "0") {
      in_HMR = true;
      hmr_sum_cpgs = domains[i].get_score();
      hmr_count++;
      hmrs.push_back(GenomicRegion(domains[i]));
      hmrs.back().set_name(toa(hmr_count));
    } else if (domains[i].get_name() == bg_class) {
      if (in_HMR) {
        in_HMR = false;
        hmrs.back().set_end(prev_end);
        hmrs.back().set_score(hmr_sum_cpgs);
      }
    } else {
      hmr_sum_cpgs += domains[i].get_score();
    }
    prev_end = domains[i].get_end();
  }
}

int
main(int argc, const char **argv) {

  try {
    
    size_t desert_size = 1000;
    size_t fg_mode = 3;

    string outfile;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "Program for identifying "
			   "HMRs in methylation data", "<cpg-BED-file>");
    opt_parse.add_opt("out", 'o', "output hmr file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("fg_mode", 'F', "how many inner states in foreground",
                      false, fg_mode);
    opt_parse.add_opt("desert", 'd', "max dist btwn cpgs with reads in HMR",
                      false, desert_size);
    
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
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
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string segments_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    
    vector<GenomicRegion> domains;
    load_domains(segments_file, domains);
    
    vector<GenomicRegion> hmrs;
    build_hmr_domains(domains, hmrs, toa(fg_mode));
    
    
    // output HMRs
    std::ostream *out = outfile.empty() ?
    &std::cout : new std::ofstream(outfile.c_str());
    
    for (size_t i = 0; i < hmrs.size(); ++i) {
      *out << hmrs[i] << '\t' << 0 << endl;
    }
    
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
