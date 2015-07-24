/*
 *  Continuous-time Variable-duration HMM
 *
 * Copyright (C) 2009-2012 University of Southern California
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
#include "NBVDHMM.hpp"
#include "distribution.hpp"


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
load_cpgs(const string &cpgs_file, vector<SimpleGenomicRegion> &cpgs,
          vector<pair<double, double> > &meth, vector<size_t> &reads)
{

  string chrom, prev_chrom;
  size_t pos, prev_pos = 0;
  string strand, seq;
  double level;
  size_t coverage;
  
  std::ifstream in(cpgs_file.c_str());
  while (in >> chrom >> pos >> strand >> seq >> level >> coverage) {
    // sanity check
    if (chrom.empty() || strand.empty() || seq.empty()
        || level < 0.0 || level > 1.0) {
      std::ostringstream oss;
      oss << chrom << "\t" << pos << "\t" << strand << "\t"
      << seq << "\t" << level << "\t" << coverage << "\n";
      throw SMITHLABException("Invalid input line:" + oss.str());
    }
    // order check
    if (prev_chrom > chrom || (prev_chrom == chrom && prev_pos > pos)) {
      throw SMITHLABException("CpGs not sorted in file \"" + cpgs_file + "\"");
    }
    prev_chrom = chrom;
    prev_pos = pos;
    
    // append site
    cpgs.push_back(SimpleGenomicRegion(chrom, pos, pos+1));
    reads.push_back(coverage);
    meth.push_back(std::make_pair(0.0, 0.0));
    meth.back().first = static_cast<size_t>(round(level * coverage));
    meth.back().second = static_cast<size_t>(coverage  - meth.back().first);
  }
}


template <class T, class U> static void
separate_regions(const bool VERBOSE, const size_t desert_size,
                 vector<SimpleGenomicRegion> &cpgs,
                 vector<T> &meth, vector<U> &reads,
                 vector<size_t> &reset_points) {
  if (VERBOSE)
    cerr << "[SEPARATING BY CPG DESERT]" << endl;
  // eliminate the zero-read cpgs
  size_t j = 0;
  for (size_t i = 0; i < cpgs.size(); ++i)
    if (reads[i] > 0) {
      cpgs[j] = cpgs[i];
      meth[j] = meth[i];
      reads[j] = reads[i];
      ++j;
    }
  cpgs.erase(cpgs.begin() + j, cpgs.end());
  meth.erase(meth.begin() + j, meth.end());
  reads.erase(reads.begin() + j, reads.end());
  
  // segregate cpgs
  size_t prev_cpg = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    const size_t dist = (i > 0 && cpgs[i].same_chrom(cpgs[i - 1])) ?
    cpgs[i].get_start() - prev_cpg : numeric_limits<size_t>::max();
    if (dist > desert_size)
      reset_points.push_back(i);
    prev_cpg = cpgs[i].get_start();
  }
  reset_points.push_back(cpgs.size());
  if (VERBOSE)
    cerr << "CPGS RETAINED: " << cpgs.size() << endl
    << "DESERTS REMOVED: " << reset_points.size() - 2 << endl << endl;
}


static void
get_domain_scores(const vector<int> &classes,
                  const vector<pair<double, double> > &meth,
                  const vector<size_t> &reset_points,
                  vector<double> &scores) {
  size_t n_cpgs = 0, reset_idx = 1;
  bool in_domain = false;
  double score = 0;
  int prev_class = classes[0];
  
  for (size_t i = 0; i < classes.size(); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
        in_domain = false;
        scores.push_back(score);
        score = 0;
      }
      ++reset_idx;
    }
    if (classes[i] == prev_class) {
      in_domain = true;
      score += 1.0 - (meth[i].first/(meth[i].first + meth[i].second));
      ++n_cpgs;
    }
    else if (in_domain) {
      in_domain = false;
      scores.push_back(score);
      score = 0;
    }
    prev_class = classes[i];
  }
}


static void
build_domains(const bool VERBOSE,
              const vector<SimpleGenomicRegion> &cpgs,
              const vector<double> &post_scores,
              const vector<size_t> &reset_points,
              const vector<int> &classes,
              vector<GenomicRegion> &domains) {
  
  size_t n_cpgs = 0, n_domains = 0, reset_idx = 1, prev_end = 0;
  bool new_domain = true;
  double score = 0;
  int prev_class = classes[0];
  
  for (size_t i = 0; i < classes.size(); ++i) {
    if (classes[i] != prev_class) {
      new_domain = true;
    }
    if (reset_points[reset_idx] == i) {
      new_domain = true;
      ++reset_idx;
    }
    if (new_domain){
      if (i != 0) {
        domains.back().set_end(prev_end);
        domains.back().set_score(n_cpgs);
      }
      domains.push_back(GenomicRegion(cpgs[i]));
      domains.back().set_name(toa(classes[i]));
      n_cpgs = 1;
      score = post_scores[i];
      new_domain = false;
    } else {
      ++n_cpgs;
      score += post_scores[i];
    }
    prev_end = cpgs[i].get_end();
    prev_class = classes[i];
  }
}

static void
build_hmr_domains(const bool VERBOSE, const vector<GenomicRegion> &domains,
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



static void
shuffle_cpgs(TwoVarHMM &hmm, vector<pair<double, double> > meth,
             vector<size_t> reset_points,
             vector<double> &domain_scores) {
  srand(time(0) + getpid());
  random_shuffle(meth.begin(), meth.end());
  vector<int> classes;
  vector<double> scores;
  hmm.PosteriorDecoding(meth, reset_points, classes, scores);
  get_domain_scores(classes, meth, reset_points, domain_scores);
  sort(domain_scores.begin(), domain_scores.end());
}

static void
assign_p_values(const vector<double> &random_scores,
                const vector<double> &observed_scores,
                vector<double> &p_values) {
  const double n_randoms =
      random_scores.size() == 0 ? 1 : random_scores.size();
  for (size_t i = 0; i < observed_scores.size(); ++i) {
    p_values.push_back((random_scores.end() -
                        upper_bound(random_scores.begin(),
                                    random_scores.end(),
                                    observed_scores[i]))/n_randoms);
  }
}


static void
write_params_file(const string &outfile, TwoVarHMM &hmm) {

  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  out.precision(30);
  out << "FG_ALPHA\t" << hmm.get_fg_emission().alpha << endl
      << "FG_BETA\t" << hmm.get_fg_emission().beta << endl
      << "BG_ALPHA\t" << hmm.get_bg_emission().alpha << endl
      << "BG_BETA\t" << hmm.get_bg_emission().beta << endl
      << "FG_MODE\t" << hmm.get_fg_mode() << endl
      << "S_F\t" << hmm.get_p_sf() << endl
      << "S_B\t" << hmm.get_p_sb() << endl
      << "F_IN\t" << 1 - hmm.get_fg_p() << endl
      << "F_OUT\t" << hmm.get_fg_p() << endl
      << "B_IN\t" << 1 - hmm.get_bg_p() << endl
      << "B_OUT\t" << hmm.get_bg_p() << endl
      << "F_E\t" << hmm.get_p_ft() << endl
      << "B_E\t" << hmm.get_p_bt() << endl;
}



int
main(int argc, const char **argv) {

  try {
    
    size_t desert_size = 1000;
    size_t max_iterations = 10;
    size_t fg_mode = 3;
    size_t bg_mode = 1;
    size_t mode_search_k = 3;
    
    // run mode flags
    bool VERBOSE = false;
    bool PARTIAL_METH = false;
    
    // corrections for small values (not parameters):
    double tolerance = 1e-10;
    double min_prob  = 1e-10;

    string params_in_file, params_out_file;
    string outfile, segments_file, scores_file;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "Program for identifying "
			   "HMRs in methylation data", "<cpg-BED-file>");
    opt_parse.add_opt("out", 'o', "output hmr file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("segment", 'S', "output hmr_seg file (default: stdout)",
                      false, segments_file);
    opt_parse.add_opt("scores", 's', "scores file (WIG format)",
                      false, scores_file);
    opt_parse.add_opt("desert", 'd', "max dist btwn cpgs with reads in HMR",
		      false, desert_size);
    opt_parse.add_opt("fg_mode", 'F', "how many inner states in foreground",
                      false, fg_mode);
    opt_parse.add_opt("bg_mode", 'B', "how many inner states in background",
                      false, bg_mode);
    opt_parse.add_opt("mode_search_k", 'R', "mode search range",
                      false, mode_search_k);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("partial", '\0', "identify PMRs instead of HMRs", 
		      false, PARTIAL_METH);
    opt_parse.add_opt("params-in", 'P', "HMM parameters file (no training)", 
		      false, params_in_file);
    opt_parse.add_opt("params-out", 'p', "write HMM parameters to this file", 
		      false, params_out_file);
    
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
    const string cpgs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    // separate the regions by chrom and by desert
    vector<SimpleGenomicRegion> cpgs;
    vector< pair<double, double> > meth;
    vector<size_t> reads;
    if (VERBOSE)
      cerr << "[READING CPGS AND METH PROPS]" << endl;
    load_cpgs(cpgs_file, cpgs, meth, reads);
    if (VERBOSE)
      cerr << "TOTAL CPGS: " << cpgs.size() << endl
      << "MEAN COVERAGE: "
      << accumulate(reads.begin(), reads.end(), 0.0)/reads.size()
      << endl << endl;
  
    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    vector<size_t> reset_points;
    
    separate_regions(VERBOSE, desert_size, cpgs, meth, reads, reset_points);
    
    // set-up distributions
    double fg_alpha = 0;
    double fg_beta = 0;
    double bg_alpha = 0;
    double bg_beta = 0;

    const double n_reads =
      accumulate(reads.begin(), reads.end(), 0.0)/reads.size();
    fg_alpha = 0.33*n_reads;
    fg_beta = 0.67*n_reads;
    bg_alpha = 0.67*n_reads;
    bg_beta = 0.33*n_reads;
    
    
    BetaBin fg_emission = BetaBin(fg_alpha, fg_beta);
    BetaBin bg_emission = BetaBin(bg_alpha, bg_beta);
    

    double fg_p = 0.25;
    double bg_p = 0.25;
    
    double p_sf = 0.5;
    double p_sb = 0.5;
    
    double p_ft = 1e-10;
    double p_bt = 1e-10;
   
  // HMM initialization
    
    double prev_score = -std::numeric_limits<double>::max();
    
    
    TwoVarHMM hmm(min_prob, tolerance, max_iterations, VERBOSE);
    hmm.set_parameters(fg_emission, bg_emission, fg_mode, bg_mode,
                       fg_p, bg_p, p_sf, p_sb, p_ft, p_bt);
    double score = hmm.BaumWelchTraining(meth, reset_points);
    
    // HMM training
    //if (max_iterations > 0)
    //  hmm.BaumWelchTraining(meth, reset_points);
/*
    
    if (!params_out_file.empty()) {
      // WRITE ALL THE HMM PARAMETERS:
      write_params_file(params_out_file, hmm);
    }
 */
 
    /***********************************
     * STEP 5: DECODE THE DOMAINS
     */
    
    vector<int> classes;
    vector<double> scores;
    hmm.PosteriorDecoding(meth, reset_points, classes, scores);
    
    
    vector<double> domain_scores;
    get_domain_scores(classes, meth, reset_points, domain_scores);
    
    
    vector<double> random_scores;
    shuffle_cpgs(hmm, meth, reset_points, random_scores);
    
    vector<double> p_values;
    assign_p_values(random_scores, domain_scores, p_values);
    
    vector<GenomicRegion> domains;
    build_domains(VERBOSE, cpgs, scores, reset_points, classes, domains);

    vector<GenomicRegion> hmrs;
    build_hmr_domains(VERBOSE, domains, hmrs, toa(fg_mode));
    
    // output HMR segments
    if (!segments_file.empty()) {
      std::ostream *out_seg = new std::ofstream(segments_file.c_str());
      for (size_t i = 0; i < domains.size(); ++i) {
        *out_seg << domains[i] << '\t' << p_values[i] << endl;
      }
    }
    
    // output HMRs
    std::ostream *out = outfile.empty() ?
                        &std::cout : new std::ofstream(outfile.c_str());
    
    for (size_t i = 0; i < hmrs.size(); ++i) {
      *out << hmrs[i] << '\t' << 0 << endl;
    }
    
    // output posterior probabilities
    if (!scores_file.empty()) {
      std::ostream *out_scores = new std::ofstream(scores_file.c_str());
      for (size_t i = 0; i < cpgs.size(); ++i) {
        *out_scores << cpgs[i] << '\t' << scores[i] << endl;
      }
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
