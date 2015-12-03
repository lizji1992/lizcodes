/*
 * HMR simulator using continuous-time markov
 *
 * Copyright (C) 2015-2016 University of Southern California
 *                         Andrew D Smith
 *
 * Author: Andrew D. Smith, Xiaojing Ji
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

#include <ctime>
#include <numeric>
#include <cmath>
#include <fstream>
#include <iterator>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "RNG.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;

using std::ostream_iterator;
using std::ofstream;



struct cpg
{
  cpg() : chr("chrZ"), start(1), strand("+"), name("CpG"),
          meth(0), cov(0) {}
  cpg(const string c, const size_t s, const double m, const size_t cov) :
      chr(c), start(s), strand("+"), name("CpG"), meth(m), cov(cov) {}
  
  
  string tostring() {
    std::ostringstream s;
    s << chr << "\t" << start << "\t" << strand << "\t" << name << "\t"
      << meth << "\t" << cov;
    return s.str();
  }
  
  bool from_buffer(std::istream& s) {
    return (s >> chr >> start >> strand >> name >> meth >> cov);
  }
  
  string chr;
  size_t start;
  string strand;
  string name;
  double meth;
  size_t cov;
};


struct interval
{
  interval() : chr("chrZ"), start(1), end(2), name("HMR"), num_cpgs(0),
  strand("+"), score(0) {}
  interval(const string c, const size_t s, const size_t e, const size_t n) :
  chr(c), start(s), end(e), name("HMR"), num_cpgs(n), strand("+"), score(0) {}
  
  string tostring() {
    std::ostringstream s;
    s << chr << "\t" << start << "\t" << end << "\t" << name << "\t"
    << num_cpgs << "\t" << strand << "\t" << score;
    return s.str();
  }
  
  string chr;
  size_t start;
  size_t end;
  string name;
  size_t num_cpgs;
  string strand;
  double score;
};



static size_t
sample_read_depth(const gsl_rng *rng, const double n, const double p) {
  return static_cast<size_t>(gsl_ran_negative_binomial(rng, p, n));
}


static void
load_cpgs(const string &cpgs_file, vector<cpg> &cpgs,
          const gsl_rng *rng, const double n, const double p) {
  
  std::ifstream in(cpgs_file.c_str());
  cpg read_one_cpg;
  while (read_one_cpg.from_buffer(in)) {
    size_t sample_cov = sample_read_depth(rng, n, p);
    read_one_cpg.cov = sample_cov > 0 ? sample_cov : 1;
    cpgs.push_back(read_one_cpg);
  }
}


inline double
p00(const double a, const double b, const cpg i, const cpg j) {
  double p = 0;
  if (i.chr == i.chr) {
    p = 1-a + a*exp( -b * (j.start - i.start) );
  }
  return p;
}

inline double
p11(const double a, const double b, const cpg i, const cpg j) {
  double p = 0;
  if (i.chr == j.chr) {
    p = a + (1-a)*exp( -b * (j.start - i.start) );
  }
  return p;
}



double
noisy_meth(const size_t count, const double meth)
{
  static Runif rng(time(NULL));
  double n_meth = floor(meth * count);
  n_meth += static_cast<double>(meth * count - n_meth
                                > rng.runif(0.0, 1.0));
  double n_unmeth = count - n_meth;
  
  return n_meth / (n_meth + n_unmeth);
}


static void
sample_bg_cpgs(vector<cpg> &cpgs, size_t &idx,
               const gsl_rng *r, const gsl_rng *r_meth,
               const double a, const double b,
               const double alpha, const double beta) {
  size_t L = cpgs.size();
  if (idx < L) {
    // sample the first cpg
    double meth = gsl_ran_beta(r_meth, alpha, beta);
    cpgs[idx].meth = noisy_meth(cpgs[idx].cov, meth);
    if (idx < L-1) {
      bool stay = gsl_ran_bernoulli(r, p00(a, b, cpgs[idx], cpgs[idx+1])) == 1;
      idx++;
      while (idx < L-1 && stay) {
        meth = gsl_ran_beta(r_meth, alpha, beta);
        cpgs[idx].meth = noisy_meth(cpgs[idx].cov, meth);
        stay = gsl_ran_bernoulli(r, p00(a, b, cpgs[idx], cpgs[idx+1])) == 1;
        idx++;
      }
      if (stay) { // the last cpg?
        meth = gsl_ran_beta(r_meth, alpha, beta);
        cpgs[idx].meth = noisy_meth(cpgs[idx].cov, meth);
      }
    }
  }
}


static void
sample_fg_cpgs(vector<cpg> &cpgs, size_t &idx,
               const gsl_rng *r, const gsl_rng *r_meth,
               const double a, const double b,
               const double alpha, const double beta) {
  size_t L = cpgs.size();
  if (idx < L) {
    // sample the first cpg
    double meth = gsl_ran_beta(r_meth, alpha, beta);
    cpgs[idx].meth = noisy_meth(cpgs[idx].cov, meth);
    if (idx < L-1) {
      bool stay = gsl_ran_bernoulli(r, p11(a, b, cpgs[idx], cpgs[idx+1])) == 1;
      idx++;
      while (idx < L-1 && stay) {
        meth = gsl_ran_beta(r_meth, alpha, beta);
        cpgs[idx].meth = noisy_meth(cpgs[idx].cov, meth);
        stay = gsl_ran_bernoulli(r, p11(a, b, cpgs[idx], cpgs[idx+1])) == 1;
        idx++;
      }
      if (stay) { // the last cpg?
        meth = gsl_ran_beta(r_meth, alpha, beta);
        cpgs[idx].meth = noisy_meth(cpgs[idx].cov, meth);
      }
    }
  }
}



int
main(int argc, const char **argv)
{
  try
  {
    string outfile;
    string segment_file;
    
    double readdepth_distro_n = 20;
    double readdepth_distro_p = 0.7;
    
    double fg_rate = 0.0003;
    double bg_rate = 0.005;
    
    double bg_alpha = 2.4;
    double bg_beta = 0.6;
    double fg_alpha = 0.5;
    double fg_beta = 5.4;
    
    size_t num_hmr = std::numeric_limits<int>::max();
    
    bool VERBOSE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "A program to simulate BS-Seq data",
                           "<cpg-BED-file>");
    opt_parse.add_opt("out", 'o', "output file (BED format)",
                      true, outfile);
    opt_parse.add_opt("segment", 's', "output file of segments",
                      true, segment_file);
    opt_parse.add_opt("nhmr", 'n', "num of hmrs", false, num_hmr);
    opt_parse.add_opt("readdepth_n", '\0',
                      "Read depth negbin distribution n",
                      false, readdepth_distro_n);
    opt_parse.add_opt("readdepth_p", '\0',
                      "Read depth negbin distribution p",
                      false, readdepth_distro_p);
    opt_parse.add_opt("fg_rate", 'F', "foreground transition rate",
                      false, fg_rate);
    opt_parse.add_opt("bg_rate", 'B', "background transition rate",
                      false, bg_rate);
    opt_parse.add_opt("bg_alpha", '\0', "background alpha", false, bg_alpha);
    opt_parse.add_opt("bg_beta", '\0', "background beta", false, bg_beta);
    opt_parse.add_opt("fg_alpha", '\0', "foreground alpha", false, fg_alpha);
    opt_parse.add_opt("fg_beta", '\0', "foreground beta", false, fg_beta);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested())
    {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested())
    {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing())
    {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    
    const string cpgs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    if (VERBOSE)
      cerr << "READCOV N : " << readdepth_distro_n << endl
      << "READCOV P : " << readdepth_distro_p << endl
      << "F RATE : " << fg_rate << endl
      << "B RATE : " << bg_rate << endl
      << "F ALPHA : " << fg_alpha << endl
      << "F BETA : " << fg_beta << endl
      << "B ALPHA : " << bg_alpha << endl
      << "B BETA : " << bg_beta << endl;
    
    Runif rng(time(NULL));
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;

    gsl_rng *readdepth_rng = gsl_rng_alloc(T);
    gsl_rng_set(readdepth_rng, rng.runif(1, std::numeric_limits<int>::max()));
    
    
    vector<cpg> cpgs;
    load_cpgs(cpgs_file, cpgs,
              readdepth_rng, readdepth_distro_n, readdepth_distro_p);
    
    
    gsl_rng *r = gsl_rng_alloc (T);
    gsl_rng_set(r, rng.runif(1, std::numeric_limits<int>::max()));
    
    gsl_rng *r_fgmeth = gsl_rng_alloc(T);
    gsl_rng_set(r_fgmeth, rng.runif(1, std::numeric_limits<int>::max()));
    
    gsl_rng *r_bgmeth = gsl_rng_alloc(T);
    gsl_rng_set(r_bgmeth, rng.runif(1, std::numeric_limits<int>::max()));
    
    
    // simulation;
    vector<interval> hmrs;
    size_t idx = 0;
    
    double b = fg_rate + bg_rate;
    double a = fg_rate/b;
    // first region on each chromosome is always background
    sample_bg_cpgs(cpgs, idx, r, r_bgmeth, a, b, bg_alpha, bg_beta);
    
    size_t L = cpgs.size();
    for (size_t i = 0; i < num_hmr && idx < L; ++i) {
      const size_t begin_idx = idx;
      sample_fg_cpgs(cpgs, idx, r, r_fgmeth, a, b, fg_alpha, fg_beta);
      hmrs.push_back(interval(cpgs[begin_idx].chr, cpgs[begin_idx].start,
                              cpgs[idx-1].start, idx - begin_idx));
      
      sample_bg_cpgs(cpgs, idx, r, r_fgmeth, a, b, bg_alpha, bg_beta);
    }
    
    std::ostream *out = new std::ofstream(outfile.c_str());
    std::ostream *out_segments = new std::ofstream(segment_file.c_str());
    
    // write cpgs and intervals
    for (size_t i = 0; i < cpgs.size() && i < idx; ++i) {
      *out << cpgs[i].tostring() << endl;
    }
    
    for (size_t k = 0; k < hmrs.size(); ++k) {
      *out_segments << hmrs[k].tostring() << endl;
    }
    
    
  }
  catch (SMITHLABException &e)
  {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba)
  {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

