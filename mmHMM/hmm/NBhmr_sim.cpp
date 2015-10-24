/* HMR simulator
 * Song Qiang <qiang.song@usc.edu> 2015-2016
 *
 * No coordinates of CpGs are taken into consideration. The simualted CpGs are 
 * put onto "1,2,3,4,5 ..." of "chromosome Z" in the output.
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
#include "Distro.hpp"
#include "RNG.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::max;
using std::min;
using std::pair;
using std::make_pair;

using std::ostream_iterator;
using std::ofstream;

struct cpg
{
  cpg() : chr("chrZ"), start(1), strand("+"), name("CpG"),
          meth(0), coverage(0) {}
  cpg(const string c, const size_t s, const double m, const size_t cov) :
      chr(c), start(s), strand("+"), name("CpG"), meth(m), coverage(cov) {}
  
  string tostring() {
    std::ostringstream s;
    s << chr << "\t" << start << "\t" << strand << "\t" << name << "\t"
      << meth << "\t" << coverage;
    return s.str();
  }
  
  string chr;
  size_t start;
  string strand;
  string name;
  double meth;
  size_t coverage;
};


struct interval
{
  interval() : chr("chrZ"), start(1), end(2), name("0"), num_cpgs(0),
               strand("+"), score(0) {}
  interval(const string c, const size_t e, const size_t l) :
           chr(c), start(e-l), end(e),
           name("0"), num_cpgs(0), strand("+"), score(0) {}
  
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


void
sample_interval_cpgs(vector<cpg> &cpgs,
                     size_t &interval_len, Distro &readdepth_distro,
                     Distro &meth_distro, string &chrome_name, size_t &loci) {
  
  for (size_t k = 0; k < interval_len; ++k) {
    size_t count = static_cast<size_t>(readdepth_distro.sample());
    count = count > 0 ? count : 1;
    double meth = meth_distro.sample();
    cpgs.push_back(cpg(chrome_name, loci++, noisy_meth(count, meth), count));
  }
}


int
main(int argc, const char **argv)
{
  try
  {
    static const bool REQUIRED = true;
    static const bool OPTIONAL = false;
    
    string outfile;
    string segment_file;
    
    string readdepth = "nbd 20 0.1";
    string hypo_meth = "beta 1 9";
    string hyper_meth = "beta 9 1";
    
    string fg_duration_str = "nbd 76 0.5";
    string bg_duration_str = "geo 0.005";
    //double ff_prob = 0.95;
    //double bb_prob = 0.99;
    
    size_t hmrn = 10;
    string chrome_name = "chrZ";
    
    bool VERBOSE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "A program to simulate BS-Seq data");
    opt_parse.add_opt("out", 'o', "output file (BED format)",
                      OPTIONAL, outfile);
    opt_parse.add_opt("segment", 's', "output file of segments",
                      REQUIRED, segment_file);
    opt_parse.add_opt("hmrn", 'n', "Number of HMRs",
                      OPTIONAL, hmrn);
    opt_parse.add_opt("fg_len_distro", '\0', "Length distribution of HMRs",
                      OPTIONAL, fg_duration_str);
    opt_parse.add_opt("bg_len_distro", '\0', "Length distribution of non-HMRs",
                       OPTIONAL, bg_duration_str);
    //opt_parse.add_opt("hypo-meth-prior", '\0',
    //                  "Methylation probability prior for hypomethylated region",
    //                  OPTIONAL, hypo_meth);
    //opt_parse.add_opt("hyper-meth-prior", '\0',
    //                  "Methylation probability prior for hypermethylated region",
    //                  OPTIONAL, hyper_meth);
    //opt_parse.add_opt("readdepth", '\0',
    //                  "Read depth distribution",
    //                  OPTIONAL, readdepth);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      OPTIONAL, VERBOSE);
    
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
    /****************** END COMMAND LINE OPTIONS *****************/
    
    
    vector<cpg> cpgs;
    vector<interval> hmrs;
    //vector<inverval> nonhmrs;
    
    Distro hypo_meth_distro(hypo_meth);
    Distro hyper_meth_distro(hyper_meth);
    Distro fg_duration(fg_duration_str);
    Distro bg_duration(bg_duration_str);
    Distro readdepth_distro(readdepth);
    
    cerr << fg_duration << "\t" << bg_duration << endl;
    
    Runif rng(time(NULL));
    readdepth_distro.seed(rng.runif(1, std::numeric_limits<int>::max()));
    hypo_meth_distro.seed(rng.runif(1, std::numeric_limits<int>::max()));
    hyper_meth_distro.seed(rng.runif(1, std::numeric_limits<int>::max()));
    fg_duration.seed(rng.runif(1, std::numeric_limits<int>::max()));
    bg_duration.seed(rng.runif(1, std::numeric_limits<int>::max()));
    
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, rng.runif(1, std::numeric_limits<int>::max()));
    
    //    const size_t max_cpgs =
    //         static_cast<size_t>(hmrn / (1 - ff_prob) + hmrn / (1 - bb_prob)) * 4;

    // simulation; always assume the first and the last region are
    // non-HMRs. Simulate foregrounds and backgrounds by turns.
    
    
    std::ostream *out = outfile.empty() ?
                        &std::cout : new std::ofstream(outfile.c_str());
    std::ostream *out_segments = segment_file.empty() ?
                                 0 : new std::ofstream(segment_file.c_str());
    
    
    
    // sample the first region(bg)
    size_t interval_len = bg_duration.sample();
    size_t loci = 1;

    sample_interval_cpgs(cpgs, interval_len, readdepth_distro,
                         hyper_meth_distro, chrome_name, loci);
    
    // sample HMRs
    vector<double> paras = fg_duration.get_params();
    size_t n = 1/paras[1];
    
    for (size_t i = 0; i < hmrn; ++i) {
      // sample foreground
      interval_len = fg_duration.sample() + n;
      sample_interval_cpgs(cpgs, interval_len, readdepth_distro,
                           hypo_meth_distro, chrome_name, loci);
      hmrs.push_back(interval(chrome_name, loci, interval_len));
      
      // sample background
      interval_len = bg_duration.sample() + 1;
      sample_interval_cpgs(cpgs, interval_len, readdepth_distro,
                           hyper_meth_distro, chrome_name, loci);
      
    }
    
    // write cpgs and intervals
    for (size_t i = 0; i < cpgs.size(); ++i) {
      *out << cpgs[i].tostring() << endl;
    }
    
    for (size_t i = 0; i < hmrs.size(); ++i) {
      *out_segments << hmrs[i].tostring() << endl;
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

