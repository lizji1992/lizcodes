/* HMR simulator
 * Song Qiang <qiang.song@usc.edu> 2011
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

int
main(int argc, const char **argv) 
{
    try 
    {
//      static const bool REQUIRED = true;
        static const bool OPTIONAL = false;

        string outfile;
        string segment_file;
        string readdepth = "nbd 20 0.1";
        
        string hypo_meth = "beta 1 9";
        string hyper_meth = "beta 9 1";

        // string hmr_len_str = "nbd 20 0.2";
        // string nonhmr_len_str = "nbd 100 0.5";
        double ff_prob = 0.95;
        double bb_prob = 0.99;
        
        size_t hmrn = 10;
        string chrome_name = "chrZ";

        bool ALL_SEG = false;
        bool VERBOSE = false;
    
        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(argv[0], "A program to simulate BS-Seq data"
                               "hmrfile nonhmrfile");
        opt_parse.add_opt("out", 'o', "output file (BED format)", 
                          OPTIONAL, outfile);
        opt_parse.add_opt("segment", 's', "output file of segments", 
                          OPTIONAL, segment_file);
        opt_parse.add_opt("hmrn", 'n', "Number of HMRs", 
                          OPTIONAL, hmrn);
        // opt_parse.add_opt("hmr_len_distro", '\0', "Length distribution of HMRs", 
        //                   OPTIONAL, hmr_len_str);
        // opt_parse.add_opt("nonhmr_len_distro", '\0', "Length distribution of non-HMRs", 
        //                   OPTIONAL, nonhmr_len_str);
        opt_parse.add_opt("ff-trans", 'F',
                          "Foreground to foreground transition probability", 
                          OPTIONAL, ff_prob);
        opt_parse.add_opt("bb-trans", 'B',
                          "Backward to backward transition probability", 
                          OPTIONAL, bb_prob);
        opt_parse.add_opt("hypo-meth-prior", '\0',
                          "Methylation probability prior for hypomethylated region", 
                          OPTIONAL, hypo_meth);
        opt_parse.add_opt("hyper-meth-prior", '\0',
                          "Methylation probability prior for hypermethylated region", 
                          OPTIONAL, hyper_meth);
        opt_parse.add_opt("readdepth", '\0',
                          "Read depth distribution", 
                          OPTIONAL, readdepth);
        opt_parse.add_opt("all-segments", 'A',
                          "Output both hypo- and hyper-methylated segments ",
                          OPTIONAL, VERBOSE);
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

        
        vector<GenomicRegion> cpgs;
        Distro hypo_meth_distro(hypo_meth);
        Distro hyper_meth_distro(hyper_meth);
        Distro readdepth_distro(readdepth);

		Runif rng(time(NULL));
        readdepth_distro.seed(rng.runif(1, std::numeric_limits<int>::max()));
        hypo_meth_distro.seed(rng.runif(1, std::numeric_limits<int>::max()));
        hyper_meth_distro.seed(rng.runif(1, std::numeric_limits<int>::max()));

        const gsl_rng_type * T;
        gsl_rng * r;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        gsl_rng_set(r, rng.runif(1, std::numeric_limits<int>::max()));

        const size_t max_cpgs =
            static_cast<size_t>(hmrn / (1 - ff_prob) + hmrn / (1 - bb_prob)) * 4;


        // simulation; always assume the first and the last region are
        // non-HMRs
        std::ostream *out = outfile.empty() ?
            &std::cout : new std::ofstream(outfile.c_str());
        std::ostream *out_segments = segment_file.empty() ?
            0 : new std::ofstream(segment_file.c_str());
        
        // sample the first region
        bool is_fg = false;
        size_t switches = 0;
        size_t idx = 0;
        size_t seg_start = 0;
        
        while (switches < 2 * hmrn + 1 && cpgs.size() < max_cpgs)
        {
            if (is_fg)
            {
                size_t count = static_cast<size_t>(readdepth_distro.sample());
                double meth = hypo_meth_distro.sample();
                cpgs.push_back(GenomicRegion(chrome_name, idx, idx + 1,
                                             "CpG:" + smithlab::toa(count),
                                             noisy_meth(count, meth), '+'));
                idx += 2;
                
                const bool switch_to_bg = (gsl_ran_bernoulli (r, ff_prob) == 0);
                if (switch_to_bg)
                {
                    switches += 1;
                    is_fg = false;
                    if (out_segments)
                        *out_segments << cpgs[seg_start].get_chrom() << "\t"
                                      << cpgs[seg_start].get_start() << "\t"
                                      << cpgs.back().get_end() << "\t"
                                      << "HYPO" << endl;
                    seg_start = cpgs.size();
                }
            }
            else
            {
                size_t count = static_cast<size_t>(readdepth_distro.sample());
                double meth = hyper_meth_distro.sample();
                cpgs.push_back(GenomicRegion(chrome_name, idx, idx + 1,
                                             "CpG:" + smithlab::toa(count),
                                             noisy_meth(count, meth), '+'));
                idx += 2;
                
                const bool switch_to_fg = (gsl_ran_bernoulli (r, bb_prob) == 0);
                if (switch_to_fg)
                {
                    switches += 1;
                    is_fg = true;
                    if (out_segments  && ALL_SEG)
                        *out_segments << cpgs[seg_start].get_chrom() << "\t"
                                      << cpgs[seg_start].get_start() << "\t"
                                      << cpgs.back().get_end() << "\t"
                                      << "HYPER" << endl;
                    seg_start = cpgs.size();
                }
            }
        } // end while
        
        // output simulation result
        std::copy(cpgs.begin(), cpgs.end(),
                  std::ostream_iterator<GenomicRegion>(*out, "\n"));
        
        delete out;
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

