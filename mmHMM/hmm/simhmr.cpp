/* HMR simulator
 * Song Qiang <qiang.song@usc.edu> 2011
 */

#include <numeric>
#include <cmath>
#include <fstream>
#include <iterator>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

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
        static const bool REQUIRED = true;
        static const bool OPTIONAL = false;

        string outfile;
        string segment_file;
        string readdepth = "nbd 20 0.1";
        
        string hypo_meth = "beta 1 9";
        string hyper_meth = "beta 9 1";

        string hmr_len_str = "nbd 20 0.2";
        string nonhmr_len_str = "nbd 100 0.5";
        string nonhmr_len_file = "";
        
        size_t hmrn = 10;
        string chrome_name = "chrZ";

        bool VERBOSE = false;
        bool ALL_SEG = false;
    
        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(argv[0], "A program to simulate BS-Seq data"
                               "");
        opt_parse.add_opt("out", 'o', "output file (BED format)", 
                          OPTIONAL, outfile);
        opt_parse.add_opt("segment", 's', "output file of segments", 
                          OPTIONAL, segment_file);
        opt_parse.add_opt("hmrn", 'n', "Number of HMRs", 
                          OPTIONAL, hmrn);
        opt_parse.add_opt("hmr_len_distro", '\0', "Length distribution of HMRs", 
                          OPTIONAL, hmr_len_str);
        opt_parse.add_opt("nonhmr_len_distro", '\0', "Length distribution of non-HMRs", 
                          OPTIONAL, nonhmr_len_str);
        opt_parse.add_opt("nonhmr_len_file", '\0', "Filename of lengths  of non-HMRs", 
                          OPTIONAL, nonhmr_len_file);
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
        Distro hmr_len_distro(hmr_len_str);
        Distro nonhmr_len_distro(nonhmr_len_str);
        Distro hypo_meth_distro(hypo_meth);
        Distro hyper_meth_distro(hyper_meth);
        Distro readdepth_distro(readdepth);

        if (!nonhmr_len_file.empty())
        {
            nonhmr_len_distro = Distro("discemp");
            vector<double> lengths;
            double len;
            std::ifstream lenfile(nonhmr_len_file.c_str());
            while (lenfile >> len) lengths.push_back(len);
            nonhmr_len_distro.estimate_params_ml(lengths);
            lenfile.close();
            lengths.clear();
        }

		Runif rng(time(NULL) + getpid());
        hmr_len_distro.seed(rng.runif(1, std::numeric_limits<int>::max()));
        nonhmr_len_distro.seed(rng.runif(1, std::numeric_limits<int>::max()));
        readdepth_distro.seed(rng.runif(1, std::numeric_limits<int>::max()));
        hypo_meth_distro.seed(rng.runif(1, std::numeric_limits<int>::max()));
        hyper_meth_distro.seed(rng.runif(1, std::numeric_limits<int>::max()));

        
        // simulation; always assume the first and the last region are
        // non-HMRs
        std::ostream *out = outfile.empty() ?
            &std::cout : new std::ofstream(outfile.c_str());
        std::ostream *out_segments = segment_file.empty() ?
            0 : new std::ofstream(segment_file.c_str());
        size_t nloci;
        
        // sample the first region
        size_t idx = 0;
        nloci = static_cast<size_t>(nonhmr_len_distro.sample());
        for (size_t i = 0; i < nloci; ++i)
        {
            size_t count = static_cast<size_t>(readdepth_distro.sample());
            double meth = hyper_meth_distro.sample();
            cpgs.push_back(GenomicRegion(chrome_name, idx, idx + 1,
                                         "CpG:" + smithlab::toa(count),
                                         noisy_meth(count, meth), '+'));
            idx += 2;
        }
        if (out_segments && ALL_SEG)
            *out_segments << cpgs[idx / 2 - nloci].get_chrom() << "\t"
                          << cpgs[idx / 2 - nloci].get_start() << "\t"
                          << cpgs.back().get_end() << "\t"
                          << "HYPER" << endl;


        for (size_t i = 0; i < hmrn; ++i)
        {
            // sample HMR
            while((nloci = static_cast<size_t>(hmr_len_distro.sample())) == 0);
            for (size_t i = 0; i < nloci; ++i)
            {
                size_t count = static_cast<size_t>(readdepth_distro.sample());
                double meth = hypo_meth_distro.sample();
                cpgs.push_back(GenomicRegion(chrome_name, idx, idx + 1,
                                             "CpG:" + smithlab::toa(count),
                                             noisy_meth(count, meth), '+'));
                idx += 2;
            }
            if (out_segments)
                *out_segments << cpgs[idx / 2 - nloci].get_chrom() << "\t"
                              << cpgs[idx / 2 - nloci].get_start() << "\t"
                              << cpgs.back().get_end() << "\t"
                              << "HYPO" << endl;
            
            // sample non-HMR
            while((nloci = static_cast<size_t>(nonhmr_len_distro.sample())) == 0);
            for (size_t i = 0; i < nloci; ++i)
            {
                size_t count = static_cast<size_t>(readdepth_distro.sample());
                double meth = hyper_meth_distro.sample();
                cpgs.push_back(GenomicRegion(chrome_name, idx, idx + 1,
                                             "CpG:" + smithlab::toa(count),
                                             noisy_meth(count, meth), '+'));
                idx += 2;
            }
            if (out_segments && ALL_SEG)
                *out_segments << cpgs[idx / 2 - nloci].get_chrom() << "\t"
                              << cpgs[idx / 2 - nloci].get_start() << "\t"
                              << cpgs.back().get_end() << "\t"
                              << "HYPER" << endl;
        }


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

