/* HMR simulator
 * Song Qiang <qiang.song@usc.edu> 2011
 */

#include <numeric>
#include <cmath>
#include <fstream>
#include <iterator>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "Distro.hpp"
#include "RNG.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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

static size_t
read_count(const double c, const double ratio)
{
    static bool first  = true;
    if (first)
    {
        gsl_rng_env_setup();
        first = false;
    }
    static gsl_rng * r = gsl_rng_alloc(gsl_rng_default);
    
    const size_t i = static_cast<size_t>(floor(c * ratio));
    const double f = c - i;
    size_t fi = gsl_ran_bernoulli(r, f);
    return i + fi;
}

static double
noisy_meth(const size_t count, double meth)
{
    static Distro noise("beta 10 10");
    static bool noise_generator_initialized = false;
    if (!noise_generator_initialized)
    {
        noise.seed(time(NULL) + getpid());
        noise_generator_initialized = true;
    }
    static Runif rng(time(NULL) + getpid());
    meth += noise.sample() - 0.5;
    if (meth < 0) meth = 0.01;
    if (meth > 1) meth = 0.99;
    
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
        string outfile;

        string cpg_file;
        string template_hmr_file;
        string hmr_file;
        string hmr_len_str;

        string hypo_meth = "beta 0.2 1.8";
        string hyper_meth = "beta 1.8 0.2";
        double coverage = -1;
        size_t hmrn = 0;

        bool VERBOSE = false;
    
        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(argv[0],
                               "A program to simulate BS-Seq data",
                               "");
        opt_parse.add_opt("out", 'o', "output file (BED format)", 
                          OptionParser::OPTIONAL, outfile);
        opt_parse.add_opt("template-hmr", '\0', "File contains HMRs", 
                          OptionParser::OPTIONAL, template_hmr_file);
        opt_parse.add_opt("hmr", '\0', "HMR output filename", 
                          OptionParser::OPTIONAL, hmr_file);
        opt_parse.add_opt("cpg", '\0',
                          "File contains coverage and methylation for each CpG", 
                          OptionParser::OPTIONAL, cpg_file);
        opt_parse.add_opt("hmrn", '\0', "Number of HMRs", 
                          OptionParser::OPTIONAL, hmrn);
        opt_parse.add_opt("hmr_len_str", '\0', "HMR length distribution", 
                          OptionParser::OPTIONAL, hmr_len_str);
        opt_parse.add_opt("coverage", 'c', "Desired mean coverage", 
                          OptionParser::OPTIONAL, coverage);
        opt_parse.add_opt("hypo-meth-prior", '\0',
                          "Methylation probability prior for hypomethylated region", 
                          OptionParser::OPTIONAL, hypo_meth);
        opt_parse.add_opt("hyper-meth-prior", '\0',
                          "Methylation probability prior for hypermethylated region", 
                          OptionParser::OPTIONAL, hyper_meth);
        opt_parse.add_opt("verbose", 'v', "print more run info",
                          OptionParser::OPTIONAL, VERBOSE);
        
        vector<string> leftover_args;
        opt_parse.parse(argc, argv, leftover_args);
        if (opt_parse.help_requested()) 
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
    
        // Read in genomic locis inside and outside of HMRs
        vector<GenomicRegion> template_cpgs;
        if (VERBOSE)
            cerr << "Reading template CpG profile from " << cpg_file << " ... ";
        ReadBEDFile(cpg_file, template_cpgs);
        check_sorted(template_cpgs);
        if (VERBOSE)
            cerr << "Done." << endl;
        
        vector<double> coverages(template_cpgs.size());
        for (size_t i = 0; i < template_cpgs.size(); ++i)
        {
            const string name = template_cpgs[i].get_name();
            coverages[i] = atof(name.substr(name.find(":") + 1).c_str());
        }


        // read in location of HMRs
        vector<SimpleGenomicRegion> hmrs;
        if (VERBOSE)
            cerr << "Building HMRs: " ;
        if (!template_hmr_file.empty())
        {
            if (VERBOSE)
                cerr << "reading HMRs from " << template_hmr_file << " ... ";
            ReadBEDFile(template_hmr_file, hmrs);
        }
        else
        {
            if (VERBOSE)
                cerr << "simulating HMRs based on length distribution "
                     << hmr_len_str << " ... ";
            Runif rng(time(NULL) + getpid());
            Distro hmr_len_distro(hmr_len_str);
            hmr_len_distro.seed(rng.runif(1, std::numeric_limits<int>::max()));

            // assume start of HMRs are uniformly distributed in the genome
            vector<size_t> idxs(template_cpgs.size());
            for (size_t i = 0; i < idxs.size(); ++i)
                idxs[i] = i;
            srand(std::time(NULL) + getpid());
            std::random_shuffle(idxs.begin(), idxs.end());
            idxs.erase(idxs.begin() + hmrn, idxs.end());
            std::sort(idxs.begin(), idxs.end());
            for (size_t i = 0; i < idxs.size(); ++i)
            {
                const size_t len = static_cast<size_t>(hmr_len_distro.sample());
                const size_t start_idx = idxs[i];
                size_t end_idx = start_idx + len - 1;
                while (end_idx > template_cpgs.size() - 1
                       || template_cpgs[end_idx].get_chrom()
                       != template_cpgs[start_idx].get_chrom())
                    --end_idx;
                
                hmrs.push_back(
                    SimpleGenomicRegion(template_cpgs[start_idx].get_chrom(),
                                        template_cpgs[start_idx].get_start(),
                                        template_cpgs[end_idx].get_end()));
            }
            
            // remove overlap
            size_t j = 0;
            for (size_t i = 0; i < hmrs.size(); ++i)
                if (hmrs[j].overlaps(hmrs[i]))
                    hmrs[j].set_end(hmrs[i].get_end());
                else
                    ++j;
            hmrs.erase(hmrs.begin() + j, hmrs.end());
        }
        if (hmrn != 0 && hmrn < hmrs.size())
        {
            srand(std::time(0) + getpid());
            std::random_shuffle(hmrs.begin(), hmrs.end());
            hmrs.erase(hmrs.begin() + hmrn, hmrs.end());
            std::sort(hmrs.begin(), hmrs.end());
        }
        check_sorted(hmrs);
        if (!hmr_file.empty()) WriteBEDFile(hmr_file, hmrs);
        if (VERBOSE) cerr << "Done." << endl;
        
        // begin simulation
        if (VERBOSE)
            cerr << "Simulating CpG methylation profile ...";
		Runif rng(time(NULL) + getpid());
        Distro hypo_meth_distro(hypo_meth);
        Distro hyper_meth_distro(hyper_meth);
        hypo_meth_distro.seed(rng.runif(1, std::numeric_limits<int>::max()));
        hyper_meth_distro.seed(rng.runif(1, std::numeric_limits<int>::max()));
        
        vector<GenomicRegion> simulated_cpgs(template_cpgs);
        const double coverage_ratio =
            coverage < 0 ? 1 :
            coverage / std::accumulate(coverages.begin(), coverages.end(), 0.0) 
            *  coverages.size();
        
        size_t hmr_idx = 0;
        for (size_t i = 0; i < simulated_cpgs.size(); ++i)
        {
            while (hmr_idx < hmrs.size()
                   && ((hmrs[hmr_idx].get_chrom() < simulated_cpgs[i].get_chrom())
                       || (hmrs[hmr_idx].get_chrom() == simulated_cpgs[i].get_chrom()
                           && hmrs[hmr_idx].get_end() < simulated_cpgs[i].get_start())))
                ++hmr_idx;

            const string old_name = simulated_cpgs[i].get_name();
            const size_t count = read_count(coverages[i], coverage_ratio);
            
            simulated_cpgs[i].set_name(old_name.substr(0, old_name.find(":") + 1)
                                       + smithlab::toa(count));

            if (hmrs[hmr_idx].contains(simulated_cpgs[i]))
            {
                const double meth = noisy_meth(count, hypo_meth_distro.sample());
                simulated_cpgs[i].set_score(meth);
            }
            else
            {
                const double meth = noisy_meth(count, hyper_meth_distro.sample());
                simulated_cpgs[i].set_score(meth);
            }
        }

        std::ostream *out = outfile.empty() ?
            &std::cout : new std::ofstream(outfile.c_str());
        for (size_t i = 0; i < simulated_cpgs.size(); ++i)
            *out << simulated_cpgs[i] << endl;

        if (!outfile.empty()) delete out;
        if (VERBOSE)
            cerr << "Done" << endl;
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

