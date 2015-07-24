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

int
main(int argc, const char **argv) 
{
    try 
    {
        static const bool REQUIRED = true;
        static const bool OPTIONAL = false;

        string outfile;
        string segment_file;
        string hmrloci_file;
        string nonhmrloci_file;

        string hmr_len_str = "nbd 20 0.2";
        string nonhmr_len_str = "nbd 100 0.5";
        
        size_t hmrn = 10;
        string chrome_name = "chrZ";
        bool KEEP_ORIGAL_COORDINATES = false;

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
        opt_parse.add_opt("hmrloci_file", '\0', "File of CpG loci inside HMRs", 
                          REQUIRED, hmrloci_file);
        opt_parse.add_opt("nonhmrloci_file", '\0', "File of CpG locis outside HMRs", 
                          REQUIRED, nonhmrloci_file);
        opt_parse.add_opt("hmr_len_distro", '\0', "Length distribution of HMRs", 
                          OPTIONAL, hmr_len_str);
        opt_parse.add_opt("nonhmr_len_distro", '\0', "Length distribution of non-HMRs", 
                          OPTIONAL, nonhmr_len_str);
        opt_parse.add_opt("verbose", 'v', "print more run info",
                          OPTIONAL, VERBOSE);
        opt_parse.add_opt("keep-original-coordinate", 'k', "Keep orginal coordinates",
                          OPTIONAL, KEEP_ORIGAL_COORDINATES);
    
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
    
        // Read in genomic locis inside and outside of HMRs
        vector<GenomicRegion> hmrloci;
        vector<GenomicRegion> nonhmrloci;
        ReadBEDFile(hmrloci_file, hmrloci);
        ReadBEDFile(nonhmrloci_file, nonhmrloci);


        // region length distribution
        Distro hmr_len_distro(hmr_len_str);
        Distro nonhmr_len_distro(nonhmr_len_str);

        // simulation; always assume the first and the last region are
        // non-HMRs
        std::ostream *out = outfile.empty() ?
            &std::cout : new std::ofstream(outfile.c_str());
        std::ostream *out_segments = segment_file.empty() ?
            0 : new std::ofstream(segment_file.c_str());
        size_t nloci;
        
        // sample the first region
        nloci = static_cast<size_t>(nonhmr_len_distro.sample());
        std::random_shuffle(nonhmrloci.begin(), nonhmrloci.end());
        vector<GenomicRegion> cpgs;
        std::copy(nonhmrloci.begin(), nonhmrloci.begin() + nloci,
                  std::back_inserter(cpgs));
        if (out_segments) *out_segments << nloci << endl;

        for (size_t i = 0; i < hmrn; ++i)
        {
            // sample HMR
            nloci = static_cast<size_t>(hmr_len_distro.sample());
            std::random_shuffle(hmrloci.begin(), hmrloci.end());
            std::copy(hmrloci.begin(), hmrloci.begin() + nloci,
                      std::back_inserter(cpgs)); 
            if (out_segments) *out_segments << nloci << endl;

            // sample non-HMR
            nloci = static_cast<size_t>(nonhmr_len_distro.sample());
            std::random_shuffle(nonhmrloci.begin(), nonhmrloci.end());
            std::copy(nonhmrloci.begin(), nonhmrloci.begin() + nloci,
                      std::back_inserter(cpgs)); 
            if (out_segments) *out_segments << nloci << endl;
        }

        // normalize locations
        if (!KEEP_ORIGAL_COORDINATES)
        {
            size_t loc = 1;
            for (size_t i = 0; i < cpgs.size(); ++i)
            {
                cpgs[i].set_chrom(chrome_name);
                cpgs[i].set_start(loc);
                cpgs[i].set_end(loc + 1);
                loc += 2;
            }
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

