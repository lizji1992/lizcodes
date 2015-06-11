/* Copyright (C) 2009 University of Southern California
 *                    Andrew D Smith
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

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "TwoStateCTHMM.hpp"
#include "Distro.hpp"
#include "false_discovery_rate.hpp"


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

static void
load_cpgs(const bool VERBOSE, 
          const string &cpgs_file, vector<SimpleGenomicRegion> &cpgs,
          vector<pair<double, double> > &meth, vector<size_t> &reads,
          vector<size_t> &mytime) 
{
    if (VERBOSE)
        cerr << "[READING CPGS AND METH PROPS]" << endl;

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
    cpgs.push_back(GenomicRegion(chrom, pos, pos+1, seq, 0, strand[0]));       
    reads.push_back(coverage);
    meth.push_back(std::make_pair(0.0, 0.0));
    meth.back().first = static_cast<size_t>(round(level * coverage));
    meth.back().second = static_cast<size_t>(coverage  - meth.back().first); 
    mytime.push_back(cpgs.back().get_start());
    }
    if (VERBOSE)
      cerr << "TOTAL CPGS: " << cpgs.size() << endl
           << "MEAN COVERAGE: " 
           << accumulate(reads.begin(), reads.end(), 0.0)/reads.size() << endl
           << endl;
}

template <class T, class U, class C> static void
separate_regions(const bool VERBOSE, const size_t desert_size, 
                 vector<SimpleGenomicRegion> &cpgs,
                 vector<T> &meth, vector<U> &reads,
                 vector<C> &mytime,
                 vector<size_t> &reset_points) 
{
    if (VERBOSE)
        cerr << "[SEPARATING BY CPG DESERT]" << endl;
    // eliminate the zero-read cpgs
    size_t j = 0;
    for (size_t i = 0; i < cpgs.size(); ++i)
        if (reads[i] > 0) 
        {
            cpgs[j] = cpgs[i];
            meth[j] = meth[i];
            reads[j] = reads[i];
            mytime[j] = mytime[i];
            ++j;
        }
    cpgs.erase(cpgs.begin() + j, cpgs.end());
    meth.erase(meth.begin() + j, meth.end());
    reads.erase(reads.begin() + j, reads.end());
    mytime.erase(mytime.begin() + j, mytime.end());
  
    // segregate cpgs
    size_t prev_cpg = 0;
    for (size_t i = 0; i < cpgs.size(); ++i) 
    {
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
read_params_file(const string &params_file, 
                 betabin &fg_emission,
                 betabin &bg_emission,
                 Distro &fg_duration,
                 Distro &bg_duration) 
{
    std::ifstream in(params_file.c_str());
    string fg_emission_str, bg_emission_str, fg_duration_str, bg_duration_str;
    std::getline(in, fg_emission_str);
    std::getline(in, bg_emission_str);
    std::getline(in, fg_duration_str);
    std::getline(in, bg_duration_str);

    fg_emission = betabin(fg_emission_str);
    bg_emission = betabin(bg_emission_str);
    fg_duration = Distro(fg_duration_str);
    bg_duration = Distro(bg_duration_str);
}

static void
pick_sample(const vector<pair<double, double> > &meth,
            const vector<size_t> &mytime,
            const vector<size_t> &reset_points,
            const size_t training_size,
            vector<pair<double, double> > &meth_sample,
            vector<size_t> &mytime_sample,
            vector<size_t> &reset_points_sample)
{
    // random training sample
    vector<size_t> idxs(reset_points.size() - 1);
    for (size_t i = 0; i < idxs.size(); ++i) idxs[i] = i;
    srand(time(0) + getpid());
    std::random_shuffle(idxs.begin(), idxs.end());

    size_t sample_size = 0;
    size_t i = 0;
    reset_points_sample.push_back(sample_size);
    while (i < idxs.size() && sample_size < training_size)
    {	
        const size_t idx = idxs[i];
        const size_t start = reset_points[idx];
        const size_t end = reset_points[idx + 1];

        std::copy(meth.begin() + start, meth.begin() + end,
                  std::back_inserter(meth_sample));
        std::copy(mytime.begin() + start, mytime.begin() + end,
                  std::back_inserter(mytime_sample));
        sample_size += end - start;
        reset_points_sample.push_back(sample_size);
        ++i;
    }
    assert(meth_sample.size() == reset_points_sample.back());
}

static void
build_domains(const bool VERBOSE, 
              const vector<SimpleGenomicRegion> &cpgs,
              const vector<double> &post_scores,
              const vector<size_t> &reset_points,
              const vector<bool> &classes,
              vector<GenomicRegion> &domains) 
{
    static const bool CLASS_ID = true;
    size_t n_cpgs = 0, n_domains = 0, reset_idx = 1, prev_end = 0;
    bool in_domain = false;
    double score = 0;
    for (size_t i = 0; i < classes.size(); ++i) 
    {
        if (reset_points[reset_idx] == i) 
        {
            if (in_domain) 
            {
                in_domain = false;
                domains.back().set_end(prev_end);
                domains.back().set_score(score);
                n_cpgs = 0;
                score = 0;
            }
            ++reset_idx;
        }
        if (classes[i] == CLASS_ID) 
        {
            if (!in_domain) 
            {
                in_domain = true;
                domains.push_back(GenomicRegion(cpgs[i]));
                domains.back().set_name("HYPO" + toa(n_domains++));
            }
            ++n_cpgs;
            score += post_scores[i];
        }
        else if (in_domain) 
        {
            in_domain = false;
            domains.back().set_end(prev_end);
            domains.back().set_score(score);//n_cpgs);
            n_cpgs = 0;
            score = 0;
        }
        prev_end = cpgs[i].get_end();
    }
    // Do we miss the final domain?????
}

static void
get_domain_scores(const vector<bool> &classes,
		  const vector<pair<double, double> > &meth,
		  const vector<size_t> &reset_points,
		  vector<double> &scores) {
  static const bool CLASS_ID = true;
  size_t n_cpgs = 0, reset_idx = 1;
  bool in_domain = false;
  double score = 0;
  for (size_t i = 0; i < classes.size(); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
	in_domain = false;
	scores.push_back(score);
	score = 0;
      }
      ++reset_idx;
    }
    if (classes[i] == CLASS_ID) {
      in_domain = true;
      score += 1.0 - (meth[i].first/(meth[i].first + meth[i].second));
      ++n_cpgs;
    }
    else if (in_domain) {
      in_domain = false;
      scores.push_back(score);
      score = 0;
    }
  }
}

static void
calculate_random_scores_by_shuffle_cpgs(
    const vector<pair<double, double> > &meth, 
    const vector<size_t> &mytime, 
    const vector<size_t> &reset_points, 
    const size_t training_size,
    const betabin &fg_emission, const betabin &bg_emission,
    const Distro &fg_duration, const Distro &bg_duration,
    const size_t MAX_LEN, const double min_prob,
    const double tolerance, const size_t max_iterations,
    const bool VERBOSE,
    vector<double> &random_scores) 
{
    // get the sample used to compute background domain scores
    vector<pair<double, double> > meth_sample;
    vector<size_t> mytime_sample;
    vector<size_t> reset_points_sample;
    pick_sample(meth, mytime, reset_points,
                training_size > 0 ? training_size : meth.size(),
                meth_sample, mytime_sample, reset_points_sample);
    srand(time(0) + getpid());
    random_shuffle(meth_sample.begin(), meth_sample.end());

    TwoStateCTHMM hmm(meth_sample, mytime_sample, reset_points_sample,
                      MAX_LEN, min_prob, tolerance, max_iterations, VERBOSE);

    hmm.set_parameters(fg_emission, bg_emission, fg_duration, bg_duration);
    hmm.PosteriorDecoding();

    vector<bool> classes;
    vector<double> scores;
    hmm.get_posterior_scores(scores, classes);

    random_scores.clear();
    get_domain_scores(classes, meth, reset_points, random_scores);
    std::sort(random_scores.begin(), random_scores.end());
}

int
main(int argc, const char **argv) 
{
    try 
    {
        string outfile;
        string scores_file;
    
        size_t desert_size = 1000;
        size_t max_iterations = 10;
        size_t training_size = 0;
        
        // corrections for small values (not parameters):
        double tolerance = 1e-10;
        double min_prob  = 1e-10;
        size_t MAX_LEN = 200;
        double fdr = 0.05;
        double fdr_cutoff = std::numeric_limits<double>::max();

        string params_in_file;
        string params_out_file;

        // run mode flags
        bool VERBOSE = false;
    
        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(argv[0], "A program for segmenting DNA "
                               "methylation data"
                               "<cpg-BED-file>");
        opt_parse.add_opt("out", 'o', "output file (BED format)", 
                          false, outfile);
        opt_parse.add_opt("scores", 's', "scores file (WIG format)", 
                          false, scores_file);

        opt_parse.add_opt("desert", 'd', "desert size", false, desert_size);
        opt_parse.add_opt("training-size", '\0',
                          "Maximum number of data points for HMM training",
                          OptionParser::OPTIONAL, training_size);
        opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations); 
        opt_parse.add_opt("max-len", 'L', "max foreground length", false, MAX_LEN); 
        opt_parse.add_opt("fdr", 'F', "False discovery rate (default 0.05)",
                          OptionParser::OPTIONAL, fdr); 
        opt_parse.add_opt("fdr-cutoff", '\0',
                          "P-value cutoff based on false discovery rate",
                          OptionParser::OPTIONAL, fdr_cutoff); 

        opt_parse.add_opt("params-in", 'P', "HMM parameters file", false, params_in_file);
        opt_parse.add_opt("params-out", 'p', "HMM parameters file", false, params_out_file);
        opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    
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
        if (leftover_args.empty()) 
        {
            cerr << opt_parse.help_message() << endl;
            return EXIT_SUCCESS;
        }
        const string cpgs_file = leftover_args.front();
        /****************** END COMMAND LINE OPTIONS *****************/
    
        // separate the regions by chrom and by desert
        vector<SimpleGenomicRegion> cpgs;
        // vector<double> meth;
        vector<pair<double, double> > meth;
        vector<size_t> mytime;
        vector<size_t> reads;
        load_cpgs(VERBOSE, cpgs_file, cpgs, meth, reads, mytime);
  
        // separate the regions by chrom and by desert, and eliminate
        // those isolated CpGs
        vector<size_t> reset_points;
        separate_regions(VERBOSE, desert_size, cpgs,
                         meth, reads, mytime, reset_points);

        TwoStateCTHMM hmm(meth, mytime, reset_points, MAX_LEN,
                          min_prob, tolerance, max_iterations, VERBOSE);

        /***********************************
         * STEP: Parameter initializing and model training
         */
        betabin fg_emission, bg_emission;
        Distro fg_duration, bg_duration;
        if (!params_in_file.empty()) 
        {
            read_params_file(params_in_file, fg_emission, bg_emission,
                             fg_duration, bg_duration);
        }
        else 
        {
            const double n_reads =
                accumulate(reads.begin(), reads.end(), 0.0)/reads.size();
            double fg_alpha = 0.33*n_reads;
            double fg_beta = 0.67*n_reads;
            double bg_alpha = 0.67*n_reads;
            double bg_beta = 0.33*n_reads;
            
            fg_emission = betabin(fg_alpha, fg_beta);
            bg_emission = betabin(bg_alpha, bg_beta);

            fg_duration = Distro("exp 50");
            bg_duration = Distro("exp 500");
        }

        // model training
        if (training_size > 0)
        {
            // train with part of the dataset
            vector<pair<double, double> >  meth_sample;
            vector<size_t>  mytime_sample;
            vector<size_t> reset_points_sample;
            pick_sample(meth, mytime, reset_points, training_size,
                        meth_sample, mytime_sample, reset_points_sample);
            TwoStateCTHMM hmm_training(
                meth_sample, mytime_sample, reset_points_sample, MAX_LEN,
                min_prob, tolerance, max_iterations, VERBOSE);

            hmm_training.set_parameters(fg_emission, bg_emission,
                                        fg_duration, bg_duration);
            hmm_training.BaumWelchTraining();
            hmm_training.get_parameters(fg_emission, bg_emission,
                                        fg_duration, bg_duration);
            hmm.set_parameters(fg_emission, bg_emission,
                               fg_duration, bg_duration);
        }
        else
        {
            // train with the whole dataset
            hmm.set_parameters(fg_emission, bg_emission,
                               fg_duration, bg_duration);
            if (max_iterations > 0) 
                hmm.BaumWelchTraining();
            hmm.get_parameters(fg_emission, bg_emission,
                               fg_duration, bg_duration);
        }
      
        /***********************************
         * STEP 5: DECODE THE DOMAINS
         */
        vector<bool> classes;
        vector<double> scores;
        hmm.PosteriorDecoding();
        hmm.get_posterior_scores(scores, classes);
    
        // get domains and domain scores
        vector<double> domain_scores;
        get_domain_scores(classes, meth, reset_points, domain_scores);

        vector<GenomicRegion> domains;
        build_domains(VERBOSE, cpgs, scores, reset_points, classes, domains);

        // calculate cutoff
        vector<double> random_scores;
        calculate_random_scores_by_shuffle_cpgs(
            meth, mytime, reset_points, training_size, fg_emission, bg_emission,
            fg_duration, bg_duration, MAX_LEN, min_prob, tolerance,
            max_iterations, VERBOSE, random_scores);
        
        vector<double> p_values;
        FDR::assign_empirical_p_values(random_scores, domain_scores, p_values);
        if (fdr_cutoff == numeric_limits<double>::max())
            fdr_cutoff = FDR::get_fdr_cutoff(p_values, fdr);
        
        // filtering domains
        size_t j = 0;
        for (size_t i = 0; i < domains.size(); ++i)
            if (p_values[i] <= fdr_cutoff)
            {
                domains[j] = domains[i];
                domain_scores[j] = domain_scores[i];
                ++j;
            }
        domains.erase(domains.begin() + j, domains.end());
        domain_scores.erase(domain_scores.begin() + j, domain_scores.end());
      
        /***********************************
         * STEP 6: WRITE THE RESULTS
         */
        std::ostream *out = (outfile.empty()) ? &cout : 
            new std::ofstream(outfile.c_str());

        for (size_t i = 0; i < domains.size(); ++i) 
        {
            domains[i].set_score(domain_scores[i]);
            *out << domains[i] << endl;
        }
        if (out != &cout) delete out;
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


