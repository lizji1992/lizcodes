/* Copyright (C) 2009 University of Southern California
 *                    Andrew D Smith
 * Author: Andrew D. Smith
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
#include "TwoStateVDHMM.hpp"
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

static void
load_cpgs(const bool VERBOSE,
          const string &cpgs_file, vector<SimpleGenomicRegion> &cpgs,
          vector<pair<double, double> > &meth, vector<size_t> &reads)
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
  }
  if (VERBOSE)
    cerr << "TOTAL CPGS: " << cpgs.size() << endl
    << "MEAN COVERAGE: "
    << accumulate(reads.begin(), reads.end(), 0.0)/reads.size() << endl
    << endl;
}


double
get_fdr_cutoff(const vector<double> &scores, const double fdr) 
{
    if (fdr <= 0)
        return numeric_limits<double>::max();
    else if (fdr > 1)
        return numeric_limits<double>::min();
    vector<double> local(scores);
    std::sort(local.begin(), local.end());
    size_t i = 0;
    for (; i < local.size() - 1 && 
             local[i+1] < fdr*static_cast<double>(i+1)/local.size(); ++i);
    return local[i];
}


static void
write_scores_bedgraph(const string &filename,
                      const vector<SimpleGenomicRegion> &cpgs,
                      const vector<double> &scores) 
{
    std::ofstream wigout(filename.c_str());
    for (size_t i = 0; i < cpgs.size(); ++i)
        wigout << cpgs[i] << "\t" << scores[i] << "\n";
    wigout.close();
}

static void
get_domain_scores(const vector<bool> &classes,
                  const vector<pair<double, double> > &meth,
                  const vector<size_t> &reset_points,
                  vector<double> &scores) 
{
    static const bool CLASS_ID = true;
    size_t n_cpgs = 0, reset_idx = 1;
    bool in_domain = false;
    double score = 0;
    for (size_t i = 0; i < classes.size(); ++i) 
    {
        if (reset_points[reset_idx] == i) 
        {
            if (in_domain) 
            {
                in_domain = false;
                scores.push_back(score);
                score = 0;
            }
            ++reset_idx;
        }
        if (classes[i] == CLASS_ID) 
        {
            in_domain = true;
            score += 1.0 - (meth[i].first/(meth[i].first + meth[i].second));
            ++n_cpgs;
        }
        else if (in_domain) 
        {
            in_domain = false;
            scores.push_back(score);
            score = 0;
        }
    }
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


template <class T, class U> static void
separate_regions(const bool VERBOSE, const size_t desert_size, 
                 vector<SimpleGenomicRegion> &cpgs,
                 vector<T> &meth, vector<U> &reads,
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
            ++j;
        }
    cpgs.erase(cpgs.begin() + j, cpgs.end());
    meth.erase(meth.begin() + j, meth.end());
    reads.erase(reads.begin() + j, reads.end());
  
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


// static void
// shuffle_cpgs(const TwoStateHMMB &hmm,
//              vector<pair<double, double> > meth, 
//              vector<size_t> reset_points, 
//              const vector<double> &start_trans,
//              const vector<vector<double> > &trans,
//              const vector<double> &end_trans,
//              const double fg_alpha, const double fg_beta, 
//              const double bg_alpha, const double bg_beta,
//              vector<double> &domain_scores) 
// {
//     srand(time(0) + getpid());
//     random_shuffle(meth.begin(), meth.end());
//     vector<bool> classes;
//     vector<double> scores;
//     hmm.PosteriorDecoding(meth, reset_points, start_trans, trans, end_trans,
//                           fg_alpha, fg_beta, bg_alpha, bg_beta, classes, scores);
//     random_shuffle(meth.begin(), meth.end());
//     classes.clear();
//     scores.clear();
//     hmm.PosteriorDecoding(meth, reset_points, start_trans, trans, end_trans,
//                           fg_alpha, fg_beta, bg_alpha, bg_beta, classes, scores);
//     get_domain_scores(classes, meth, reset_points, domain_scores);
//     sort(domain_scores.begin(), domain_scores.end());
// }


// static void
// assign_p_values(const vector<double> &random_scores, 
//                 const vector<double> &observed_scores, 
//                 vector<double> &p_values) 
// {
//     const double n_randoms =
//         random_scores.size() == 0 ? 1 : random_scores.size();
//     for (size_t i = 0; i < observed_scores.size(); ++i)
//         p_values.push_back((random_scores.end() - 
//                             upper_bound(random_scores.begin(),
//                                         random_scores.end(),
//                                         observed_scores[i]))/n_randoms);
// }

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

    Distro fg(fg_emission_str);
    fg_emission = betabin(fg.get_params().front(), fg.get_params().back());

    Distro bg(bg_emission_str);
    bg_emission = betabin(bg.get_params().front(), bg.get_params().back());

    fg_duration = Distro(fg_duration_str);
    bg_duration = Distro(bg_duration_str);
}

// static void
// write_params_file(const string &outfile, 
//                   const double fg_alpha, 
//                   const double fg_beta, 
//                   const double bg_alpha, 
//                   const double bg_beta,
//                   const vector<double> &start_trans, 
//                   const vector<vector<double> > &trans, 
//                   const vector<double> &end_trans) 
// {
//     std::ostream *out = (outfile.empty()) ? &cerr : 
//         new std::ofstream(outfile.c_str());
//     out->precision(30);
//     *out << "FG_ALPHA\t" << fg_alpha << endl
//          << "FG_BETA\t" << fg_beta << endl
//          << "BG_ALPHA\t" << bg_alpha << endl
//          << "BG_BETA\t" << bg_beta << endl
//          << "S_F\t" << start_trans[0] << endl
//          << "S_B\t" << start_trans[1] << endl
//          << "F_F\t" << trans[0][0] << endl
//          << "F_B\t" << trans[0][1] << endl
//          << "B_F\t" << trans[1][0] << endl
//          << "B_B\t" << trans[1][1] << endl
//          << "F_E\t" << end_trans[0] << endl
//          << "B_E\t" << end_trans[1] << endl
//         ;
//     if (out != &cerr) delete out;
// }


int
main(int argc, const char **argv) 
{

    try 
    {

        string outfile;
        string scores_file;
        string trans_file;
        string dataset_name;
    
        size_t desert_size = 1000;
        size_t max_iterations = 10;
    
        // run mode flags
        bool VERBOSE = false;
    
        // corrections for small values (not parameters):
        double tolerance = 1e-10;
        double min_prob  = 1e-10;
        size_t MAX_LEN = 200;

        string params_in_file;
        string params_out_file;
    
        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(argv[0], "A program for segmenting DNA "
                               "methylation data"
                               "<cpg-BED-file>");
        opt_parse.add_opt("out", 'o', "output file (BED format)", 
                          false, outfile);
        opt_parse.add_opt("scores", 's', "scores file (WIG format)", 
                          false, scores_file);
        opt_parse.add_opt("tolerance", 't', "Tolerance", 
                          false, tolerance);
        opt_parse.add_opt("desert", 'd', "desert size", false, desert_size);
        opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations); 
        opt_parse.add_opt("max-len", 'L', "max foreground length", false, MAX_LEN); 
        opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
        opt_parse.add_opt("name", 'N', "data set name", false, dataset_name);
    
        opt_parse.add_opt("params-in", 'P', "HMM parameters file", false, params_in_file);
        opt_parse.add_opt("params-out", 'p', "HMM parameters file", false, params_out_file);
    
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
        vector<size_t> reads;
        load_cpgs(VERBOSE, cpgs_file, cpgs, meth, reads);
    
        // separate the regions by chrom and by desert, and eliminate
        // those isolated CpGs
        vector<size_t> reset_points;
        separate_regions(VERBOSE, desert_size, cpgs, meth, reads, reset_points);
    
        TwoStateVDHMM hmm(min_prob, tolerance, max_iterations, VERBOSE, MAX_LEN);
    
        betabin fg_emission, bg_emission;
        Distro fg_duration, bg_duration;
    
        // double fdr_cutoff = std::numeric_limits<double>::max();

        if (!params_in_file.empty()) 
        {
            read_params_file(params_in_file,
                             fg_emission, bg_emission,
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

            fg_duration = Distro("nbd 30 0.5");
            bg_duration = Distro("geo 0.005");
        }
    
        hmm.init(meth, reset_points);
        hmm.initialize_parameters(fg_emission, bg_emission,
                                  fg_duration, bg_duration);

        if (max_iterations > 0)
            hmm.BaumWelchTraining();
    
        // if (!params_out_file.empty()) 
        // {
        //     // WRITE ALL THE HMM PARAMETERS:
        //     write_params_file(params_out_file, fg_alpha, fg_beta, bg_alpha, bg_beta,
        //                       start_trans, trans, end_trans);
        // }
    
        /***********************************
         * STEP 5: DECODE THE DOMAINS
         */
        vector<bool> classes;
        vector<double> scores;
        hmm.get_posterior_scores(scores, classes);
    
        vector<double> domain_scores;
        get_domain_scores(classes, meth, reset_points, domain_scores);
    
        // // vector<double> random_scores;
        // // shuffle_cpgs(hmm, meth, reset_points, start_trans, trans, end_trans,
        // //              fg_alpha, fg_beta, bg_alpha, bg_beta, random_scores);
    
        // vector<double> p_values;
        // assign_p_values(random_scores, domain_scores, p_values);
    
        // if (fdr_cutoff == numeric_limits<double>::max())
        //     fdr_cutoff = get_fdr_cutoff(p_values, 0.01);

        // if (!params_out_file.empty()) 
        // {
        //     std::ofstream out(params_out_file.c_str(), std::ios::app);
        //     out.precision(30);
        //     out << "FDR_CUTOFF\t" << fdr_cutoff << endl;
        //     out.close();
        // }
    
        /***********************************
         * STEP 6: WRITE THE RESULTS
         */
        if (!scores_file.empty())
            write_scores_bedgraph(scores_file, cpgs, scores);
        // if (!trans_file.empty()) 
        // {
        //     vector<double> fg_to_bg_scores;
        //     hmm.TransitionPosteriors(meth, reset_points, start_trans, trans, end_trans, 
        //                              fg_alpha, fg_beta, bg_alpha, bg_beta, 
        //                              1, fg_to_bg_scores);
        //     vector<double> bg_to_fg_scores;
        //     hmm.TransitionPosteriors(meth, reset_points, start_trans, trans, end_trans, 
        //                              fg_alpha, fg_beta, bg_alpha, bg_beta, 
        //                              2, bg_to_fg_scores);
        //     for (size_t i = 0; i < fg_to_bg_scores.size(); ++i)
        //         fg_to_bg_scores[i] = max(fg_to_bg_scores[i], bg_to_fg_scores[i]);
        //     write_scores_bedgraph(trans_file, cpgs, fg_to_bg_scores);
        // }
    
        vector<GenomicRegion> domains;
        build_domains(VERBOSE, cpgs, scores,
                      reset_points, classes, domains);
      
        std::ostream *out = (outfile.empty()) ? &cout : 
            new std::ofstream(outfile.c_str());

        for (size_t i = 0; i < domains.size(); ++i) 
        {
            // if (p_values[i] < fdr_cutoff) 
            // {
                domains[i].set_score(domain_scores[i]);
                *out << domains[i] << '\n';
            // }
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


