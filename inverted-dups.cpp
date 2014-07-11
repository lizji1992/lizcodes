/*    inverted-dups: a program for doing quality control relative to
 *    the inverted duplicates reads problem and masking "bad" reads if needed.
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith, Liz Ji and Jenny Qu
 *
 *    Authors: Andrew D. Smith, Liz Ji and Jenny Qu
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <sstream>
#include <algorithm>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::min;
using std::accumulate;
using std::tr1::unordered_map;
using std::ptr_fun;


// store each read from one end
struct FASTQRecord {
  string name;
  string seq;
  string seqtag;
  string score;

  string tostring() const {
    string tmp;
    std::ostringstream s;
    s << '@' << name << '\n' << seq << '\n' << tmp << '\n' << score;  
    return s.str();
  }
    
  void revcomp();
  
  // see if two reads from two ends match to each other
  // (they should have the same name)
  static bool 
  mates(const size_t to_ignore_at_end, // in case names have #0/1 name ends
  const FASTQRecord &a, const FASTQRecord &b) {
    const string namefield1 = a.name.substr(0, a.name.find_first_of(' '));
    const string namefield2 = b.name.substr(0, b.name.find_first_of(' '));

    const bool same_name = 
               equal(namefield1.begin(), 
               namefield1.begin() + namefield1.length() - to_ignore_at_end,
               namefield2.begin());
    return same_name && a.seq.length() == b.seq.length();
  }
  
};

// Reverse the sequence
void
FASTQRecord::revcomp() {
  string tmp(seq);
  reverse(tmp.begin(), tmp.end());
  revcomp_inplace(seq);
  for (size_t i = 0; i < seq.length(); ++i)
    if (tmp[i] == '_') seq[i] = '_';
  reverse(score.begin(), score.end());
}

std::ostream& 
operator<<(std::ostream& s, const FASTQRecord &r) {
  return s << r.tostring();
}

// Read 4 lines one time from fastq and fill in the FASTQRecord structure
std::istream& 
operator>>(std::istream& s, FASTQRecord &r) {
  if (getline(s, r.name)) {
    if (r.name.empty() || r.name[0] != '@') 
      throw SMITHLABException("FASTQ file out of sync at '@'");
    
    if (!getline(s, r.seq))
      throw SMITHLABException("FASTQ file truncated expecting seq");
    
    if (!getline(s, r.seqtag))
      throw SMITHLABException("FASTQ file truncated expecting '+' line");
    
    if (r.seqtag.empty() || r.seqtag[0] != '+') 
      throw SMITHLABException("FASTQ file out of sync [missing '+']");
    
    if (!getline(s, r.score))
      throw SMITHLABException("FASTQ file truncated expecting score");
  }
  else s.setstate(std::ios::badbit);

  return s;
}

////////////////////////////////////////////////////////////////////////////////

static bool
similar_letters_bisulfite(const char a, const char b) {
  return (a == b) || (a == 'T' && b == 'C') || (a == 'G' && b == 'A');
}

static bool
similar_letters_bisulfite_TCwildcard(const char a, const char b) {
  return (a == b) || (a == 'T' && b == 'C');
}

static bool
similar_letters_bisulfite_GAwildcard(const char a, const char b) {
  return (a == b) || (a == 'G' && b == 'A');
}


// Compare two reads to detect the overlapped region
vector<size_t>
invdup_similarity(FASTQRecord &r1, FASTQRecord &r2,
                  string &brick1, string &brick2,
                  vector<size_t> &pos_count_overlap) {
  vector<size_t> sim (2, 0);
  size_t sim_TC = 0;
  size_t sim_GA = 0;

  if (pos_count_overlap.empty()) pos_count_overlap.resize(r1.seq.size());
  
  string::const_iterator it1(r1.seq.begin());
  string::const_iterator it2(r2.seq.begin());
  vector<size_t>::iterator it_pos(pos_count_overlap.begin());
  while (it1 < r1.seq.end() && it2 < r2.seq.end()) {
    if (similar_letters_bisulfite(*it1, *it2)) {
      sim[0]++;
      *it_pos += 1;
      if (*it1 == *it2) {
        brick1 += '*'; brick2 += '*';
      } else {
        brick1 += *it1; brick2 += *it2;
      }
    } else {
      brick1 += '-'; brick2 += '-';
    }
    sim_TC += similar_letters_bisulfite_TCwildcard(*it1, *it2);
    sim_GA += similar_letters_bisulfite_GAwildcard(*it1, *it2);
    it1++; it2++; it_pos++;
  }
  sim[1] = max(sim_TC, sim_GA);
  return sim;
}

int 
main(int argc, const char **argv) {

  try {

    string fp_repo;
    string fp_stat;
    string fp_proc_fq;
    double cutoff = 0.95;
    size_t to_ignore_at_end_of_name = 0;
    bool VERBOSE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "count the invdup reads "
                           "in the input files", "<end1-fastq> <end2-fastq>");
    opt_parse.add_opt("output", 'o', 
                      "Name of the scanning results (default: stdout)", 
                      false, fp_repo);
    opt_parse.add_opt("stat", 's',
                      "Name of the output stats file (default: stdout)", 
                      false, fp_stat);
    opt_parse.add_opt("masking", 'm',
                      "Name of the new second-end fastq file "
                      "if you want to mask the invdup reads", 
                      false, fp_proc_fq);
    opt_parse.add_opt("cutoff", 'c',
                      "The cutoff for invdup reads (default: 0.95)",
                      false, cutoff);
    opt_parse.add_opt("ignore", 'i', "Ignore this number of letters "
                      "at end of name", false, to_ignore_at_end_of_name);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested() || leftover_args.size() != 2) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    const string reads_file_one = leftover_args.front();
    const string reads_file_two = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    // Input: paired-end reads with end1 and end2
    std::ifstream if_read1(reads_file_one.c_str());
    if(!if_read1)
      throw SMITHLABException("cannot open input file " + reads_file_one);
    
    std::ifstream if_read2(reads_file_two.c_str());
    if(!if_read2)
      throw SMITHLABException("cannot open input file " + reads_file_two);

    // Output scanning results: 
    std::ofstream of_repo;
    if (!fp_repo.empty()) of_repo.open(fp_repo.c_str());
    std::ostream os_repo(fp_repo.empty() ? cout.rdbuf() : of_repo.rdbuf());   

    string fp_repo_brick = fp_repo + ".brick";
    std::ofstream of_repo_brick(fp_repo_brick.c_str());
    
    // Output the proccessed fastq files:
    std::ofstream of_proc_fq(fp_proc_fq.c_str());
    if(!fp_proc_fq.empty() && !of_proc_fq)
      throw SMITHLABException("cannot write new fastq file " + fp_proc_fq);

    //------------SCAN THE READS------------//
    size_t num_read = 0;
    vector<size_t> num_bad_read (2, 0);
    vector<double> sum_percent_overlap (2, 0);
    vector<double> sum_bad_percent_overlap (2, 0);
    vector<size_t> pos_count_overlap;
    FASTQRecord end_one, end_two;

    while (if_read1 >> end_one && if_read2 >> end_two) {
      
      // Two reads should be in paired-ends
      if (!FASTQRecord::mates(to_ignore_at_end_of_name, end_one, end_two)) 
        throw SMITHLABException("expected mates, got:" + 
              end_one.tostring() + "\n" + end_two.tostring());
      
      // See if inverted duplicates emerge
      string brick_end1;
      string brick_end2;
      const vector<size_t> sim = invdup_similarity(end_one, end_two,
                                 brick_end1, brick_end2, pos_count_overlap);
      vector<double> percent_overlap (2, 0);
      // percent_overlap[0] = two wildcards are allowed at the same time
      // percent_overlap[1] = two wildcards are not allowed at the same time
      percent_overlap[0] = static_cast<double>(sim[0])/end_one.seq.length();
      percent_overlap[1] = static_cast<double>(sim[1])/end_one.seq.length();
      
      if (percent_overlap[0] > cutoff) {
        num_bad_read[0]++;
        sum_bad_percent_overlap[0] += percent_overlap[0];
      }
      if (percent_overlap[1] > cutoff) {
        num_bad_read[1]++;
        sum_bad_percent_overlap[1] += percent_overlap[1];
      }

      // Write new fastq file if -n 
      if (of_proc_fq) {
        const string masked_seq(end_two.seq.length(), 'N');
        of_proc_fq << end_two.name << "\n";
        if (percent_overlap[0] > cutoff) 
          of_proc_fq << masked_seq << "\n";
        else {
          of_proc_fq << end_two.seq << "\n";
          of_proc_fq << end_two.seqtag << "\n"
                    << end_two.score << "\n";
        }
      }

      num_read++;
      sum_percent_overlap[0] += percent_overlap[0];
      sum_percent_overlap[1] += percent_overlap[1];
      os_repo << sim[0] << "," << sim[1] << '\t' 
              << percent_overlap[0] << "," << percent_overlap[1] << "\n";
      if (percent_overlap[0] > cutoff) 
        of_repo_brick << brick_end1 << "\n"
                      << brick_end2 << "\n" << endl;
    }
    //------------WRITE STAT INFORMATION------------//
    std::ofstream of_stat;
    if (!fp_stat.empty()) of_stat.open(fp_stat.c_str());
    std::ostream os_stat(fp_stat.empty() ? cout.rdbuf() : of_stat.rdbuf()); 
    os_stat << "CUTOFF:\t" << cutoff << "\n" 
            << "TOTAL READ PAIRS:\t" << num_read << "\n"  
            << "SUSPECT INVERTED-DUPLICATED READ PAIRS:\t" 
            << num_bad_read[0] << "," << num_bad_read[1] << "\n"  
            << "PERCENTAGE OF GOOD READS:\t" 
            << 1 - static_cast<double>(num_bad_read[0])/num_read << ","  
            << 1 - static_cast<double>(num_bad_read[1])/num_read << "\n"  
            << "MEAN OVERLAP PERCENTAGE:\t" 
            << sum_percent_overlap[0]/num_read << "," 
            << sum_percent_overlap[1]/num_read << "\n" 
            << "MEAN OVERLAP PERCENTAGE OF INVERTED-DUPLICATES:\t" 
            << sum_bad_percent_overlap[0]/num_bad_read[0] << ","
            << sum_bad_percent_overlap[1]/num_bad_read[1] << "\n"
            << endl;

    size_t count_pos = 0;
    vector<size_t>::iterator it(pos_count_overlap.begin());
    while (it < pos_count_overlap.end()) {
      os_stat << count_pos << ' '
              << static_cast<double>(*it++)/num_read << endl;
      count_pos++;
    }
  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
