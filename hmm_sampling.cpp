/* Copyright (C) 2019 University of Southern California
 *                    Xiaojing Ji, Jianghan Qu and Andrew D Smith
 *
 * Author: Andrew D. Smith, Jianghan Qu and Xiaojing Ji
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
#include <string>
#include <stdexcept>
#include <random>
#include <algorithm>
#include <iterator>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"

using std::istream_iterator;
using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::runtime_error;
using std::begin;
using std::end;
using std::pair;
using std::make_pair;
using std::setw;

using std::abs;
using std::max;
using std::min;

//typedef double two_by_two[2][2];
typedef vector<vector<double> > two_by_two;


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////
///////        ##     ##    ########     ####     ##            ######
///////        ##     ##       ##         ##      ##           ##
///////        ##     ##       ##         ##      ##            ######
///////        ##     ##       ##         ##      ##                 ##
///////         #######        ##        ####     ########      ######
////////
static void
report_param_header_for_verbose() {
  cerr << setw(3) << "ITR" << setw(14) << "F PARAM" << setw(14) << "B PARAM"
  << setw(10) << "P_FB" << setw(10) << "P_BF" << setw(15) << "DELTA" << endl;
}

static void
report_params_for_verbose(const size_t i, const double fg_p, const double bg_p,
                          const double p_fb_est, const double p_bf_est,
                          const double delta) {
  cerr.precision(2);
  cerr << setw(3) << i << setw(12) << fg_p << setw(13) << bg_p
  << setw(14) << p_fb_est << setw(10) << p_bf_est << setw(15) << delta << endl;
}

inline double
get_delta(const double a, const double b) {
  return (b - a)/max(abs(a), abs(b));
}

inline double
log_sum_log(const double p, const double q) {
  if (p == 0.0) {return q;}
  else if (q == 0.0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}

template <class T> void
one_minus(T a, const T a_end, T b) {
    while (a != a_end)
        *b++ = 1.0 - *a++;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////
///////        ########     ##    ##
///////        ##     ##    ####  ##
///////        ########     ## ## ##
///////        ##     ##    ##   ###
///////        ########     ##    ##
////////
struct Bernoulli {
  Bernoulli() : p(0.5) {}
  Bernoulli(const double _p) : p(_p) {}
  double operator()(const bool val) const;
  void fit(const vector<bool> &vals, const vector<double> &pp);
  double p;
};

double
Bernoulli::operator()(const bool val) const {
  return val ? p : (1-p);
}

void
Bernoulli::fit(const vector<bool> &vals, const vector<double> &pp) {
  const double epsilon = 0.01;
  const double denom = std::accumulate(begin(pp), end(pp), 0.0);
  const double nom = inner_product(begin(vals), end(vals), begin(pp), 0.0);
  assert(denom > 0);
  
  p = min(max(nom / denom, epsilon), 1 - epsilon);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////
///////        ##     ##    ##     ##    ##     ##
///////        ##     ##    #### ####    #### ####
///////        #########    ## ### ##    ## ### ##
///////        ##     ##    ##     ##    ##     ##
///////        ##     ##    ##     ##    ##     ##
////////

struct TwoStateHMM {
  
  TwoStateHMM(const double tol, const size_t max_itr,
              const bool v, const bool e) :
  tolerance(tol), max_iterations(max_itr), VERBOSE(v), FIX_EMIT(e) {}
  void initialize(const vector<bool> &obs);
  void initialize(const string params_file);
  double single_iteration(const vector<bool> &obs,
                          vector<pair<double, double> > &forward,
                          vector<pair<double, double> > &backward,
                          vector<pair<double, double> > &emit,
                          vector<two_by_two> &joint);
  double BaumWelchTraining(const vector<bool> &obs);
  void StatesSampling(const vector<bool> &obs, vector<bool> &x,
                      std::mt19937 &gen) const;
  
  double tolerance;
  size_t max_iterations;
  bool VERBOSE;
  bool FIX_EMIT;

  double p_fb;
  double p_bf;
  Bernoulli fg_distr;
  Bernoulli bg_distr;
  
  double llh; // log likelihood of observed data
};


void
TwoStateHMM::initialize(const vector<bool> &obs) {
  two_by_two N_ij(2, vector<double>(2, 0.0));
  for (size_t i = 0; i < obs.size() - 1; i++)
    N_ij[obs[i]][obs[i+1]]++;
  
  p_fb = N_ij[1][0] / ( N_ij[1][0] +  N_ij[1][1]);
  p_bf = N_ij[0][1] / ( N_ij[0][0] +  N_ij[0][1]);
  fg_distr.p = 0.8;
  bg_distr.p = 0.1;
}

void
TwoStateHMM::initialize(const string params_file) {
  string jnk;
  std::ifstream in(params_file);
  if (!in)
    throw runtime_error("failed to parse params file: " + params_file);
  in >> jnk >> fg_distr.p >> jnk >> bg_distr.p >> jnk >> p_fb >> jnk >> p_bf;
  
  max_iterations = 0;

  if (VERBOSE)
    cerr << "[LOAD PARAMETERS]" << endl
    << "FG_P\t" << fg_distr.p << endl << "BG_P\t" << bg_distr.p << endl
    << "F_B\t" << p_fb << endl << "B_F\t" << p_bf << endl;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


static void
get_log_emissions(const vector<bool> &v, vector<pair<double, double> > &emit,
                  const Bernoulli &fg_distr, const Bernoulli &bg_distr) {
  for(size_t i = 0; i < v.size(); i++)
    emit[i] = make_pair(log(bg_distr(v[i])), log(fg_distr(v[i])));
}

inline static double
get_posterior(const pair<double, double> &f, const pair<double, double> &b) {
  const double fg = f.second + b.second;
  return exp(fg - log_sum_log(fg, f.first + b.first));
}

inline static void
get_posteriors(const vector<pair<double, double> > &forward,
               const vector<pair<double, double> > &backward,
               vector<double> &posteriors) {
  posteriors.resize(forward.size());
  for (size_t i = 0; i < forward.size(); ++i)
    posteriors[i] = get_posterior(forward[i], backward[i]);
}

static double
forward_algorithm(const vector<double> &ls, const two_by_two &lt,
                  const vector<pair<double, double> > &emit,
                  vector<pair<double, double> > &f) {
  f[0] = make_pair(emit[0].first + ls[0], emit[0].second + ls[1]);
  for (size_t j = 1; j < f.size(); ++j) {
    const size_t i = j - 1;
    f[j].first = emit[j].first + log_sum_log(f[i].first + lt[0][0],
                                             f[i].second + lt[1][0]);
    f[j].second = emit[j].second + log_sum_log(f[i].first + lt[0][1],
                                               f[i].second + lt[1][1]);
  }
  return log_sum_log(f.back().first, f.back().second);
}

static double
backward_algorithm(const vector<double> &ls, const two_by_two &lt,
                   const vector<pair<double, double> > &emit,
                   vector<pair<double, double> > &b) {
  b.back() = make_pair(0.0, 0.0);
  for (size_t j = b.size() - 1; j > 0; --j) {
    const size_t i = j - 1;
    const double bg_a = emit[j].first + b[j].first;
    const double fg_a = emit[j].second + b[j].second;
    b[i].first = log_sum_log(lt[0][0] + bg_a, lt[0][1] + fg_a);
    b[i].second = log_sum_log(lt[1][0] + bg_a, lt[1][1] + fg_a);
  }
  return log_sum_log(b[0].first + emit[0].first + ls[0],
                     b[0].second + emit[0].second + ls[1]);
}

static void
backward_sampling(const two_by_two &lt,
                  const vector<pair<double, double> > &emit,
                  vector<pair<double, double> > &f,
                  vector<bool> &x, std::mt19937 &gen) {

  std::uniform_real_distribution<double> unif(0.0, 1.0);
  
  pair<double, double> b = make_pair(0.0, 0.0);
  double p1 = get_posterior(f.back(), b);
  x.back() = (unif(gen) < p1);

  for (size_t j = f.size() - 1; j > 0; --j) {
    const size_t i = j - 1;
    const double em = x[j] ? emit[j].second : emit[j].first;
    const double bg = f[i].first + lt[0][x[j]] + em;
    const double fg = f[i].second + lt[1][x[j]] + em;
    
    p1 = exp(fg - log_sum_log(fg, bg));
    x[i] = (unif(gen) < p1);
  }
}

static void
summarize_transitions(const vector<pair<double, double> > &f,
                      const vector<pair<double, double> > &b,
                      const double total,
                      const vector<pair<double, double> > &emit,
                      const two_by_two &lt, vector<two_by_two> &joint) {
  
  for (size_t j = 1; j < joint.size(); ++j) {
    const size_t i = j - 1;
    const double left_bg = f[i].first;
    const double left_fg = f[i].second;
    const double right_bg = b[j].first + emit[j].first - total;
    const double right_fg = b[j].second + emit[j].second - total;

    joint[i][0][0] = exp(left_bg + lt[0][0] + right_bg);
    joint[i][0][1] = exp(left_bg + lt[0][1] + right_fg);
    joint[i][1][0] = exp(left_fg + lt[1][0] + right_bg);
    joint[i][1][1] = exp(left_fg + lt[1][1] + right_fg);
  }
}

double
TwoStateHMM::single_iteration(const vector<bool> &obs,
                              vector<pair<double, double> > &log_forward,
                              vector<pair<double, double> > &log_backward,
                              vector<pair<double, double> > &emit,
                              vector<two_by_two> &joint) {
  
  const vector<double> ls = {log(p_fb/(p_bf + p_fb)), log(p_bf/(p_bf + p_fb))};
  const two_by_two lt { {log(1.0 - p_bf), log(p_bf)},
    {log(p_fb), log(1.0 - p_fb)}};
  
  assert(isfinite(ls[0]) && isfinite(ls[1]) && isfinite(lt[0][0]) &&
         isfinite(lt[0][1]) && isfinite(lt[1][0]) && isfinite(lt[1][1]));
  
  get_log_emissions(obs, emit, fg_distr, bg_distr);
  
  const double new_llh = forward_algorithm(ls, lt, emit, log_forward);
  const double backward_llh = backward_algorithm(ls, lt, emit, log_backward);
   
  assert(fabs(get_delta(new_llh, backward_llh)) < tolerance);
  summarize_transitions(log_forward, log_backward, new_llh, emit, lt, joint);
  
  if (get_delta(llh, new_llh) > tolerance) { // not converged
    two_by_two sum_joint = two_by_two(2, vector<double> (2, 0.0));
    for (size_t i = 0; i < joint.size(); ++i) {
      sum_joint[0][0] += joint[i][0][0];
      sum_joint[0][1] += joint[i][0][1];
      sum_joint[1][0] += joint[i][1][0];
      sum_joint[1][1] += joint[i][1][1];
    }
    
    // Update transition probabilities
    const double p_bf_est = sum_joint[0][1] / (sum_joint[0][0] + sum_joint[0][1]);
    const double p_fb_est = sum_joint[1][0] / (sum_joint[1][0] + sum_joint[1][1]);
    assert(p_bf_est > tolerance);
    p_bf = p_bf_est;
    assert(p_fb_est > tolerance);
    p_fb = p_fb_est;
    
    if (!FIX_EMIT) { // Update Emission parameters
      vector<double> posteriors;
      get_posteriors(log_forward, log_backward, posteriors);
      fg_distr.fit(obs, posteriors);
      one_minus(begin(posteriors), end(posteriors), begin(posteriors));
      bg_distr.fit(obs, posteriors);
    }
  }
    return new_llh;
}


double
TwoStateHMM::BaumWelchTraining(const vector<bool> &obs) {
  
  const size_t n_vals = obs.size();
  vector<pair<double, double> > log_forward(n_vals, make_pair(0.0, 0.0));
  vector<pair<double, double> > log_backward(n_vals, make_pair(0.0, 0.0));
  vector<pair<double, double> > emit(n_vals, make_pair(0.0, 0.0));
  vector<two_by_two> joint(n_vals, two_by_two(2, vector<double>(2, 0.0)));
  
  llh = - std::numeric_limits<double>::max();
  double delta = std::numeric_limits<double>::max();
  
  if (VERBOSE) {
    report_param_header_for_verbose();
    report_params_for_verbose(0, fg_distr.p, bg_distr.p, p_fb, p_bf, delta);
  }
  
  for (size_t i = 0; i < max_iterations && (delta > tolerance); ++i) {
    const double new_llh =
    single_iteration(obs, log_forward, log_backward, emit, joint);
    delta = get_delta(llh, new_llh);
    
    if (delta < tolerance) {
      if (VERBOSE)
        cerr << "CONVERGED" << endl;
    }
    else {
      if (VERBOSE)
        report_params_for_verbose(i+1, fg_distr.p, bg_distr.p,
                                  p_fb, p_bf, delta);
      llh = new_llh;
    }
  }
  return llh;
}


void
TwoStateHMM::StatesSampling(const vector<bool> &obs, vector<bool> &x,
                            std::mt19937 &gen) const {
  
  const vector<double> ls = {log(p_fb/(p_bf + p_fb)), log(p_bf/(p_bf + p_fb))};
  const two_by_two lt { {log(1.0 - p_bf), log(p_bf)},
    {log(p_fb), log(1.0 - p_fb)}};
  
  assert(isfinite(ls[0]) && isfinite(ls[1]) && isfinite(lt[0][0]) &&
         isfinite(lt[0][1]) && isfinite(lt[1][0]) && isfinite(lt[1][1]));
  
  const size_t n_vals = obs.size();
  vector<pair<double, double> > log_forward(n_vals, make_pair(0.0, 0.0));
  vector<pair<double, double> > emit(n_vals, make_pair(0.0, 0.0));
  x.resize(n_vals, false);
  
  get_log_emissions(obs, emit, fg_distr, bg_distr);
  
  forward_algorithm(ls, lt, emit, log_forward);
  backward_sampling(lt, emit, log_forward, x, gen);
}


////////////////////////////////////////////////////////////////////////

static void
write_params_file(const string &outfile, const double fg_p, const double bg_p,
                  const double p_fb, const double p_bf) {

  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out.precision(30);
  out << "FG_P\t" << fg_p << "\tBG_P\t" << bg_p
      << "\tF_B\t" << p_fb << "\tB_F\t" << p_bf << endl;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


int
main(int argc, const char **argv) {

  try {

    string outfile;

    const static double tolerance = 1e-10;
    size_t max_iterations = 100;
    size_t rng_seed = std::numeric_limits<size_t>::max();

    // run mode flags
    bool VERBOSE = false;
    bool FIX_EMIT = false;

    string params_in_file;
    string params_out_file;

    double fg_p = 0.8;
    double bg_p = 0.1;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "Sample 2-state HMM path", "<state-file>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("fixemit", 'e', "do not update emit", false, FIX_EMIT);
    opt_parse.add_opt("fgemit", 'F', "foreground emission", false, fg_p);
    opt_parse.add_opt("bgemit", 'B', "background emission", false, bg_p);
    opt_parse.add_opt("params-in", 'P', "HMM parameter file "
                      "(override training)", false, params_in_file);
    opt_parse.add_opt("params-out", 'p', "write HMM parameters to this "
                      "file (default: none)", false, params_out_file);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.set_show_defaults();
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
    const string states_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    // READ OBSERVED DATA
    if (VERBOSE)
      cerr << "[OBTAINING OBSERVED SEQUENCE]" << endl;
    std::ifstream states_in(states_file);
    vector<bool> obs;
    copy(istream_iterator<bool>(states_in), istream_iterator<bool>(),
         std::back_inserter(obs));
    
    // HMM INITIALIZATION
    if (VERBOSE)
      cerr << "[HMM INITIALIZATION]" << endl;
    TwoStateHMM hmm(tolerance, max_iterations, VERBOSE, FIX_EMIT);

    if (!params_in_file.empty()) { // load parameters file
      hmm.initialize(params_in_file);
      max_iterations = 0;
    } else {
      hmm.initialize(obs);
      hmm.fg_distr.p = fg_p;
      hmm.bg_distr.p = bg_p;
    }
    
    // HMM TRAINING
    if (VERBOSE)
      cerr << "[HMM TRAINING]" << endl;
    
    if (max_iterations > 0)
      hmm.BaumWelchTraining(obs);

    // HMM OUTPUT
    if (!params_out_file.empty())
      write_params_file(params_out_file, hmm.fg_distr.p, hmm.bg_distr.p,
                        hmm.p_fb, hmm.p_bf);
    
    // HMM SAMPLE HIDDEN STATES
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    std::mt19937 gen(rng_seed);

    if (VERBOSE)
      cerr << "[HMM SAMPLING]" << endl;
    
    vector<bool> states;
    hmm.StatesSampling(obs, states, gen);
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    copy(begin(states), end(states), std::ostream_iterator<double>(out, "\n"));
  

  }
  catch (runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
