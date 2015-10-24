//
//  distribution.h
//  
//
//  Created by Liz Ji on 6/24/15.
//
//

#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

#include <utility>
#include <string>
#include <vector>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

using std::vector;
using std::string;

typedef vector< vector<double> > matrix;

// struct BetaBin;
struct BetaBin
{
  BetaBin() : alpha(1), beta(1), lnbeta_helper(gsl_sf_lnbeta(1, 1)) {}
  BetaBin(const double a, const double b) :
  alpha(a), beta(b), lnbeta_helper(gsl_sf_lnbeta(a, b)) {}
  
  double operator()(const std::pair<double, double> &val) const;
  void fit(const std::vector<double> &vals_a,
           const std::vector<double> &vals_b,
           const std::vector<double> &p);
  string tostring() const;
  
  double alpha;
  double beta;
  double lnbeta_helper;
  static const double tolerance = 1e-10;
};


struct NegBin
{
  NegBin() : r(2), p(0.5) {set_helpers();}
  NegBin(const size_t a, const double b) : r(a), p(b) {set_helpers();}
  
  void set_helpers();
  
  //static const double max_allowed_alpha = 100;
  //static const double min_allowed_alpha = 1e-20;
  //static const double alpha_allowed_error = 1e-10;
  
  size_t r; // num of failure
  double p; // prob of failure
  double r_log_p_minus_lngamma_r; // helper
};


//////////////////////////////////////////////
//////       struct CTHMM duration      //////
//////////////////////////////////////////////


class ExpTransEstimator {
public:
  
  ExpTransEstimator() : a(0.02), b(0.002),
                        step_size(0.01), tolerance(1e-10) {}
  ExpTransEstimator(const double _a,
                    const double _b) : a(_a), b(_b),
                                       step_size(0.01), tolerance(1e-10) {}
  ExpTransEstimator(const double _a, const double _b, const double s,
                    const double t) : a(_a), b(_b),
                                      step_size(s), tolerance(t) {}
  
  void set_stepsize(const double s) {step_size = s;}
  void set_tolerance(const double t) {tolerance = t;}
  
  double get_a() const {return a;}
  double get_b() const {return b;}
  
  
  // GRADIENT ASCENT METHOD
  void mle_GradAscent(const matrix &r, const vector<size_t> &t);
  void GA_stepforward(const double grad_a, const double grad_b,
                      const double &old_llh, double &new_llh,
                      const matrix &r, const vector<size_t> &t);

private:
  
  void calc_internal_data(const vector<size_t> &t, const double u,
                          const double v, vector<double> &log_a,
                          vector<double> &log_1_a,
                          vector<double> &neg_bt) const;
  double calc_log_llh(const matrix &r, const vector<double> &log_a,
                      const vector<double> &log_1_a,
                      const vector<double> &neg_bt) const;
  
  double llh_grad_a(const matrix &lr, const vector<double> &neg_bt) const;
  double llh_grad_b(const matrix &lr, const vector<double> &neg_bt) const;

  
  //  parameters
  double a, b;
  double step_size;
  double tolerance;
  
};

//////////////////////////////////////////////
//////       numerical                  //////
//////////////////////////////////////////////

inline double
log_sum_log(const double p, const double q);

inline double
log_sub_log(const double p, const double q);

/*
class NegBinomDistro : public Distro_ {
public:
  
  NegBinomDistro() : Distro_(std::vector<double>(2, 0)) {}
  NegBinomDistro(std::vector<double> p) : Distro_(p) {set_helpers();}
  NegBinomDistro(const NegBinomDistro &rhs);
  NegBinomDistro& operator=(const NegBinomDistro &rhs);
  ~NegBinomDistro() {}
  void set_helpers();
  double sample() const;
  size_t required_params() const {return 2;}
  void set_params(const std::vector<double> &p);
  void estimate_params_ml(const std::vector<double> &vals);
  void estimate_params_ml(const std::vector<double> &vals,
                          const std::vector<double> &probs);
  void estimate_params_ml(const std::vector<double> &vals,
                          const std::vector<double> &scales,
                          const std::vector<double> &probs);
  
  double log_likelihood(double val) const;
  double log_likelihood(const double &val, const double &scale) const;

};
*/
#endif
