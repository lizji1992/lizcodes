/*
  Copyright (C) 2011 University of Southern California
  Authors: Andrew D. Smith, Song Qiang

  This file is part of rmap.

  rmap is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  rmap is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with rmap; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "TwoStateCTHMM.hpp"

#include <iomanip>
#include <numeric>
#include <limits>
#include <cmath>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

#include "numerical-utils.hpp"

// #define DEBUG

using std::vector;
using std::pair;
using std::setw;
using std::max;
using std::min;
using std::cerr;
using std::endl;
using std::string;
using std::setprecision;

TwoStateCTHMM::TwoStateCTHMM(
    const std::vector<std::pair<double, double> > &_observations,
    const std::vector<size_t> &_time,
    const std::vector<size_t> &_reset_points,
    const size_t _MAX_LEN, const double mp, const double tol,
    const size_t max_itr, const bool v) : 
    observations(_observations), time(_time), reset_points(_reset_points),
    meth_lp(_observations.size()),
    unmeth_lp(_observations.size()),
    fg_log_likelihood(_observations.size()),
    bg_log_likelihood(_observations.size()),
    forward(_observations.size()),
    backward(_observations.size()),
    fg_posteriors(_observations.size()),
    bg_posteriors(_observations.size()),
    MAX_LEN(_MAX_LEN), MIN_PROB(mp), tolerance(tol),
    max_iterations(max_itr), VERBOSE(v) 
{
    for (size_t i = 0; i < observations.size(); ++i)
    {
        const double m = observations[i].first;
        const double u = observations[i].second;
        
        meth_lp[i] = 
            log(std::min(std::max(m/(m + u), 1e-2), 1.0 - 1e-2));
        unmeth_lp[i] = 
            log(std::min(std::max(u/(m + u), 1e-2), 1.0 - 1e-2));
    }
}

void
TwoStateCTHMM::set_parameters(
    const betabin &_fg_emission,
    const betabin &_bg_emission,
    const Distro &_fg_duration,
    const Distro &_bg_duration)
{
    fg_emission = _fg_emission;
    bg_emission = _bg_emission;
    fg_duration = _fg_duration;
    bg_duration = _bg_duration;
    update_observation_likelihood();

    lp_sf = log(0.5);
    lp_sb = log(0.5);
    lp_ft = log(1e-10);
    lp_bt = log(1e-10);
}

void
TwoStateCTHMM::get_parameters(
    betabin &_fg_emission,
    betabin &_bg_emission,
    Distro &_fg_duration,
    Distro &_bg_duration)
{
    _fg_emission = fg_emission;
    _bg_emission = bg_emission;
    _fg_duration = fg_duration;
    _bg_duration = bg_duration;
}


//////////////////////////////////////////////
////// forward and backward algorithms  //////
//////////////////////////////////////////////
void
TwoStateCTHMM::update_observation_likelihood()
{
    fg_log_likelihood.front() = fg_emission(observations.front());
    bg_log_likelihood.front() = bg_emission(observations.front());
    
    for (size_t i = 1; i < observations.size(); ++i)
    {
        fg_log_likelihood[i] =
            fg_log_likelihood[i - 1] + fg_emission(observations[i]);
        bg_log_likelihood[i] =
            bg_log_likelihood[i - 1] + bg_emission(observations[i]);
    }
}

double
TwoStateCTHMM::fg_segment_log_likelihood(
    const size_t start, const size_t end)
{
    return
        (start == 0)
        ? fg_log_likelihood[end - 1]
        : fg_log_likelihood[end - 1] - fg_log_likelihood[start - 1];
}

double
TwoStateCTHMM::bg_segment_log_likelihood(
    const size_t start, const size_t end)
{
    return
        (start == 0)
        ? bg_log_likelihood[end - 1]
        : bg_log_likelihood[end - 1] - bg_log_likelihood[start - 1];
}


double
TwoStateCTHMM::forward_algorithm(const size_t start, const size_t end) 
{
#ifdef DEBUG
    cerr << "check enter forward_algorithm: "<< "OK" << endl;
#endif

    const double fg_rate = 1.0 / fg_duration.get_params().front();
    const double bg_rate = 1.0 / bg_duration.get_params().front();

    std::fill(forward.begin() + start, forward.begin() + end,
              std::make_pair(0.0, 0.0)); 

    forward[start].first =
        lp_sf + fg_segment_log_likelihood(start, start + 1);
    forward[start].second =
        lp_sb + bg_segment_log_likelihood(start, start + 1);

    for (size_t i = start + 1; i < end; ++i)
    {
        const double ff =
            (bg_rate +
             fg_rate * exp(-(fg_rate + bg_rate) * (time[i] - time[i - 1])))
            / (fg_rate + bg_rate);
        const double bb =
            (fg_rate
             + bg_rate * exp(-(fg_rate + bg_rate) * (time[i] - time[i - 1])))
            / (fg_rate + bg_rate);

        const double lp_ff = log(ff);
        const double lp_fb = log(1 - ff);
        const double lp_bf = log(1 - bb);
        const double lp_bb = log(bb);
        

        //  in foreground segment
        forward[i].first =
            log_sum_log(forward[i - 1].first + lp_ff,
                        forward[i - 1].second + lp_bf)
            + fg_segment_log_likelihood(i, i + 1);
        

        //  in background segment
        forward[i].second =
            log_sum_log(forward[i - 1].first + lp_fb,
                        forward[i - 1].second + lp_bb)
            + bg_segment_log_likelihood(i, i + 1);
    }

#ifdef DEBUG
    cerr << "check forward_algorithm: "<< "OK" << endl;
#endif

    return log_sum_log(forward[end - 1].first + lp_ft,
                       forward[end - 1].second + lp_bt);
}

double
TwoStateCTHMM::backward_algorithm(const size_t start, const size_t end)
{
#ifdef DEBUG
    cerr << "check backward_algorithm: "<< "OK" << endl;
#endif

    const double fg_rate = 1.0 / fg_duration.get_params().front();
    const double bg_rate = 1.0 / bg_duration.get_params().front();

    std::fill(backward.begin() + start, backward.begin() + end,
              std::make_pair(0.0, 0.0)); 

    backward[end - 1].first = lp_ft;
    backward[end - 1].second = lp_bt;

    for (size_t i = end - 2; i >= start; --i)
    {
        const double ff =
            (bg_rate +
             fg_rate * exp(-(fg_rate + bg_rate) * (time[i + 1] - time[i])))
            / (fg_rate + bg_rate);
        const double bb =
            (fg_rate
             + bg_rate * exp(-(fg_rate + bg_rate) * (time[i + 1] - time[i])))
            / (fg_rate + bg_rate);

        const double lp_ff = log(ff);
        const double lp_fb = log(1 - ff);
        const double lp_bf = log(1 - bb);
        const double lp_bb = log(bb);

        // observesvation i is foreground segment end
        backward[i].first =
            log_sum_log(lp_ff + fg_segment_log_likelihood(i + 1, i + 2)
                        + backward[i + 1].first,
                        lp_fb + bg_segment_log_likelihood(i + 1, i + 2)
                        + backward[i + 1].second);

        backward[i].second =
            log_sum_log(lp_bf + fg_segment_log_likelihood(i + 1, i + 2)
                        + backward[i + 1].first,
                        lp_bb + bg_segment_log_likelihood(i + 1, i + 2)
                        + backward[i + 1].second);
 
        if (i == 0) break;
   }

    const double llh =
        log_sum_log(lp_sf + fg_segment_log_likelihood(start, start + 1)
                    + backward[start].first,
                    lp_sb + bg_segment_log_likelihood(start, start + 1)
                    + backward[start].second);
    
#ifdef DEBUG
    cerr << "check backward_algorithm: "<< "OK" << endl;
#endif
    return llh;
}

//////////////////////////////////////////////
//////       Baum-Welch Training        //////
//////////////////////////////////////////////
// Expectation
void
TwoStateCTHMM::estimate_state_posterior(const size_t start, const size_t end)  
{
#ifdef DEBUG
    cerr << "check enter estimate_state_posterior: "<< "OK" << endl;
#endif

    for (size_t i = start; i < end; ++i)
    {
        const double fg_evidence = forward[i].first + backward[i].first;
        const double bg_evidence = forward[i].second + backward[i].second;

        const double denom = log_sum_log(fg_evidence, bg_evidence);
        fg_posteriors[i] = exp(fg_evidence - denom);
        bg_posteriors[i] = exp(bg_evidence - denom);

        assert(fabs(fg_posteriors[i] + bg_posteriors[i] - 1.0) < 1e-6);
    }

#ifdef DEBUG
    cerr << "check estimate_state_posterior: "<< "OK" << endl;
#endif
}

// Maximization
void 
TwoStateCTHMM::estimate_parameters()
{
// /////
//     cerr << "check enter estimate_parameters: "<< "OK" << endl;
// /////
    static const bool FG_STATE = true;
//    static const bool BG_STATE = false;
    
    fg_emission.fit(meth_lp, unmeth_lp, fg_posteriors);
    bg_emission.fit(meth_lp, unmeth_lp, bg_posteriors);
    update_observation_likelihood();
    
    vector<double> fg_lengths, bg_lengths;
    for (size_t idx = 0; idx < reset_points.size() - 1; ++idx)
    {
        const size_t start = reset_points[idx];
        const size_t end = reset_points[idx + 1];
        bool prev_state = fg_posteriors[start] > bg_posteriors[start];
        size_t domain_start = start;
        for (size_t i = start + 1; i < end; ++i)
        {
            const bool state = fg_posteriors[i] > bg_posteriors[i];
            if (state != prev_state)
            {
                const size_t domain_end = i - 1;
                if (prev_state == FG_STATE)
                    fg_lengths.push_back(time[domain_end] - time[domain_start]);
                else
                    bg_lengths.push_back(time[domain_end] - time[domain_start]);

                prev_state = state;
                domain_start = i;
            }
        }
        const size_t domain_end = end - 1;
        if (prev_state == FG_STATE)
            fg_lengths.push_back(time[domain_end] - time[domain_start]);
        else
            bg_lengths.push_back(time[domain_end] - time[domain_start]);
    }

    if (fg_lengths.size() > 0)
        fg_duration.estimate_params_ml(fg_lengths);

    if (bg_lengths.size() > 0)
        bg_duration.estimate_params_ml(bg_lengths);

// /////
//     cerr << "check estimate_parameters: "<< "OK" << endl;
// /////
}

double
TwoStateCTHMM::single_iteration()
{
// /////
//     cerr << "check enter single_iteration: "<< "OK" << endl;
// /////

    double total_score = 0;

    for (size_t i = 0; i < reset_points.size() - 1; ++i)
    {
        const double forward_score =
            forward_algorithm(reset_points[i], reset_points[i + 1]);
        const double backward_score =
            backward_algorithm(reset_points[i], reset_points[i + 1]);
        
        assert(fabs((forward_score - backward_score)
                    / max(forward_score, backward_score))
               < 1e-10);
        estimate_state_posterior(reset_points[i], reset_points[i + 1]);
        total_score += forward_score;
    }
    
    estimate_parameters();
    
// /////
//     cerr << "check single_iteration: "<< "OK" << endl;
// /////
    
    return total_score;
}

double
TwoStateCTHMM::BaumWelchTraining()  
{
#ifdef DEBUG
    cerr << "check enter BaumWelchTraining: " << "OK" << endl;
#endif

    if (VERBOSE)
        cerr << setw(5)  << "ITR"
             << setw(16) << "FG Emission"
             << setw(18) << "FG Duration"
             << setw(16) << "BG Emission"
             << setw(16) << "BG Duration"
             << setw(14) << "Likelihood"
             << setw(14) << "DELTA"
             << endl;
  
    double prev_total = -std::numeric_limits<double>::max();
  
    for (size_t i = 0; i < max_iterations; ++i) 
    {
        const betabin old_fg_emission = fg_emission;
        const betabin old_bg_emission = bg_emission;
        const Distro old_fg_duration = fg_duration;
        const Distro old_bg_duration = bg_duration;
        
        double total = single_iteration();
        
        if (VERBOSE)
            cerr << setw(5) << i + 1
                 << setw(16) << old_fg_emission.tostring()
                 << setw(18) << old_fg_duration.tostring()
                 << setw(16) << old_bg_emission.tostring()
                 << setw(16) << old_bg_duration.tostring()
                 << setw(14) << total
                 << setw(14) << (total - prev_total)/std::fabs(total)
                 << endl;
        
        if ((total - prev_total)/std::fabs(total) < tolerance)
        {
            fg_emission = old_fg_emission;
            bg_emission = old_bg_emission;
            update_observation_likelihood();
            
            fg_duration = old_fg_duration;
            bg_duration = old_bg_duration;
            
            if (VERBOSE)
                cerr << "CONVERGED" << endl << endl;
            break;
        }
        prev_total = total;
    }

#ifdef DEBUG
     cerr << "check exit BaumWelchTraining: "<< "OK" << endl;
#endif
    return prev_total;
}

//////////////////////////////////////////////
//////          Posterior decoding      //////
//////////////////////////////////////////////
double
TwoStateCTHMM::PosteriorDecoding()
{
// /////
//     cerr << "check enter PosteriorDecoding: "<< "OK" << endl;
// /////

    double total_score = 0;

    for (size_t i = 0; i < reset_points.size() - 1; ++i)
    {
        const double forward_score =
            forward_algorithm(reset_points[i], reset_points[i + 1]);
        const double backward_score =
            backward_algorithm(reset_points[i], reset_points[i + 1]);
        
        assert(fabs((forward_score - backward_score)
                    / max(forward_score, backward_score))
                    < 1e-10);
        estimate_state_posterior(reset_points[i], reset_points[i + 1]);
        total_score += forward_score;
    }

// /////
//     cerr << "check exit PosteriorDecoding: "<< "OK" << endl;
// /////

    return total_score;
}



//////////////////////////////////////////////
//////          export result           //////
//////////////////////////////////////////////

void
TwoStateCTHMM::get_posterior_scores(
    std::vector<double> & scores,
    std::vector<bool> & classes)
{
    scores = fg_posteriors;

    classes.resize(observations.size());
    for (size_t i = 0; i < observations.size(); ++i)
        classes[i] = fg_posteriors[i] > bg_posteriors[i];
}



