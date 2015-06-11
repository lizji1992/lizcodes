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

#include "vdhmm-utils.hpp"
#include "TwoStateHDHMM-linear.hpp"

#include <iomanip>
#include <numeric>
#include <limits>
#include <cmath>
#include <queue>
#include <utility>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

using std::vector;
using std::pair;
using std::setw;
using std::max;
using std::min;
using std::cerr;
using std::endl;
using std::string;
using std::setprecision;
using std::priority_queue;

TwoStateHDHMMLinear::TwoStateHDHMMLinear(
    const std::vector<std::pair<double, double> > &_observations,
    const std::vector<size_t> &_reset_points,
    const double mp, const double tol,
    const size_t max_itr, const bool v,
    const size_t _MAX_LEN) : 
    observations(_observations), reset_points(_reset_points),
    meth_lp(_observations.size()), unmeth_lp(_observations.size()),
    fg_log_likelihood(_observations.size()),
    bg_log_likelihood(_observations.size()),
    forward(_observations.size()), backward(_observations.size()), 
    fg_posteriors(_observations.size()), bg_posteriors(_observations.size()), 
    forward_support_points(_observations.size()),
    backward_support_points(_observations.size()),
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
TwoStateHDHMMLinear::init(
    const vector<pair<double, double> > & _observations,
    const vector<size_t> & _reset_points)
{
    observations = _observations;
    reset_points = _reset_points;

    meth_lp.resize(observations.size());
    unmeth_lp.resize(observations.size());

    for (size_t i = 0; i < observations.size(); ++i)
    {
        const double m = observations[i].first;
        const double u = observations[i].second;
        
        meth_lp[i] = 
            log(std::min(std::max(m/(m + u), 1e-2), 1.0 - 1e-2));
        unmeth_lp[i] = 
            log(std::min(std::max(u/(m + u), 1e-2), 1.0 - 1e-2));
    }


    fg_posteriors.resize(observations.size());
    bg_posteriors.resize(observations.size());

    forward.resize(observations.size());
    backward.resize(observations.size());

    forward_support_points.resize(observations.size());
    backward_support_points.resize(observations.size());

/////
    cerr << "FUNCTION TO BE DEPRECATED" << endl;
/////
}

void
TwoStateHDHMMLinear::initialize_parameters(
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


//////////////////////////////////////////////
////// forward and backward algorithms  //////
//////////////////////////////////////////////
void
TwoStateHDHMMLinear::update_observation_likelihood()
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
TwoStateHDHMMLinear::fg_segment_log_likelihood(
    const size_t start, const size_t end)
{
    return
        (start == 0)
        ? fg_log_likelihood[end - 1]
        : fg_log_likelihood[end - 1] - fg_log_likelihood[start - 1];
}

double
TwoStateHDHMMLinear::bg_segment_log_likelihood(
    const size_t start, const size_t end)
{
    return
        (start == 0)
        ? bg_log_likelihood[end - 1]
        : bg_log_likelihood[end - 1] - bg_log_likelihood[start - 1];
}


double
TwoStateHDHMMLinear::forward_algorithm(const size_t start, const size_t end) 
{
/////
     cerr << "check enter forward_algorithm: "<< "OK" << endl;
/////

    const double self_lp = log(1 - bg_duration.get_params().front());
    const double switch_lp = log(bg_duration.get_params().front());

    std::fill(forward.begin() + start, forward.begin() + end,
              std::make_pair(0.0, 0.0)); 
    
    std::fill(forward_support_points.begin() + start,
              forward_support_points.begin() + end,
              vector<size_t>());

    priority_queue<pair<double, size_t> > candidates, new_candidates;
    
    forward[start].first =
        lp_sf + fg_segment_log_likelihood(start, start + 1)
        + fg_duration.log_likelihood(1);
    
    forward[start].second =
        lp_sb + bg_segment_log_likelihood(start, start + 1);
    
    for (size_t i = start + 1; i < end; ++i)
    {
        // Index i terminates a foreground domain
        const size_t ending = i + 1;       // exclusive

        new_candidates = priority_queue<pair<double, size_t> >();
        candidates.push(
            std::make_pair(std::numeric_limits<double>::max(), i));
        const size_t num_candidates = candidates.size();

        for (size_t j = 0; j < num_candidates; ++j)
        {
            const size_t beginning = candidates.top().second;
            candidates.pop();

            // foreground segment
            const double fg_seg_llh = (beginning == start) ?      // segment [start, end)
                lp_sf
                + fg_segment_log_likelihood(beginning, ending)
                + fg_duration.log_likelihood(ending - beginning)
                :
                forward[beginning - 1].second + switch_lp
                + fg_segment_log_likelihood(beginning, ending)
                + fg_duration.log_likelihood(ending - beginning);

            forward[i].first = log_sum_log(forward[i].first, fg_seg_llh);

            forward_support_points[i].push_back(beginning);

            const double fg_seg_llh_adjusted = 
                fg_seg_llh
                - fg_duration.log_likelihood(ending - beginning)
                + lp_segment_longer_than(ending - beginning,
                                         fg_duration);
            new_candidates.push(std::make_pair(fg_seg_llh_adjusted, beginning));
            if (fg_seg_llh - forward[i].first < log(tolerance))
                break;
        }

// /////
//         cerr << "support points: [" << i << "] ";
//         for (size_t j = 0; j < forward_support_points[i].size(); ++j)
//             cerr << forward_support_points[i][j] << "\t";
//         cerr << endl;
// /////

        std::sort(forward_support_points[i].begin(),
                  forward_support_points[i].end(),
                  std::greater<size_t>());
        candidates = new_candidates;



        //  in background segment
        forward[i].second =
            log_sum_log(forward[i - 1].first,
                        forward[i - 1].second + self_lp)
            + bg_segment_log_likelihood(i, i + 1);
    }

    return log_sum_log(forward[end - 1].first + lp_ft,
                       forward[end - 1].second + lp_bt);

//     /////
     cerr << "check forward_algorithm: "<< "OK" << endl;
// /////
}

double
TwoStateHDHMMLinear::backward_algorithm(const size_t start, const size_t end)
{
// /////
     cerr << "check backward_algorithm: "<< "OK" << endl;
// /////

    const double self_lp = log(1 - bg_duration.get_params().front());
    const double switch_lp = log(bg_duration.get_params().front());

    std::fill(backward.begin() + start, backward.begin() + end,
              std::make_pair(0.0, 0.0)); 

    std::fill(backward_support_points.begin() + start,
              backward_support_points.begin() + end,
              vector<size_t>());

    priority_queue<pair<double, size_t> > candidates, new_candidates;

    backward[end - 1].first = lp_ft;
    backward[end - 1].second = lp_bt;

    for (size_t i = end - 2; i >= start; --i)
    {
        // observesvation i is foreground segment end
        backward[i].first =
            bg_segment_log_likelihood(i + 1, i + 2)
            + backward[i + 1].second;
        
        // observation i in a background segment
        // remain in background
        backward[i].second =
            self_lp
            + bg_segment_log_likelihood(i + 1, i + 2)
            + backward[i + 1].second; 
        
        // switch to a foreground segment
        const size_t beginning = i + 1;
        new_candidates = priority_queue<pair<double, size_t> >();
        candidates.push(
            std::make_pair(std::numeric_limits<double>::max(), beginning + 1));
        const size_t num_candidates = candidates.size();

        for (size_t j = 0; j < num_candidates; ++j)
        {
            const size_t ending = candidates.top().second;
            candidates.pop();

            // foreground segment
            const double fg_seg_llh =       // segment [start, end)
                switch_lp  
                + fg_segment_log_likelihood(beginning, ending)
                + fg_duration.log_likelihood(ending - beginning)
                + backward[ending - 1].first;
            backward[i].second = log_sum_log(backward[i].second, fg_seg_llh);

            backward_support_points[beginning].push_back(ending);

            const double fg_seg_llh_adjusted = 
                fg_seg_llh
                - fg_duration.log_likelihood(ending - beginning)
                + lp_segment_longer_than(ending - beginning,
                                         fg_duration);
            new_candidates.push(std::make_pair(fg_seg_llh_adjusted, ending));
            if (fg_seg_llh - backward[i].second < log(tolerance))
                break;
        }
        
// /////
//         cerr << "support points: [" << i << "] ";
//         for (size_t j = 0; j < backward_support_points[beginning].size(); ++j)
//             cerr << backward_support_points[beginning][j] << "\t";
//         cerr << endl;
// /////

        std::sort(backward_support_points[beginning].begin(),
                  backward_support_points[beginning].end(),
                  std::greater<size_t>());
        candidates = new_candidates;

        if (i == 0) break;
   }

    // whole likelihood
    double llh = 0;
    
    // the first segment is foreground
    const size_t beginning = start;
    candidates.push(
        std::make_pair(std::numeric_limits<double>::max(), beginning + 1));
    const size_t num_candidates = candidates.size();

    for (size_t j = 0; j < num_candidates; ++j)
    {
        const size_t ending = candidates.top().second;
        candidates.pop();

        // foreground segment
        const double fg_seg_llh =       // segment [start, end)
            lp_sf  
            + fg_segment_log_likelihood(beginning, ending)
            + fg_duration.log_likelihood(ending - beginning)
            + backward[ending - 1].first;
        llh = log_sum_log(llh, fg_seg_llh);

        backward_support_points[beginning].push_back(ending);
        if (fg_seg_llh - llh < log(tolerance))
            break;
    }


    // the first segment is background
    llh = log_sum_log(llh, 
                      lp_sb
                      + bg_segment_log_likelihood(start, start + 1)
                      + backward[start].second);

    return llh;

// /////
     cerr << "check backward_algorithm: "<< "OK" << endl;
// /////

}

//////////////////////////////////////////////
//////       Baum-Welch Training        //////
//////////////////////////////////////////////
// Expectation
void
TwoStateHDHMMLinear::estimate_state_posterior(const size_t start, const size_t end)  
{
// /////
//     cerr << "check enter estimate_state_posterior: "<< "OK" << endl;
// /////

    // const double self_lp = log(1 - bg_duration.get_params().front());
    const double switch_lp = log(bg_duration.get_params().front());

    vector<double> fg_evidence(end - start, 0), bg_evidence(end - start, 0);
    for (size_t s = start; s < end; ++s)
    {
        // foreground
        double accu_evidence = 0;
        const size_t num_supp_points = 
            backward_support_points[s].size();

        for (size_t i = 0; i < num_supp_points; ++i)
        {
            const size_t e = backward_support_points[s][i];

            const double evidence = (s == start) ?
                lp_sf
                + fg_duration.log_likelihood(e - s)
                + fg_segment_log_likelihood(s, e)
                + backward[e - 1].first
                :
                forward[s - 1].second + switch_lp
                + fg_duration.log_likelihood(e - s)
                + fg_segment_log_likelihood(s, e)
                + backward[e - 1].first;
            accu_evidence = log_sum_log(accu_evidence, evidence);

            // update observation state posterior
            const size_t first_left_cp =
                (i == num_supp_points - 1) ?
                s : backward_support_points[s][i + 1];
            for (size_t j = first_left_cp; j < e; ++j)
                fg_evidence[j] =
                    log_sum_log(fg_evidence[j], accu_evidence);
        }
        
        // background
        bg_evidence[s - start] = forward[s].second + backward[s].second;
    }

    // state posterior
    for (size_t i = start; i < end; ++i)
    {
        const double denom = log_sum_log(fg_evidence[i - start],
                                         bg_evidence[i - start]);
        fg_posteriors[i] = exp(fg_evidence[i - start] - denom);
        bg_posteriors[i] = exp(bg_evidence[i - start] - denom);

        assert(fabs(fg_posteriors[i] + bg_posteriors[i] - 1.0) < 1e-6);
    }

// /////
//     cerr << "check estimate_state_posterior: "<< "OK" << endl;
// /////

}

// Maximization
void 
TwoStateHDHMMLinear::estimate_parameters()
{
// /////
//     cerr << "check enter estimate_parameters: "<< "OK" << endl;
// /////
    static const bool FG_STATE = true;
//    static const bool BG_STATE = false;
    
    fg_emission.fit(meth_lp, unmeth_lp, fg_posteriors);
    bg_emission.fit(meth_lp, unmeth_lp, bg_posteriors);

// /////
//     cerr << "check fg_emission: "<< fg_emission.tostring() << endl;
//     cerr << "check bg_emission: "<< bg_emission.tostring() << endl;
// /////

    
    vector<double> fg_lengths, bg_lengths;
    for (size_t idx = 0; idx < reset_points.size() - 1; ++idx)
    {
        const size_t start = reset_points[idx];
        const size_t end = reset_points[idx + 1];
        bool prev_state = fg_posteriors[start] > bg_posteriors[start];
        size_t len = 1;
        for (size_t i = start + 1; i < end; ++i)
        {
            const bool state = fg_posteriors[i] > bg_posteriors[i];
            if (state == prev_state)
                ++len;
            else
            {
                if (prev_state == FG_STATE)
                    fg_lengths.push_back(len);
                else
                    bg_lengths.push_back(len);

                prev_state = state;
                len = 1;
            }
        }
    }

    if (fg_lengths.size() > 0)
        fg_duration.estimate_params_ml(fg_lengths);

    if (bg_lengths.size() > 0)
        bg_duration.estimate_params_ml(bg_lengths);

    update_observation_likelihood();
// /////
//     cerr << "check estimate_parameters: "<< "OK" << endl;
// /////

}

double
TwoStateHDHMMLinear::single_iteration()
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
        
        if (fabs((forward_score - backward_score)
                 / max(forward_score, backward_score))
            > 1e-10)
            /////
            cerr << "check single_iteration: "<< forward_score << "\t"
                 << backward_score << endl;
/////

        // assert(fabs((forward_score - backward_score)
        //             / max(forward_score, backward_score))
        //             < 1e-10);
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
TwoStateHDHMMLinear::BaumWelchTraining()  
{
// /////
     cerr << "check enter BaumWelchTraining: "<< "OK" << endl;
// /////

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
            fg_duration = old_fg_duration;
            bg_duration = old_bg_duration;

            if (VERBOSE)
                cerr << "CONVERGED" << endl << endl;
            break;
        }
        prev_total = total;
    }

// /////
     cerr << "check BaumWelchTraining: "<< "OK" << endl;
// /////

    return prev_total;
}


//////////////////////////////////////////////
//////          export result           //////
//////////////////////////////////////////////

void
TwoStateHDHMMLinear::get_posterior_scores(
    std::vector<double> & scores,
    std::vector<bool> & classes)
{
    scores = fg_posteriors;

    classes.resize(observations.size());
    for (size_t i = 0; i < observations.size(); ++i)
        classes[i] = fg_posteriors[i] > bg_posteriors[i];

// /////
//     cerr << "check get_posterior_scores: "<< "OK" << endl;
// /////

}




