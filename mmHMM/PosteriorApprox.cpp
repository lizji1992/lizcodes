#include <string>
#include <vector>
#include <iostream>
#include <iomanip> // std::setw
#include <numeric> //accumulate
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTree.hpp"
#include "EpiPhyloHMM.hpp"
#include "BetaBin.hpp"


using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::setw;
using std::pair;
using std::make_pair;
using std::setprecision;
using std::max;
using std::min;


inline double
log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}


struct val_ind{
  double val;
  size_t ind;
};

struct by_val { 
    bool operator()(val_ind const &left, val_ind const &right) { 
        return left.val < right.val;
    }
};

//////////// Defined in HistMethEvo.cpp/////////////
// void
// rate_to_prob(const vector<double> &Q, const double branch,
//              const bool anc, vector<double> &p);
// //only inheritance
// double
// get_single_firstprob(const bool root, const bool anc, 
//                      const bool cur, const vector<double> &G, 
//                      const vector<double> &Q, const double branch);
//
// //combined processes
// double
// get_single_transprob(const bool root, const bool anc, 
//                      const bool prev, const bool cur,
//                      const vector<double> &G, const vector<double> &Q,
//                      const double branch);
///////////////////////////////////////////////////

double
single_update(const size_t nodeidx, size_t pos, 
              const vector<vector<size_t> > &childidx,
              const vector<double> &G, const vector<double> &Q, 
              const vector<double> &branches,
              const vector<vector<double> > &node_posterior){
  vector<double> pi(2,0);
  pi[0] = (1-G[1])/(2-G[0]-G[1]);
  pi[1] = (1-G[0])/(2-G[0]-G[1]);
  
  size_t lid = childidx[nodeidx][0];
  size_t rid = childidx[nodeidx][1];
  vector<vector<double> > joint_children(2, vector<double>(2,0));
  joint_children[0][0] = node_posterior[pos][lid]*node_posterior[pos][rid];
  joint_children[0][1] = node_posterior[pos][lid]*(1-node_posterior[pos][rid]);
  joint_children[1][0] = (1-node_posterior[pos][lid])*node_posterior[pos][rid];
  joint_children[1][1] = (1.0-node_posterior[pos][lid])*(1- node_posterior[pos][rid]); 
  vector<double> p0,p1;
  rate_to_prob(Q, branches[nodeidx], true, p0); 
  rate_to_prob(Q, branches[nodeidx], false, p1); 
  vector<double> tmp(2, 0.0);
  for(size_t k = 0; k <1; ++k){
    for(size_t l = 0; l <1; ++l){
      tmp[0] += pi[0]*p0[k]*p0[l]/(pi[0]*p0[k]*p0[l] + pi[1]*p1[k]*p1[l] )*joint_children[k][l];
      tmp[1] += pi[1]*p1[k]*p1[l]/(pi[0]*p0[k]*p0[l] + pi[1]*p1[k]*p1[l] )*joint_children[k][l];
    }
  }
  assert(tmp[0]+tmp[1] > 0);
  return ( tmp[0]/(tmp[0]+tmp[1]) );
}

double 
iterate_update(const bool root, const bool leaf, 
               const bool start, const bool end, 
               const size_t pos, const size_t nodeid, 
               const vector<vector<size_t> > &childidx,
               const vector<size_t> &parentidx,
               const vector<size_t> &leafidx,
               const vector<vector<size_t> > &clade_leaf_idx,
               const vector<double> &G, const vector<double> &Q, 
               const vector<double> &branches,
               const vector<pair<double, double> > &Fparams,
               const vector<pair<double, double> > &Bparams,
               const HistoryME &history,
               const vector<vector<pair<double, double> > > meth,
               vector<vector<double> > &node_posterior){
  double oldscore = node_posterior[pos][nodeid];
  vector<bool> binarystates(2, false);
  binarystates[0] = true; //true = in foreground
  bool anc = false; bool next_anc= false; bool prev = false; bool next= false;
  bool lc = false; bool rc = false; bool lc_prev =false; bool rc_prev= false;

  double tol = 1e-5;
  double tmp = 0; 
  vector<double> pu, pm, mar;
  //  double pu0 = min(max( tol, node_posterior[pos][nodeid]), 1- tol) ; // weight the posterior with old estimates???
  //  double pm0 = 1 - pu0;// weight the posterior with old estimates???

  double pu0 = 0.5, pm0=0.5;
  /*when initiating pu0 and pm0, use more info from the leaves*/
  if(clade_leaf_idx[nodeid].size()>1){
    pu0 = 0; pm0 = 0;
    for(size_t lr = 0; lr <1; ++lr){
      size_t clade_leaf_size = clade_leaf_idx[childidx[nodeid][lr]].size();
      for(size_t i = 0; i < clade_leaf_size; ++i){
        size_t lrleaf = clade_leaf_idx[childidx[nodeid][lr]][i];
        vector<size_t>::const_iterator it = std::find(leafidx.begin(), leafidx.end(), lrleaf);
        size_t whichleaf = std::distance(leafidx.begin(), it);
        betabin FG(Fparams[whichleaf].first, Fparams[whichleaf].second);
        betabin BG(Bparams[whichleaf].first, Bparams[whichleaf].second);
        double height = history.get_distance(nodeid, lrleaf);
        pu0 += get_single_firstprob(false, true, true, G, Q, height)*exp(FG(meth[whichleaf][pos]))/clade_leaf_size; //cur = leaf
        pu0 += get_single_firstprob(false, true, false, G, Q, height)*exp(BG(meth[whichleaf][pos]))/clade_leaf_size; //cur = leaf
        pm0 += get_single_firstprob(false, false, true, G, Q, height)*exp(FG(meth[whichleaf][pos]))/clade_leaf_size; //cur = leaf
        pm0 += get_single_firstprob(false, false, false, G, Q, height)*exp(BG(meth[whichleaf][pos]))/clade_leaf_size; //cur = leaf
      }
    }
    double psum = pu0+pm0; 
    assert (psum>0);
    pu0 = pu0/psum;
    pm0 = pm0/psum;
  }

  double pmtmp, putmp;
  double null_branch = 0;
  if(root && start){// cerr << "case1" << endl;
    /////////----Case 1----//////////
    for(size_t i = 0; i < 2; ++i){
      next = binarystates[i];
      tmp = (i==0)? node_posterior[pos+1][nodeid]:1-node_posterior[pos+1][nodeid];
      double nextp = min(max(tol, tmp), 1-tol);
      for(size_t j = 0; j < 2; ++j){
        lc = binarystates[j];
        tmp = (j==0)? node_posterior[pos][childidx[nodeid][0]]:(1-node_posterior[pos][childidx[nodeid][0]]); 
        double lcp = min(max(tol, tmp), 1- tol);
        for(size_t k = 0; k < 2; ++k){
          rc = binarystates[k];
          tmp = (k==0)? node_posterior[pos][childidx[nodeid][1]]:(1-node_posterior[pos][childidx[nodeid][1]]);
          double rcp =  min(max(tol, tmp), 1- tol);
          mar.push_back( nextp*lcp*rcp );
          putmp = pu0; pmtmp = pm0;
          putmp *= get_single_firstprob(false, true, lc, G, Q, branches[childidx[nodeid][0]]); //cur = lc
          pmtmp *= get_single_firstprob(false, false, lc, G, Q, branches[childidx[nodeid][0]]); //cur = lc
          putmp *= get_single_firstprob(false, true, rc, G, Q, branches[childidx[nodeid][0]]);//cur = rc
          pmtmp *= get_single_firstprob(false, false, rc, G, Q, branches[childidx[nodeid][1]]); //cur = rc
          putmp *= get_single_transprob(true, false, true, next, G, Q, null_branch);//cur = next
          pmtmp *= get_single_transprob(true, false, false, next, G, Q, null_branch);//cur = next
          pu.push_back(putmp); pm.push_back(pmtmp);
          assert(putmp + pmtmp > 0); 
          //if (root) cerr << pos << "\t" << pu0 << "\t" << pm0 << "\t" << putmp/(putmp+pmtmp) << "\t" << putmp/(putmp+pmtmp) << endl;
        }
      }
    }
  }else if (root && end){// cerr << "case2" << endl;
    /////////----Case 2----//////////
    for(size_t i = 0; i < 2; ++i){
      prev = binarystates[i]; 
      tmp = (i==0)? node_posterior[pos-1][nodeid]:1-node_posterior[pos-1][nodeid];
      double prevp = min(max(tol, tmp), 1-tol);
      for( size_t j = 0; j < 2; ++ j){ 
        lc_prev = binarystates[j];
        tmp = (j==0)? node_posterior[pos-1][childidx[nodeid][0]]:1-node_posterior[pos-1][childidx[nodeid][0]];
        double lc_prevp =  min(max(tol, tmp), 1-tol);
        for(size_t k = 0; k < 2; ++ k){
          rc_prev = binarystates[k];
          tmp = (k==0)? node_posterior[pos-1][childidx[nodeid][1]]:1-node_posterior[pos-1][childidx[nodeid][1]];
          double rc_prevp= min(max(tol, tmp), 1-tol);
          for(size_t l = 0; l < 2; ++l){
            lc = binarystates[l];
            tmp = (l==0)? node_posterior[pos][childidx[nodeid][0]]:1-node_posterior[pos][childidx[nodeid][0]];
            double lcp = min(max(tol, tmp), 1-tol);
            for(size_t m = 0; m <2; ++m){
              rc = binarystates[m];
              tmp = (m==0)? node_posterior[pos][childidx[nodeid][1]]:1-node_posterior[pos][childidx[nodeid][1]];
              double rcp = min(max(tol, tmp), 1-tol);
              mar.push_back( prevp*lc_prevp*rc_prevp*lcp*rcp);
              putmp = pu0; pmtmp = pm0;
              putmp *= get_single_transprob(true, false, prev, true, G, Q, null_branch);//cur=cur
              pmtmp *= get_single_transprob(true, false, prev, false, G, Q, null_branch);//cur=cur
              putmp *= get_single_transprob(false, true, lc_prev, lc, G, Q, branches[childidx[nodeid][0]]); //cur=lc
              pmtmp *= get_single_transprob(false, false, lc_prev, lc, G, Q, branches[childidx[nodeid][0]]); //cur=lc
              putmp *= get_single_transprob(false, true, rc_prev, rc, G, Q, branches[childidx[nodeid][1]]); //cur=rc
              pmtmp *= get_single_transprob(false, false, rc_prev, rc, G, Q, branches[childidx[nodeid][1]]); //cur=rc
              pu.push_back(putmp); pm.push_back(pmtmp);
              assert(putmp + pmtmp > 0); 
              //if (root) cerr << pos << "\t" << pu0 << "\t" << pm0 << "\t" << putmp/(putmp+pmtmp) << "\t" << putmp/(putmp+pmtmp) << endl;
            }
          }
        }
      }
    }
  }else if (root){ //cerr << "case3" << endl;
    /////////----Case 3----//////////
    for(size_t i = 0; i < 2; ++i){
      prev=binarystates[i]; 
      tmp = (i==0)? node_posterior[pos-1][nodeid]:1-node_posterior[pos-1][nodeid];
      double prevp = min(max(tol, tmp), 1-tol);
      for(size_t j = 0; j < 2; ++j){
        next = binarystates[j]; 
        tmp = (j==0)? node_posterior[pos+1][nodeid]:1-node_posterior[pos+1][nodeid];
        double nextp = min(max(tol, tmp), 1-tol);
        for(size_t k = 0; k < 2; ++k){
          lc = binarystates[k]; 
          tmp = (k==0)? node_posterior[pos][childidx[nodeid][0]]:1-node_posterior[pos][childidx[nodeid][0]];
          double lcp = min(max(tol, tmp), 1-tol);
          for(size_t l = 0; l < 2; ++l){
            rc = binarystates[l]; 
            tmp= (l==0)? node_posterior[pos][childidx[nodeid][1]]:1-node_posterior[pos][childidx[nodeid][1]];
            double rcp = min(max(tol, tmp), 1-tol);
            for(size_t m = 0; m < 2; ++m){
              lc_prev = binarystates[m]; 
              tmp = (m==0)? node_posterior[pos-1][childidx[nodeid][0]]:1-node_posterior[pos-1][childidx[nodeid][0]];
              double lc_prevp = min(max(tol, tmp), 1-tol);
              for(size_t n = 0; n < 2; ++n){
                rc_prev = binarystates[n]; 
                tmp = (n==0)? node_posterior[pos-1][childidx[nodeid][1]]:1-node_posterior[pos-1][childidx[nodeid][1]];
                double rc_prevp = min(max(tol, tmp), 1-tol);
                mar.push_back( exp(log(prevp)+log(nextp)+log(lc_prevp)+log(rc_prevp)+log(lcp)+log(rcp)));
                putmp = pu0; pmtmp = pm0;
                putmp *= get_single_transprob(true, false, prev, true, G, Q, null_branch);//cur=cur
                pmtmp *= get_single_transprob(true, false, prev, false, G, Q, null_branch);//cur=cur
                putmp *= get_single_transprob(true, false, true, next, G, Q, null_branch);//cur=next
                pmtmp *= get_single_transprob(true, false, false, next, G, Q, null_branch);//cur=next
                putmp *= get_single_transprob(false, true, lc_prev, lc, G, Q, branches[childidx[nodeid][0]]); //cur=lc
                pmtmp *= get_single_transprob(false, false, lc_prev, lc, G, Q, branches[childidx[nodeid][0]]); //cur=lc
                putmp *= get_single_transprob(false, true, rc_prev, rc, G, Q, branches[childidx[nodeid][1]]); //cur=rc
                pmtmp *= get_single_transprob(false, false, rc_prev, rc, G, Q, branches[childidx[nodeid][1]]); //cur=rc
                pu.push_back(putmp); pm.push_back(pmtmp);
                assert(putmp + pmtmp > 0); 
                // if (root) cerr << pos << "\t" << pu0 << "\t" << pm0 << "\t" << putmp/(putmp+pmtmp) << "\t" << pmtmp/(putmp+pmtmp) << endl;
          
              }
            }
          }
        }
      }
    }
  }else if (leaf && start){ //cerr<< "case4" << endl;
    /////////----Case 4----//////////
    for(size_t i = 0; i< 2; ++i){
      anc = binarystates[i]; 
      tmp = (i==0)? node_posterior[pos][parentidx[nodeid]]:1-node_posterior[pos][parentidx[nodeid]];
      double ancp = min(max(tol, tmp), 1- tol);
      for(size_t j =0; j <2; ++j){
        next = binarystates[j];
        tmp = (j==0)? node_posterior[pos+1][nodeid]:1-node_posterior[pos+1][nodeid];
        double nextp = min(max(tol, tmp), 1- tol);
        for(size_t k = 0; k <2; ++k ){
          next_anc = binarystates[k];
          tmp = (k==0)? node_posterior[pos+1][parentidx[nodeid]]:1-node_posterior[pos+1][parentidx[nodeid]];
          double next_ancp = min(max(tol, tmp), 1- tol);
          mar.push_back(ancp*nextp*next_ancp);
          vector<size_t>::const_iterator it = std::find(leafidx.begin(), leafidx.end(), nodeid);
          size_t whichleaf = std::distance(leafidx.begin(), it);
          betabin FG(Fparams[whichleaf].first, Fparams[whichleaf].second);
          betabin BG(Bparams[whichleaf].first, Bparams[whichleaf].second);
          putmp = pu0; pmtmp = pm0;
          putmp *= exp(FG(meth[whichleaf][pos]));
          pmtmp *= exp(BG(meth[whichleaf][pos]));
          putmp *= get_single_firstprob(false, anc, true, G, Q, branches[nodeid]);
          pmtmp *= get_single_firstprob(false, anc, false, G, Q, branches[nodeid]);
          putmp *= get_single_transprob(false, next_anc, true, next, G, Q, branches[nodeid]); 
          pmtmp *= get_single_transprob(false, next_anc, false, next, G, Q, branches[nodeid]); 
          pu.push_back(putmp); pm.push_back(pmtmp);
          assert(putmp + pmtmp > 0); 
        }
      }
    }
  }else if (leaf && end) { //cerr<< "case5" << endl;
    /////////----Case 5----//////////
    for(size_t i = 0; i<2; ++i){
      anc = binarystates[i]; 
      tmp = (i==0)? node_posterior[pos][parentidx[nodeid]]: 1 - node_posterior[pos][parentidx[nodeid]];
      double ancp = min(max( tol, tmp), 1- tol); 
      for(size_t j = 0; j <2 ; ++j){
        prev = binarystates[j]; 
        tmp  = (j==0)? node_posterior[pos-1][nodeid]: 1 - node_posterior[pos-1][nodeid];
        double prevp = min(max( tol, tmp), 1- tol);
        mar.push_back(ancp*prevp);
        vector<size_t>::const_iterator it = std::find(leafidx.begin(), leafidx.end(), nodeid);
        size_t whichleaf = std::distance(leafidx.begin(), it);
        betabin FG(Fparams[whichleaf].first, Fparams[whichleaf].second);
        betabin BG(Bparams[whichleaf].first, Bparams[whichleaf].second);
        putmp = pu0; pmtmp = pm0;
        putmp *= exp(FG(meth[whichleaf][pos]));
        pmtmp *= exp(BG(meth[whichleaf][pos]));
        putmp *= get_single_transprob(false, anc, prev, true, G, Q, branches[nodeid]); //cur=cur 
        pmtmp *= get_single_transprob(false, anc, prev, false, G, Q, branches[nodeid]); //cur=cur
        pu.push_back(putmp); pm.push_back(pmtmp);
        assert(putmp + pmtmp > 0);  
     }
    }
  }else if (leaf){ //cerr<< "case6" << endl;
    /////////----Case 6----//////////
    for(size_t i = 0; i < 2; ++i){
      prev = binarystates[i];
      tmp = (i==0)? node_posterior[pos-1][nodeid]: 1 - node_posterior[pos-1][nodeid];
      double prevp = min(max(tol, tmp), 1- tol);
      for(size_t j = 0; j < 2; ++j){
        anc = binarystates[j];
        tmp = (j==0)? node_posterior[pos][parentidx[nodeid]]:1 - node_posterior[pos][parentidx[nodeid]];
        double ancp = min(max(tol, tmp), 1- tol);
        for(size_t k = 0; k < 2; ++k){
          next = binarystates[k];
          tmp = (k==0)? node_posterior[pos+1][nodeid]:1 - node_posterior[pos+1][nodeid];
          double nextp = min(max(tol, tmp), 1- tol);
          for(size_t l = 0; l < 2; ++l){
            next_anc = binarystates[l]; 
            tmp = (l==0)? node_posterior[pos][parentidx[nodeid]]:1 - node_posterior[pos][parentidx[nodeid]];
            double next_ancp = min(max(tol, tmp), 1- tol);
            mar.push_back(ancp*prevp*next_ancp*nextp);
            vector<size_t>::const_iterator it;
            it = std::find(leafidx.begin(), leafidx.end(), nodeid);
            size_t whichleaf = std::distance(leafidx.begin(), it);
            betabin FG(Fparams[whichleaf].first, Fparams[whichleaf].second);
            betabin BG(Bparams[whichleaf].first, Bparams[whichleaf].second);
            putmp = pu0; pmtmp = pm0;
            putmp *= exp(FG(meth[whichleaf][pos]));
            pmtmp *= exp(BG(meth[whichleaf][pos]));
            putmp *= get_single_transprob(false, anc, prev, true, G, Q, branches[nodeid]); //cur=cur 
            pmtmp *= get_single_transprob(false, anc, prev, false, G, Q, branches[nodeid]); //cur=cur
            putmp *= get_single_transprob(false, next_anc, true, next, G, Q, branches[nodeid]); //cur=next 
            pmtmp *= get_single_transprob(false, next_anc, false, next, G, Q, branches[nodeid]); //cur=next
            pu.push_back(putmp); pm.push_back(pmtmp);
            assert(putmp + pmtmp > 0); 
          }
        }
      }
    }
  }else if (start){ //cerr << "case7" <<endl;
    /////////----Case 7----//////////
    for(size_t i = 0; i < 2; ++i){
      anc = binarystates[i];
      tmp = (i==0)? node_posterior[pos][parentidx[nodeid]]: 1 -  node_posterior[pos][parentidx[nodeid]]; 
      double ancp = min(max(tol, tmp), 1-tol);
      for(size_t j = 0; j < 2; ++j){
        next = binarystates[j];
        tmp = (j==0)? node_posterior[pos+1][nodeid]: 1- node_posterior[pos+1][nodeid];
        double nextp = min(max(tol, tmp), 1-tol);
        for(size_t k = 0; k < 2; ++k ){
          next_anc = binarystates[k];
          tmp = (k==0)? node_posterior[pos+1][parentidx[nodeid]]:1 - node_posterior[pos+1][parentidx[nodeid]];
          double next_ancp = min(max(tol, tmp), 1-tol);
          for(size_t l = 0; l < 2; ++l){
            lc = binarystates[l];
            tmp = (l==0)? node_posterior[pos][childidx[nodeid][0]]:1-node_posterior[pos][childidx[nodeid][0]];
            double lcp = min(max(tol, tmp), 1-tol);
            for(size_t m = 0; m < 2; ++m){
              rc = binarystates[m];
              tmp = (m ==0) ? node_posterior[pos][childidx[nodeid][1]]:1-node_posterior[pos][childidx[nodeid][1]];
              double rcp =  min(max(tol, tmp), 1-tol);
              mar.push_back(ancp*nextp*next_ancp*lcp*rcp);
              putmp = pu0; pmtmp = pm0;
              putmp *= get_single_firstprob(false, anc, true, G, Q, branches[nodeid]);
              pmtmp *= get_single_firstprob(false, anc, false, G, Q, branches[nodeid]);
              putmp *= get_single_firstprob(false, true, lc, G, Q, branches[childidx[nodeid][0]]);//cur=lc
              pmtmp *= get_single_firstprob(false, false, lc, G, Q, branches[childidx[nodeid][0]]);//cur=lc
              putmp *= get_single_firstprob(false, true, rc, G, Q, branches[childidx[nodeid][1]]); //cur=rc
              pmtmp *= get_single_firstprob(false, false, rc, G, Q, branches[childidx[nodeid][1]]);//cur=rc
              putmp *= get_single_transprob(false, next_anc, true, next, G, Q, branches[nodeid]); //cur=next 
              pmtmp *= get_single_transprob(false, next_anc, false, next, G, Q, branches[nodeid]); //cur=next
              pu.push_back(putmp); pm.push_back(pmtmp); 
              assert(putmp + pmtmp > 0); 
            }
          }
        }
      }
    }
  }else if (end){// cerr << "case8"  <<endl; 
    /////////----Case 8----//////////
    for(size_t i = 0; i <2; ++ i){
      anc = binarystates[i];
      tmp = (i==0)? node_posterior[pos][parentidx[nodeid]] : 1 - node_posterior[pos][parentidx[nodeid]];
      double ancp = min(max( tol, tmp), 1- tol);
      for(size_t j = 0; j <2; ++ j){
        prev = binarystates[j];
        tmp = (j==0)? node_posterior[pos-1][nodeid] : 1 - node_posterior[pos-1][nodeid]; 
        double prevp = min(max(tol, tmp), 1- tol);
        for(size_t k = 0; k < 2; ++ k){
          lc = binarystates[k];
          tmp = (k==0)? node_posterior[pos][childidx[nodeid][0]]: 1 - node_posterior[pos][childidx[nodeid][0]];
          double lcp = min(max(tol, tmp), 1- tol);
          for(size_t l = 0; l < 2; ++l){
            rc = binarystates[l];
            tmp = (l==0)? node_posterior[pos][childidx[nodeid][1]]: 1 - node_posterior[pos][childidx[nodeid][1]]; 
            double rcp = min(max(tol, tmp), 1- tol);
            for(size_t m = 0; m < 2; ++ m){
              lc_prev = binarystates[m];
              tmp = (m==0)? node_posterior[pos-1][childidx[nodeid][0]]: 1 - node_posterior[pos-1][childidx[nodeid][0]];  
              double lc_prevp = min(max(tol, tmp), 1- tol); 
              for(size_t n = 0; n < 2; ++ n){
                rc_prev = binarystates[n];
                tmp = (n==0)? node_posterior[pos-1][childidx[nodeid][1]]: 1 - node_posterior[pos-1][childidx[nodeid][1]];
                double rc_prevp = min(max(tol, tmp), 1-tol);
                mar.push_back(prevp*ancp*lcp*rcp*lc_prevp*rc_prevp);
                putmp = pu0; pmtmp = pm0;
                putmp *= get_single_transprob(false, anc, prev, true, G, Q, branches[nodeid]); //cur=cur 
                pmtmp *= get_single_transprob(false, anc, prev, false, G, Q, branches[nodeid]); //cur=cur
                putmp *= get_single_transprob(false, true, lc_prev, lc, G, Q, branches[childidx[nodeid][0]]); //cur=lc 
                pmtmp *= get_single_transprob(false, false, lc_prev, lc, G, Q, branches[childidx[nodeid][0]]); //cur=lc
                putmp *= get_single_transprob(false, true, rc_prev, rc, G, Q, branches[childidx[nodeid][1]]); //cur=rc 
                pmtmp *= get_single_transprob(false, false, rc_prev, rc, G, Q, branches[childidx[nodeid][1]]); //cur=rc
                pu.push_back(putmp); pm.push_back(pmtmp);   
                assert(putmp + pmtmp > 0);     
              }
            }
          }
        }
      }
    } 
  }else{ //cerr << "case9" << endl;
    /////////----Case 9----//////////
    for(size_t i = 0; i <2; ++ i){
      prev = binarystates[i];
      tmp = (i==0)? node_posterior[pos-1][nodeid]: 1-node_posterior[pos-1][nodeid];
      double prevp = min(max(tol, tmp), 1-tol); 
      for(size_t j = 0; j < 2; ++ j){
        anc = binarystates[j];
        tmp = (j==0)? node_posterior[pos][parentidx[nodeid]] : 1 - node_posterior[pos][parentidx[nodeid]]; 
        double ancp= min(max(tol, tmp), 1-tol);  
        for(size_t k = 0; k < 2; ++ k){
          next = binarystates[k];
          tmp = (k==0)? node_posterior[pos+1][nodeid]: 1 - node_posterior[pos+1][nodeid];
          double nextp = min(max(tol, tmp), 1-tol);  
          for(size_t l = 0; l < 2; ++l){
            next_anc = binarystates[l];
            tmp = (l==0)? node_posterior[pos+1][parentidx[nodeid]]: 1 - node_posterior[pos+1][parentidx[nodeid]]; 
            double next_ancp = min(max(tol, tmp), 1-tol);  
            for(size_t m = 0; m < 2; ++ m){
              lc = binarystates[m];
              tmp = (m==0)? node_posterior[pos][childidx[nodeid][0]]: 1 - node_posterior[pos][childidx[nodeid][0]]; 
              double lcp = min(max(tol, tmp), 1-tol);  
              for(size_t n = 0; n < 2; ++ n){
                rc = binarystates[n];
                tmp = (n==0)? node_posterior[pos][childidx[nodeid][1]]: 1 - node_posterior[pos][childidx[nodeid][1]];
                double rcp = min(max(tol, tmp), 1-tol);  
                for(size_t p = 0; p < 2; ++p){
                  lc_prev = binarystates[p];
                  tmp = (p==0)? node_posterior[pos-1][childidx[nodeid][0]]: 1 - node_posterior[pos-1][childidx[nodeid][0]];
                  double lc_prevp  = min(max(tol, tmp), 1-tol);
                  for(size_t q = 0; q < 2; ++ q){
                    rc_prev = binarystates[q]; 
                    tmp = (q==0)? node_posterior[pos-1][childidx[nodeid][1]]:1-node_posterior[pos-1][childidx[nodeid][1]];
                    double rc_prevp = min(max(tol, tmp), 1-tol); 
                    mar.push_back(prevp*ancp*nextp*next_ancp*lcp*rcp*lc_prevp*rc_prevp);
                    putmp = pu0; pmtmp =pm0;
                    putmp *= get_single_transprob(false, anc, prev, true, G, Q, branches[nodeid]); //cur=cur 
                    pmtmp *= get_single_transprob(false, anc, prev, false, G, Q, branches[nodeid]); //cur=cur
                    putmp *= get_single_transprob(false, next_anc, true, next, G, Q, branches[nodeid]); //cur=next
                    pmtmp *= get_single_transprob(false, next_anc, false, next, G, Q, branches[nodeid]); //cur=next
                    putmp *= get_single_transprob(false, true, lc_prev, lc, G, Q, branches[childidx[nodeid][0]]); //cur=lc
                    pmtmp *= get_single_transprob(false, false, lc_prev, lc, G, Q, branches[childidx[nodeid][0]]); //cur=lc
                    putmp *= get_single_transprob(false, true, rc_prev, rc, G, Q, branches[childidx[nodeid][1]]); //cur=rc
                    pmtmp *= get_single_transprob(false, false, rc_prev, rc, G, Q, branches[childidx[nodeid][1]]); //cur=rc
                    pu.push_back(putmp); pm.push_back(pmtmp); 
                    assert(putmp + pmtmp > 0); 
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  double newscore = 0; 
  for(size_t i = 0; i < pu.size(); ++i){
    newscore += pu[i]/(pu[i]+pm[i])*mar[i]; 
    assert(newscore >0);
  }
  node_posterior[pos][nodeid] = newscore;
  return (std::fabs(newscore - oldscore));
}


void 
approx_posterior( const vector<vector<pair<double, double> > > &meth,
                  const vector<size_t> &reset_points, const HistoryME &history, 
                  const size_t MAXITER,
                  vector<vector<double> > &node_posterior){ 
  //posterior of being in foreground
  double TOL = 1e-2; 
  size_t nleaf = history.get_nleaf(); 
  size_t ns = history.get_nspecies();
  size_t nsites = meth[0].size();
  vector<size_t> leafidx;
  history.get_leaf_idx(leafidx);
  vector<vector<size_t> > childidx;
  history.get_child_idx(childidx);
  vector<size_t> paidx;
  history.get_parent_idx(paidx);
  vector<double> G, Q;
  history.get_G(G);
  history.get_Q(Q);
  vector<double> branches;
  history.get_branches(branches);
  vector<vector<size_t> > clade_leaf_idx;
  history.get_clade_leaf_idx(clade_leaf_idx);

  node_posterior = vector< vector<double > >(nsites, vector<double>(ns, 0.0)) ; 
  vector<pair<double, double> > Fparams;
  vector<pair<double, double> > Bparams;
  for(size_t i = 0; i < nleaf; ++i){
    Fparams.push_back(history.get_Fparam(i));
    Bparams.push_back(history.get_Bparam(i));
  }
  //(A) Initialize leaf posteriors
  for(size_t i = 0; i < nleaf; ++i){
    betabin FG(Fparams[i].first, Fparams[i].second);
    betabin BG(Fparams[i].first, Fparams[i].second);
    for(size_t j = 0; j < nsites; ++j){
      node_posterior[j][leafidx[i]] = 1.0/(1.0 + exp(BG(meth[i][j])-FG(meth[i][j])));
    }
  }

  vector<double> pi(2,0);
  pi[0] = (1-G[1])/(2-G[0]-G[1]);
  pi[1] = (1-G[0])/(2-G[0]-G[1]);
  
  //(B) Initialize internal posteriors
  for(size_t i =0; i <nsites; ++i){
    vector<size_t> wait4update;
    for(size_t j = 0; j < ns; ++j){
      if(childidx[j].size() > 0){
        if(node_posterior[i][childidx[j][0]] !=0 && 
           node_posterior[i][childidx[j][1]] !=0 ){
          node_posterior[i][j] = single_update(j, i, childidx,G, Q,branches, node_posterior);
        }else{
          wait4update.push_back(j);
        }
      }
    }
    while(wait4update.size()>0){
      size_t j = wait4update.back();
      wait4update.pop_back();
      node_posterior[i][j] = single_update(j, i, childidx,G, Q, branches, node_posterior);
    }
  }

  //(C) Iterate until stable
  bool root = false, leftend = false, rightend = false, isleaf = false;
  double delta;
  size_t iteration;
  for(size_t i = 0; i < reset_points.size()-1; ++i){
    cerr << ".";
    iteration = 0;
    size_t start = reset_points[i];
    size_t end = reset_points[i+1];
    do{
      ++iteration ;
      delta= 0;
      for(size_t pos = start; pos < end; ++pos){
        for(size_t nodeid = 0; nodeid < ns; ++nodeid){
          root = (nodeid==0)? true: false; 
          leftend = (pos==start)?true:false; 
          rightend = (pos==end-1)?true:false;
          isleaf = (childidx[nodeid].size()>0)?false:true;
          delta += iterate_update(root, isleaf, leftend, rightend, pos, nodeid, 
                                  childidx, paidx, leafidx, clade_leaf_idx, G, Q, 
                                  branches, Fparams,Bparams, history, meth, node_posterior);
        }
      }
      // cerr << endl << iteration << "\t" << delta/(end-start)/ns << endl;
    }while(delta > TOL*(end-start)*ns && iteration < MAXITER);
    if(delta < TOL*(end-start)*ns)  cerr << "Converged at iteration" << iteration << endl;
    else cerr << delta << " hasn't reached the convergence criteria " << TOL*(end-start)*ns << endl; 
  }
}

void 
candidate_hme(const size_t NHME,
              const vector<vector<double> > &node_posterior,
              vector<vector<vector<bool> > > &hmecand,
              vector<vector<double> > &hmecand_ps){
  hmecand_ps.clear();
  hmecand.clear();
  size_t ns = node_posterior[0].size();
  size_t nsites = node_posterior.size();
  assert(NHME < std::pow(2, double(ns)));

  for(size_t i =0; i < nsites; ++i){
    vector<val_ind> mapscore;
    //sort by max(FG, BG) posterior score in ascending order
    for(size_t j =0; j < ns; ++j){
      val_ind tmp;
      tmp.val = node_posterior[i][j];
      tmp.ind = j;
      if(tmp.val < 0.5) 
        tmp.val = 1-tmp.val;
      mapscore.push_back(tmp);
    }
    std::sort(mapscore.begin(), mapscore.end(), by_val());
    //Adding HMEs one by one
    vector<vector<bool> > topHMEs;
    vector<double> top_scores;
    size_t alt = 0; 
    while(topHMEs.size()<NHME){
      vector<bool> hmetmp(ns, false); 
      double scoretmp=1;
      if (topHMEs.size()==0){
        for(size_t j = 0; j < ns; ++j){
          scoretmp *= mapscore[j].val;
          size_t sp = mapscore[j].ind ;
          if (node_posterior[i][sp] > 0.5) //foreground posterior > 0.5
            hmetmp[sp] = true; 
        }
        topHMEs.push_back(hmetmp);
        top_scores.push_back(scoretmp);
      }else{
        size_t cur_nhme = topHMEs.size();
        for(size_t j=0; j < min(cur_nhme, NHME-cur_nhme); ++j){
          vector<bool> hmetmp=topHMEs[j];
          size_t sp = mapscore[alt].ind ;
          hmetmp[sp] = !(hmetmp[sp]); 
          topHMEs.push_back(hmetmp);
          top_scores.push_back( top_scores[j]*(1- mapscore[alt].val)/( mapscore[alt].val));
        }
        ++alt;
      }
    }
    hmecand.push_back(topHMEs);
    hmecand_ps.push_back(top_scores);
    //for(size_t k =0; k < top_scores.size(); ++k) cerr << top_scores[k] << "\t";
    //cerr << endl;
  }
}


//ksi_pos(i,j)
void 
trans_posterior(const HistoryME &history,
                const size_t NHME,
                const vector<vector<vector<bool> > > &hmecand,
                const vector<vector<double> > &hmecand_ps,
                const size_t pos,
                vector<vector<double> > &constrained_trans_post ){
  constrained_trans_post  = vector<vector<double> >(NHME, vector<double>(NHME, 0)); 
  PhyloTree t(history.get_Newick_format() );
  double total=0;
  for(size_t idx1 = 0; idx1 < NHME; ++idx1){
    size_t i = history.get_HME_idx(hmecand[pos][idx1]);
    for(size_t idx2 = 0; idx2 < NHME; ++idx2){
      size_t j = history.get_HME_idx(hmecand[pos+1][idx2]);
      
      //find out which HMEs correspond to hmecand[idx1] and   hmecand[idx2]
      double aij = history.get_logTPM_element(i, j) ; 
      double tp = log(hmecand_ps[pos][idx1]) + aij + log(hmecand_ps[pos+1][idx2]);
      if(idx1==0 && idx2==0) total = tp;
      else  total = log_sum_log(total, tp);
      constrained_trans_post[idx1][idx2] = tp;
    }
  }
  for(size_t idx1 = 0; idx1 < NHME; ++idx1){
    for(size_t idx2 = 0; idx2 < NHME; ++idx2){  
      constrained_trans_post[idx1][idx2] = exp(constrained_trans_post[idx1][idx2] -total);
    }
  }
}


/******************************************************************
 *Convert posterior scores to coefficients in objective function  *
 *A version for approximation                                     *
 *Similar to the full version in EpiPhyloHMM.cpp                  * 
 */
void
post_scores_to_coeff_app(const bool VERBOSE,
                         const vector<vector<pair<double, double> > > &values,
                         const std::vector<size_t> &reset_points,
                         const HistoryME &history,
                         const size_t NHME,
                         const vector<vector<vector<bool> > > &hmecand,
                         const vector<vector<double> > &hmecand_ps,
                         const vector<bool> &unknown,
                         vector<double> &coeff){
  if(VERBOSE)
    cerr << "BEGIN---post_scores_to_coeff----[" << NHME << "] top hmes---" << endl;

  size_t n = history.get_nleaf();
  size_t ns = history.get_nspecies();
  size_t nb = ns-1; //#branches
  size_t a, b, c; 

  PhyloTree t(history.get_Newick_format());
  vector<size_t> pa_idx;
  history.get_parent_idx(pa_idx);
  size_t idx_child, idx_parent;

  vector<double> coeff2(4*nb, 0);
  vector<double> coeff1(4, 0);
  vector<double> coeff3(8*nb, 0);
  vector<double> coeff4(4*nb, 0);
  for(size_t i = 0; i < reset_points.size()-1; ++i){
    size_t start = reset_points[i];
    size_t end = reset_points[i+1];
    if(VERBOSE)
      cerr << "In block " << start << " to " << end << endl;

    for(size_t i = 0; i < NHME; ++ i){
      size_t hmeidx = history.get_HME_idx(hmecand[start][i]);
      for(size_t idx_child = 1; idx_child < ns; ++ idx_child){ 
        idx_parent = pa_idx[idx_child]; 
        assert (idx_parent < ns);
        a = history.get_state(hmeidx, idx_parent)? 0:1;
        b = history.get_state(hmeidx, idx_child)? 0:1;
        coeff2[(idx_child-1)*4 + a*2 + b] += hmecand_ps[start][i];
        assert(hmecand_ps[start][i] < 1 && hmecand_ps[start][i] >0);
      }
    }

    for(size_t pos = start; pos < end-1; ++pos){ //genomic locations
      if (VERBOSE && pos %1000 ==0)
        cerr << "=" ;

      vector<vector<double> > constrained_trans_post;
      trans_posterior(history, NHME, hmecand, hmecand_ps, pos, constrained_trans_post);

      for(size_t i = 0; i < NHME; ++i){
        size_t hmeidx1 = history.get_HME_idx(hmecand[pos][i]);
        size_t cur = history.get_state(hmeidx1, 0)? 0:1;

        for(size_t k = 0; k < NHME; ++k){
          size_t hmeidx2 = history.get_HME_idx(hmecand[pos][k]);
          size_t next = history.get_state(hmeidx2, 0) ? 0:1;
          double tp = constrained_trans_post[i][k];
          coeff1[cur*2+next] += tp;
          assert(tp > 0 && tp < 1);

          for(size_t j = 1; j < ns; ++j){ //branches
            idx_child = j;
            idx_parent = pa_idx[idx_child];
            assert (idx_parent < ns);
            a = history.get_state(hmeidx1, idx_child)? 0:1;
            b = history.get_state(hmeidx2, idx_parent)? 0:1;
            c = history.get_state(hmeidx2, idx_child)? 0:1;
            coeff3[(j-1)*8 + a*4 + b*2 + c] += tp; 
            coeff4[(j-1)*4 + a*2 + b] += tp; 
          }
        }
      }
    }
  }
  if(VERBOSE)
    cerr << endl;

  coeff.clear();
  coeff.push_back(static_cast<double>(n));
  coeff.push_back(static_cast<double>(nb));
  coeff.insert(coeff.end(), coeff1.begin(), coeff1.end());
  coeff.insert(coeff.end(), coeff2.begin(), coeff2.end());
  coeff.insert(coeff.end(), coeff3.begin(), coeff3.end());
  coeff.insert(coeff.end(), coeff4.begin(), coeff4.end());
  for(size_t i = 0; i < unknown.size();++i){
    if( unknown[i]) coeff.push_back(0);
    else coeff.push_back(1);
  }
  if(VERBOSE){
    cerr << "------post_scores_to_coeff-------------END" << endl
         << endl << "---[coeff1]---" << endl;
    for (size_t i = 0; i < coeff1.size(); ++i){
      cerr << coeff1[i] << "\t";
    }
    cerr << endl << "---[coeff2]---" << endl;
    for (size_t i = 0; i < coeff2.size(); ++i){
      cerr << coeff2[i] << "\t";
      if ((i+1)%nb ==0) cerr << endl;
    }
    cerr << endl << "---[coeff3]---" << endl;
    for (size_t i = 0; i < coeff3.size(); ++i){
      cerr << coeff3[i] << "\t";
      if ((i+1)%nb ==0) cerr << endl;
    }
    cerr << endl << "---[coeff4]---" << endl;
    for (size_t i = 0; i < coeff4.size(); ++i){
      cerr << coeff4[i] << "\t";
      if ((i+1)%nb ==0) cerr << endl;
    }
  }
}
