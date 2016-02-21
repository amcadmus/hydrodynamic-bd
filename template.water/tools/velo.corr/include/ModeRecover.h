#pragma once

#include <vector>
#include "ModeOverlapping.h"
#include "Integral1D.h"

using namespace std;

class ModeRecover 
{
public:
  ModeRecover (const vector<double> & kks,
	       const double & h_w,
	       const double cond_numb = 1e6);
  void reinit (const vector<double> & kks,
	       const double & h_w);
  void solve (const vector<double > & rhs,
	      vector<double > & out) const;
private:
  ModeOverlapping ovlp;
  Integral1D<ModeOverlapping,double > inte1d;
  unsigned nn;
  unsigned start_posi;
  double cond_numb;
  vector<double > AA;
  vector<double > eig;
private:
  double computeMatEle (const double & kk,
			const double & xi,
			const double & hw ) ;
  double computeMatEle2 (const double & kk,
			 const double & xi,
			 const double & hw ) ;
  void diagonalize ();
  void apply_Tt (vector<double > & vec) const;
  void apply_T  (vector<double > & vec) const;
  void print_eig () const;
  void print_AA  (ostream & fp = cout) const;
  void print_vec (const vector<double > & vec) const;
};


