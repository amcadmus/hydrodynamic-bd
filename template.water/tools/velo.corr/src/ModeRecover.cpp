#include "ModeRecover.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include "mkl_lapack.h"
#include <stdlib.h>
#include <fstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

double
ModeRecover::
computeMatEle (const double & kk,
	       const double & xi,
	       const double & hw ) 
{
  double period_kk = 2. * M_PI / kk;
  double period_xi = 2. * M_PI / xi;
  double n_kk = 2 * hw / period_kk;
  double n_xi = 2 * hw / period_xi;
  int nx = (n_kk < n_xi ? n_xi : n_kk);
  if (nx == 0) nx = 1;
  ovlp.reinit (kk, xi);
  return inte1d.cal_int (2,
			 ovlp,
			 -hw, hw,
			 1e-16,
			 nx);
}

double
ModeRecover::
computeMatEle2 (const double & kk_,
		const double & xi_,
		const double & hw ) 
{
  if (kk_ != xi_){
    double kk, xi;
    if (kk_ < xi_) {
      kk = xi_;
      xi = kk_;
    }
    else {
      kk = kk_;
      xi = xi_;
    }
    double diff = kk - xi;
    double summ = kk + xi;
    return sin (diff * hw) / diff - sin(summ * hw) / summ;
  }
  else{
    double summ = kk_ + xi_;
    return hw - sin(summ * hw) / summ;
  }
}

ModeRecover::
ModeRecover (const vector<double> & kks,
	     const double & h_w,
	     const double cond_numb_)
    : nn(0), cond_numb (cond_numb_)
{
  reinit (kks, h_w);
}

void 
ModeRecover::
reinit (const vector<double> & kks,
	const double & h_w)
{
  nn = kks.size();
  AA.resize (nn * nn);

  for (unsigned ii = 0; ii < nn; ++ii){
    for (unsigned jj = ii; jj < nn; ++jj){
      AA[ii*nn+jj] = computeMatEle2 (kks[ii], kks[jj], h_w);
      // cout << AA[ii*nn+jj] << "\t"
      // 	   << computeMatEle2 (kks[ii], kks[jj], h_w) << "\t"
      // 	   << fabs (AA[ii*nn+jj] - computeMatEle2 (kks[ii], kks[jj], h_w))
      // 	   << endl;
      if (jj > ii){
      	AA[jj*nn+ii] = AA[ii*nn+jj];
      }
    }
  }

  ofstream a_out;
  a_out.open ("AA.out");
  print_AA(a_out);
  a_out.close ();

  diagonalize ();
  
  print_eig();
  // print_AA();

  start_posi = 0;
  for (int ii = nn-2; ii >= 0; ii --){
    if ( fabs (eig.back() / eig[ii]) > cond_numb ) {
      start_posi = ii+1;
      break;
    }
  }
  cout << "# matrix size " << nn 
       << " for cond_numb " << cond_numb
       // << " start_posi is " << start_posi
       << " used numb eig is " << nn - start_posi
       << endl;
}

void
ModeRecover::
solve (const vector<double > & rhs,
       vector<double > & out) const
{
  out = rhs;
  apply_Tt (out);
  for (unsigned ii = 0; ii < start_posi; ++ii){
    out[ii] = 0;
  }
  for (unsigned ii = 0; ii < nn; ++ii){
    out[ii] /= eig[ii];
  }
  apply_T (out);
  // cout << "in" << endl;
  // print_vec (rhs);
  // cout << "out" << endl;
  // print_vec (out);
}

void ModeRecover::
diagonalize () 
{
  eig.resize (nn);
  char eig_v = 'V';
  char lower = 'L';
  int int_n = nn;
  int lwork = nn * (3 + nn/2);
  double * work = new double [lwork];
  int info;
  dsyev_ (&eig_v, &lower, &int_n, &AA[0], &int_n, &eig[0], &work[0], &lwork, &info);
  delete [] work;
  if (info < 0){
    cerr << "illegal " << -info << "-th argument of dsyev " << endl;
    exit (1);
  }
  else if (info > 0){
    cerr << "dsyev fails to converge " << info << "-th diagonal element " << endl;
    exit (1);
  }
}

void ModeRecover::
apply_Tt (vector<double > & vec) const
{
  vector<double > tmp(vec);
  for (unsigned ii = start_posi; ii < nn; ++ii){
    // vec[ii] = 0;
    double sum = 0;
// #pragma omp parallel for reduction (+:sum)
    for (unsigned jj = 0; jj < nn; ++jj){
      sum += AA[ii*nn+jj] * tmp[jj];
    }
    vec[ii] = sum;
  }
}

void ModeRecover::
apply_T (vector<double > & vec) const
{
  vector<double > tmp(vec);
// #pragma omp parallel for
  for (unsigned ii = 0; ii < nn; ++ii){
    vec[ii] = 0;
  }
  
  for (unsigned jj = start_posi; jj < nn; ++jj){
// #pragma omp parallel for
    for (unsigned ii = 0; ii < nn; ++ii){
      vec[ii] += AA[jj*nn+ii] * tmp[jj];
    }
  }
}

void ModeRecover::
print_AA (ostream & fp) const
{
  fp << setprecision (16) << setfill (' ');
  for (unsigned ii = 0; ii < nn; ++ii){
    for (unsigned jj = 0; jj < nn; ++jj){
      fp << setw(18) << AA[ii*nn+jj] << " " ;
    }
    fp << endl;
  }
}

void ModeRecover::
print_eig () const
{
  cout << setprecision (16) << setfill (' ');
  for (unsigned ii = 0; ii < eig.size(); ++ii){
    cout << setw(18) << eig[ii] << " " ;
    cout << endl;
  }
}

void ModeRecover::
print_vec (const vector<double > & vec) const
{
  cout << setprecision (16) << setfill (' ');
  for (unsigned ii = 0; ii < vec.size(); ++ii){
    cout << setw(18) << vec[ii] << " " ;
    cout << endl;
  }
}
