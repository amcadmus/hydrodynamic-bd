#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <cmath>
#include <time.h>

#include "GroFileManager.h"
#include <boost/program_options.hpp>
#include "Trajectory.h"
#include "BlockAverage.h"
#include "AutoCorrelCalculator.h"

namespace po = boost::program_options;

using namespace std;

int main (int argc, char * argv[])
{
  std::string ifile;
  std::string ofile;
  double kkStart, kkStep, kkEnd;
  po::options_description desc ("Allow options");
  int index;
  
  desc.add_options()
      ("help,h", "print this message")
      ("mode-start", po::value<double > (&kkStart)->default_value (0.1, "the start of mode k"))
      ("mode-step", po::value<double > (&kkStep)->default_value (0.1, "the step of mode k"))
      ("mode-end", po::value<double > (&kkEnd)->default_value (1.0, "the end of mode k"))
      ("index,n", po::value<int > (&index)->default_value (1, "the index of the data point used to interpolate"))
      ("output,o",  po::value<std::string > (&ofile)->default_value (std::string("slope.corr.out"), "output of integrated corr as a function of modes"))
      ("input,f",  po::value<std::string > (&ifile)->default_value (std::string("corr"), "output of corr, stored by modes"));  

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help") || vm.count("h")){
    std::cout << desc<< "\n";
    return 0;
  }

  unsigned nkk = int ((kkEnd - kkStep + kkStep * 0.5) / kkStep) + 1;
  vector<double > kks (nkk);
  for (unsigned ii = 0; ii < nkk; ++ii) {
    kks[ii] = kkStart + ii * kkStep;
    // cout << kks[ii] << endl;
  }

  if (index <= 0){
    cerr << "wrong index number " << index << " should be larger than 0 " << endl;
    return 1;
  }
  
  FILE * fpo = fopen (ofile.c_str(), "w");
  // int lastStop = 0;
  for (unsigned ii = 0; ii < nkk; ++ii){
    // fprintf (fpo, "%f %e\n", kks[ii], inte[ii]);
    // double myiup = (kks.back() - kks[ii]) / (kks.back() - kks.front()) * (iup - ilo) + ilo;
    // cout << "kk is " << kks[ii] << " myiup is " << myiup << endl;
    char modei[32], modef[32];
    sprintf (modei, "%02d", int(kks[ii]));
    sprintf (modef, "%04d", int(10000 * (kks[ii] - int(kks[ii]) + 0.00005)));
    string filename = ifile + string (modei) + string(".") + string(modef) + string(".out");
    FILE * fp = fopen (filename.c_str(), "r");
    if (fp == NULL){
      cerr << "cannot open file " << filename << endl;
      return 1;
    }
    vector<double > times;
    vector<double > values;
    vector<double > errors;
    double tmp0, tmp1, tmp2;
    while (3 == fscanf (fp, "%lf %lf %lf", &tmp0, &tmp1, &tmp2)){
      times.push_back (tmp0);
      values.push_back (tmp1);
      errors.push_back (tmp2);
    }
    double sum = 0.;
    // bool reachEnd = false;
    // int jj;
    // for (jj = 0; jj < int(times.size()) - 1; ++jj){
    //   if (reachEnd) {
    // 	// cout << "jj " << times[jj] << endl;
    // 	if (lastStop == 0){
    // 	  lastStop = jj;
    // 	  cout << "file " << filename << " end reached at " << times[jj] << " break at " << times[jj]  << endl;
    // 	  break;
    // 	}
    // 	else if (jj < lastStop && fabs (lastStop - jj) <= 1.5) {
    // 	  lastStop = jj;
    // 	  cout << "file " << filename << " end reached at " << times[jj] << " break at " << times[jj]  << endl;
    // 	  break;
    // 	}
    // 	else if (jj >= lastStop) {
    // 	  lastStop = jj;
    // 	  cout << "file " << filename << " end reached at " << times[jj] << " break at " << times[jj]  << endl;
    // 	  break;
    // 	}
    //   }
    //   if (!reachEnd && fabs(values[jj]) - errors[jj] < 0){
    // 	reachEnd = true;
    //   }
    //   sum += 0.5 * (values[jj] + values[jj+1]) * (times[jj+1]  - times[jj]);
    // }
    // cout << "jj " << times[jj] << endl;
    // if (!reachEnd){
    //   lastStop = times.size() - 2;
    //   cerr << "file " << filename << " end not reached!\n";
    // }
    // for (jj = 0; jj < int(times.size()) - 1; ++jj){
    //   if (myiup >= 0 && times[jj] > myiup) break;
    //   sum += 0.5 * (values[jj] + values[jj+1]) * (times[jj+1]  - times[jj]);
    // }
    // sum = (times[0] - times[index]) / (values[index] - values[0]);
    sum = - times[index] / log(values[index]);
    fprintf (fpo, "%f %e\n", kks[ii], sum);
    fclose (fp);
  }
  fclose (fpo);
  
  return 0;
}


