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
#include <pthread.h>
#include <omp.h> 

#include "GroFileManager.h"
#include <boost/program_options.hpp>
#include "Trajectory.h"
#include "BlockAverage.h"
#include "AutoCorrelCalculator.h"
#include "ModeRecover.h"

namespace po = boost::program_options;

using namespace std;

int main (int argc, char * argv[])
{
  std::string ifile;
  std::string cfile;
  std::string ofile, omfile;
  double kkStart, kkStep, kkEnd;
  double ll = 0.;
  double begin, end;
  double h_window;
  double cond_numb;
  unsigned nDataInBlock, nDataCorr;
  unsigned every;
  po::options_description desc ("Allow options");
  
  desc.add_options()
      ("help,h", "print this message")
      ("mode-start", po::value<double > (&kkStart)->default_value (0.2, "the start of mode k"))
      ("mode-step", po::value<double > (&kkStep)->default_value (0.01, "the step of mode k"))
      ("mode-end", po::value<double > (&kkEnd)->default_value (3.3, "the end of mode k"))
      ("layer,l", po::value<double > (&ll)->default_value (0.0, "the layer width l"))
      ("window,w", po::value<double > (&h_window)->default_value (-1, "the size of the window, minus value for infinitely large window."))      
      ("cond-numb", po::value<double > (&cond_numb)->default_value (1e6, "the allowed condition number for solving the linear system."))      
      ("begin,b", po::value<double > (&begin)->default_value (5000.0, "the starting time"))
      ("end,e", po::value<double > (&end)->default_value (11000, "the ending time"))
      ("every", po::value<unsigned > (&every)->default_value (1, "skip frame"))
      ("n-data-block,n", po::value<unsigned > (&nDataInBlock)->default_value (100, "number of block"))
      ("n-corr", po::value<unsigned > (&nDataCorr)->default_value (25, "number of corr data"))
      ("conf,c",  po::value<std::string > (&cfile)->default_value (std::string("conf.gro"), "input conf file name"))
      ("input,f",  po::value<std::string > (&ifile)->default_value (std::string("traj.trr"), "input traj file name"))
      ("output,o",  po::value<std::string > (&ofile)->default_value (std::string("inte.corr.out"), "output of integrated corr as a function of modes"))
      ("output-mode,m",  po::value<std::string > (&omfile)->default_value (std::string("corr"), "output of corr, stored by modes"));
  

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help") || vm.count("h")){
    std::cout << desc<< "\n";
    return 0;
  }

  unsigned nkk = int ((kkEnd - kkStart + kkStep * 0.5) / kkStep) + 1;
  vector<double > kks (nkk);
  for (unsigned ii = 0; ii < nkk; ++ii) {
    kks[ii] = kkStart + ii * kkStep;
    // cout << kks[ii] << endl;
  }

  ModeRecover mode_recv (kks, h_window, cond_numb);

  // return 0;

  vector<double > box;
  std::vector<int >  resdindex1;
  std::vector<std::string >   resdname1;
  std::vector<std::string >   atomname1;
  std::vector<int >  atomindex1;
  std::vector<std::vector<double > >  posi1;
  std::vector<std::vector<double > >  velo1;
  GroFileManager::read (cfile,
			resdindex1, resdname1,
			atomname1, atomindex1,
			posi1, velo1, box);

  double minz = 1e10;
  double maxz = -1e10;
  for (unsigned ii = 0; ii < posi1.size(); ++ii){
    if (atomname1[ii] == string("OW")){
      if (posi1[ii][2] < minz) minz = posi1[ii][2];
      if (posi1[ii][2] > maxz) maxz = posi1[ii][2];
    }
  }
  double hh = (maxz - minz) / 2.;
  double hbox = (maxz + minz)/2.; 
  printf ("# maxz: %f minz: %f   halfbox should be: %f, real half box: %f   h: %f\n",
	  maxz, minz, hbox, box[2] / 2., hh);
  
  TrrLoader tl;
  vector<vector<double > > xx;
  vector<vector<double > > vv;
  vector<vector<double > > ff;
  
  int nframes = 0;
  int nCountFrame = 0;
  vector<BlockAverage_acc > bas1 (nkk);
  vector<AutoCorrelCalculator > acc (nkk);
  vector<vector<double > > normal_acc (nkk);
  vector<vector<double > > normal_acc_err (nkk);
  double dt = 0.;
  for (unsigned ii = 0; ii < nkk; ++ii){
    bas1[ii].reinit (nDataInBlock);
    acc[ii].reinit (nDataCorr, nDataInBlock);
  }

  if (! tl.reinit (ifile.c_str())){
    cerr << "cannot init reading " << ifile << endl;
    return 1;
  }
  
  box = tl.getBox ();
  // double hz = box[2] * 0.5;
  
  while (tl.load()){
    if (tl.getTime () < begin) continue;
    if (end != 0 && tl.getTime () > end) {
      break;
    }
    if (nCountFrame == 0){
      dt = tl.getTime();
    }
    else if (nCountFrame == 1){
      dt = tl.getTime() - dt;
    }
    if (nCountFrame++ % every != 0) continue;
    if ((nCountFrame-1) % 10 == 0){
      printf ("# dt: %f load frame at time %f       \r", dt, tl.getTime());
      fflush (stdout);
    }
    
    nframes ++;

    tl.getFrame (xx, vv, ff);
    if (vv.size () != atomname1.size()){
      cerr << "inconsistent vec sizes: "
	   << "vv: " << vv.size()
	   << "conf: " << atomname1.size() << endl;
      return 1;
    }
    double comv[3];
    comv[0] = 0.;
    comv[1] = 0.;
    comv[2] = 0.;
    double comx[3];
    comx[0] = 0.;
    comx[1] = 0.;
    comx[2] = 0.;
    double liquid_x0, liquid_x1;
    liquid_x0 = hbox - h_window;
    liquid_x1 = hbox + h_window;
    double countcomv = 0.;
    double comv0 (0), comv1(0), comv2(0);
    double comx0 (0), comx1(0), comx2(0);
#pragma omp parallel for reduction (+:countcomv, comv0, comv1, comv2, comx0, comx1, comx2) 
    for (unsigned ii = 0; ii < vv.size(); ++ii){
      if (atomname1[ii] == string("OW") ) {    
	xx[ii][2 ]= (xx[ii][2]*15.99940 + xx[ii+1][2]*1.00800 + xx[ii+2][2]*1.00800)/18.0154;
        if( h_window < 0 ||
	    (xx[ii][2] > liquid_x0 && xx[ii][2] < liquid_x1)
	    ) { 
	  vv[ii][0] = (vv[ii][0]*15.99940 + vv[ii+1][0]*1.00800 + vv[ii+2][0]*1.00800)/18.0154;
	  vv[ii][1] = (vv[ii][1]*15.99940 + vv[ii+1][1]*1.00800 + vv[ii+2][1]*1.00800)/18.0154;
	  vv[ii][2] = (vv[ii][2]*15.99940 + vv[ii+1][2]*1.00800 + vv[ii+2][2]*1.00800)/18.0154;
	  xx[ii][0] = (xx[ii][0]*15.99940 + xx[ii+1][0]*1.00800 + xx[ii+2][0]*1.00800)/18.0154;
	  xx[ii][1] = (xx[ii][1]*15.99940 + xx[ii+1][1]*1.00800 + xx[ii+2][1]*1.00800)/18.0154;
	  comv0 +=vv[ii][0];
	  comv1 +=vv[ii][1];
	  comv2 +=vv[ii][2];
	  comx0 +=xx[ii][0];
	  comx1 +=xx[ii][1];
	  comx2 +=xx[ii][2];
	  countcomv += 1.;
        }
      }
    }
    comv[0] = comv0 / double(countcomv);
    comv[1] = comv1 / double(countcomv);
    comv[2] = comv2 / double(countcomv);
    comx[0] = comx0 / double(countcomv);
    comx[1] = comx1 / double(countcomv);
    comx[2] = comx2 / double(countcomv);
    // printf ("comx: %f, hz: %f\n", comx[2], hz);
    // printf ("frame: %d  comv: %e\n", nframes, comv);
    // printf ("countcomv= %f\n ", countcomv );
    vector<double > current_mode (kks.size(), 0.);
#pragma omp parallel for
    for (unsigned mm = 0; mm < kks.size(); ++mm){	
      double count = 0;
      for (unsigned ii = 0; ii < xx.size(); ++ii){
	if (atomname1[ii] == string("OW")){
	  if( h_window < 0 ||
	      (xx[ii][2] > liquid_x0 && xx[ii][2] < liquid_x1)
	      ) {
	    double tmp1 = (vv[ii][0] - comv[0]) * sin(kks[mm] * (xx[ii][2] - comx[2]));
	    current_mode[mm] += (tmp1);
	    count += 1.;
	  }
	}
      }
      current_mode[mm] /= count;
    }
    vector<double > normalized_mode;
    mode_recv.solve (current_mode, normalized_mode);
    for (unsigned mm = 0; mm < kks.size(); ++mm){
      acc[mm].push_back (normalized_mode[mm]);
      bas1[mm].deposite (normalized_mode[mm]);
    }
  }
  printf ("\n");

  vector<double > integrals (nkk);
  for (unsigned ii = 0; ii < nkk; ++ii){
    double integral = 0.5 * (hh - ll) - 1./(4. * kks[ii]) * sin(2. * kks[ii] * (hh - ll));
    integrals[ii] = (integral * 2);
  }
  // double integral = 0.5 * (hh - ll) - 1./(4. * kk) * sin(2. * kk * (hh - ll));
  // integral *= 2;
  // cout << "# integral is " << integral << endl;

  vector<double > inte;
  for (unsigned ii = 0; ii < nkk; ++ii){  
    bas1[ii].calculate ();
    acc[ii].calculate ();
    double sum = 0.;
    for (unsigned jj = 0; jj < acc[ii].nData(); ++jj){
      normal_acc[ii].push_back ((acc[ii].value(jj) - bas1[ii].getAvg() * bas1[ii].getAvg()) / bas1[ii].getVar());
      normal_acc_err[ii].push_back ((acc[ii].error(jj) - bas1[ii].getAvg() * bas1[ii].getAvg()) / bas1[ii].getVar());
    }
    sum += 0.5 * dt * normal_acc[ii].front();
    for (unsigned jj = 1; jj < acc[ii].nData()-1; ++jj){
      sum += dt * normal_acc[ii][jj];
    }
    sum += 0.5 * dt * normal_acc[ii].back();
    inte.push_back (sum);
  }

  FILE * fpo = fopen (ofile.c_str(), "w");
  for (unsigned ii = 0; ii < nkk; ++ii){
    fprintf (fpo, "%f %e\n", kks[ii], inte[ii]);
    char modei[32], modef[32];
    sprintf (modei, "%02d", int(kks[ii]));
    sprintf (modef, "%04d", int(10000 * (kks[ii] - int(kks[ii]) + 0.00005)));
    string filename = omfile + string (modei) + string(".") + string(modef) + string(".out");
    FILE * fp = fopen (filename.c_str(), "w");
    for (unsigned jj = 0; jj < acc[ii].nData(); ++jj){
      // fprintf (fp, "%f %e %e\n", jj * dt, acc[ii].value(jj), acc[ii].error(jj));
      fprintf (fp, "%f %e %e\n", jj * dt, normal_acc[ii][jj], normal_acc_err[ii][jj]);
    }
    fclose (fp);
  }
  fclose (fpo);
  
  return 0;
}


