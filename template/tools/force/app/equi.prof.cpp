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

namespace po = boost::program_options;

using namespace std;

int main (int argc, char * argv[])
{
  std::string ifile;
  std::string ofile;
  std::string cfile;
  int nbins;
  po::options_description desc ("Allow options");
  
  desc.add_options()
      ("help,h", "print this message")
      ("nbins,r", po::value<int > (&nbins)->default_value (30, "number of bins on z direction"))
      ("output,o", po::value<std::string > (&ofile)->default_value (std::string("f.prof.out"), "output file name"))
      ("conf,c",  po::value<std::string > (&cfile)->default_value (std::string("conf.gro"), "input conf file name"))
      ("input,f",  po::value<std::string > (&ifile)->default_value (std::string("traj.trr"), "input traj file name"));

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help") || vm.count("h")){
    std::cout << desc<< "\n";
    return 0;
  }

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
  
  TrrLoader tl;
  vector<vector<double > > xx;
  vector<vector<double > > vv;
  vector<vector<double > > ff;
  vector<double > values0 (nbins, 0.);
  vector<double > values1 (nbins, 0.);
  vector<double > values2 (nbins, 0.);
  vector<int > countSt (nbins, 0);
  int nframes = 0;

  if (! tl.reinit (ifile.c_str())){
    cerr << "cannot init reading " << ifile << endl;
    return 1;
  }
  
  box = tl.getBox ();
  double size = box[2] / nbins;
  double sum00 = 0.;
  double sum01 = 0.;
  double sum10 = 0.;
  double sum11 = 0.;
  double sum20 = 0.;
  double sum21 = 0.;
  // int count0 = 0;
  // int count1 = 0;

  while (tl.load()){
    nframes ++;
    tl.getFrame (xx, vv, ff);
    for (unsigned ii = 0; ii < xx.size(); ++ii){
      int idx = int(xx[ii][2] / size);
      if (idx >= 0 && idx < nbins){
	if (atomname1[ii] == string("tp03")){
	  values0[idx] += ff[ii][0];
	  values1[idx] += ff[ii][1];
	  values2[idx] += ff[ii][2];
	  countSt[idx] ++;
	  if (xx[ii][2] > box[2] * 0.5){
	    sum00 += ff[ii][0];
	    sum10 += ff[ii][1];
	    sum20 += ff[ii][2];
	    // count0 ++;
	  }
	}
	else {
	  if (xx[ii][2] > box[2] * 0.5){
	    sum01 += ff[ii][0];
	    sum11 += ff[ii][1];
	    sum21 += ff[ii][2];
	    // count1 ++;
	  }
	}
      }
      else {
	cerr << "out of range: " << idx
	     << ", z: " << xx[ii][2]
	     << ", box: " << box[2]
	     << ", need pd bd cond" << endl;
      }
    }
  }

  FILE * fp;
  int errno=1;
  if ((fp = fopen(ofile.c_str(), "w")) == NULL) {
    std::cerr << "ERROR: errno=" << errno << " opening file"
	      << " at " << __FILE__ << ":" << __LINE__
	      << std::endl << std::flush;
    exit(1);
  }

  cout << "box z is " << box[2] << endl;
  for (unsigned ii = 0; ii < values0.size(); ++ii){
    if (countSt[ii] != 0){
      values0[ii] /= double(countSt[ii]);
      values1[ii] /= double(countSt[ii]);
      values2[ii] /= double(countSt[ii]);
    }
    fprintf (fp, "%f %e %e %e\n",
	     (ii+0.5) * size,
	     values0[ii],
	     values1[ii],
	     values2[ii]);
  }
  cout << "nframes: " << nframes << endl;
  cout << "sum00: " << sum00 / nframes 
       << " sum01: " << sum01 / nframes
       << " sum: " << (sum00+sum01) / nframes
       << endl;
  cout << "sum10: " << sum10 / nframes 
       << " sum11: " << sum11 / nframes
       << " sum: " << (sum10+sum11) / nframes
       << endl;
  cout << "sum20: " << sum20 / nframes 
       << " sum21: " << sum21 / nframes
       << " sum: " << (sum20+sum21) / nframes
       << endl;

  fclose(fp);
  
  return 0;
}

  
