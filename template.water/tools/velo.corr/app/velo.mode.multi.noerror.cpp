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

namespace po = boost::program_options;

using namespace std;

int main (int argc, char * argv[])
{
  std::string ifile;
  std::string cfile;
  double kkStart, kkStep, kkEnd;
  double ll = 0.;
  double begin, end;
  unsigned nblock;
  po::options_description desc ("Allow options");
  
  desc.add_options()
      ("help,h", "print this message")
      ("mode-start", po::value<double > (&kkStart)->default_value (0.1, "the start of mode k"))
      ("mode-step", po::value<double > (&kkStep)->default_value (0.1, "the step of mode k"))
      ("mode-end", po::value<double > (&kkEnd)->default_value (1.0, "the end of mode k"))
      ("layer,l", po::value<double > (&ll)->default_value (0.0, "the layer width l"))
      ("begin,b", po::value<double > (&begin)->default_value (10.0, "the starting time"))
      ("nblock,n", po::value<unsigned > (&nblock)->default_value (10, "number of block"))
      ("end,e", po::value<double > (&end)->default_value (0.0, "the ending time"))
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

  double minz = 1e10;
  double maxz = -1e10;
  for (unsigned ii = 0; ii < posi1.size(); ++ii){
    if (atomname1[ii] == string("tp03")){
      if (posi1[ii][2] < minz) minz = posi1[ii][2];
      if (posi1[ii][2] > maxz) maxz = posi1[ii][2];
    }
  }
  double hh = (maxz - minz) / 2.;
  printf ("# maxz: %f minz: %f   halfbox should be: %f, real half box: %f   h: %f\n",
	  maxz, minz, (maxz + minz) / 2., box[2] / 2., hh);
  
  TrrLoader tl;
  vector<vector<double > > xx;
  vector<vector<double > > vv;
  vector<vector<double > > ff;
  unsigned nkk = int ((kkEnd - kkStart + kkStep * 0.5) / kkStep) + 1;
  vector<double > kks (nkk);
  for (unsigned ii = 0; ii < nkk; ++ii) {
    kks[ii] = kkStart + ii * kkStep;
    // cout << kks[ii] << endl;
  }
  int nframes = 0;
  // double sum = 0.;
  vector<double > sum(nkk, 0.);
  // double sum2 = 0.;
  

  if (! tl.reinit (ifile.c_str())){
    cerr << "cannot init reading " << ifile << endl;
    return 1;
  }
  
  box = tl.getBox ();
  double hz = box[2] * 0.5;

  while (tl.load()){
    if (tl.getTime () < begin) continue;
    if (end != 0 && tl.getTime () > end) {
      break;
    }
    if (nframes %100 == 0){
      printf ("# load frame at time %f       \r", tl.getTime());
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
    for (unsigned ii = 0; ii < vv.size(); ++ii){
      if (atomname1[ii] == string("tp03")){
	comv[0] += vv[ii][0];
	comv[1] += vv[ii][1];
	comv[2] += vv[ii][2];
      }
    }
    comv[0] /= double(vv.size());
    comv[1] /= double(vv.size());
    comv[2] /= double(vv.size());
    // printf ("frame: %d  comv: %e\n", nframes, comv);
    for (unsigned ii = 0; ii < xx.size(); ++ii){
      if (atomname1[ii] == string("tp03")){
	for (unsigned mm = 0; mm < kks.size(); ++mm){	
	  double tmp0 = (vv[ii][0] - comv[0]) * sin(kks[mm] * (xx[ii][2] - hz));
	  sum[mm] += (tmp0);
	}	
      }
    }
  }
  printf ("\n");

  vector<double > integrals (nkk);
  for (unsigned ii = 0; ii < nkk; ++ii){
    double integral = 0.5 * (hh - ll) - 1./(4. * kks[ii]) * sin(2. * kks[ii] * (hh - ll));
    integrals[ii] = (integral * 2);
    sum[ii] /= double (posi1.size() * nframes);
  }

  for (unsigned ii = 0; ii < nkk; ++ii){  
    printf ("%f  %e    %e\n", kks[ii],
	    sum[ii] * sum[ii] * integrals[ii] * box[0] * box[1],
	    sum[ii]);
  }
  
  return 0;
}


