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
  double kk, ll = 0.;
  double begin, end;
  unsigned nblock;
  po::options_description desc ("Allow options");
  
  desc.add_options()
      ("help,h", "print this message")
      ("mode,k", po::value<double > (&kk)->default_value (1.0, "the mode k"))
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
  int nframes = 0;
  double sum = 0.;
  vector<double > record;
  double sum2 = 0.;

  if (! tl.reinit (ifile.c_str())){
    cerr << "cannot init reading " << ifile << endl;
    return 1;
  }
  
  box = tl.getBox ();
  double hz = box[2] * 0.5;

  while (tl.load()){
    if (tl.getTime () < begin) continue;
    if (end != 0 && tl.getTime () > end) continue;
    printf ("# load frame at time %f       \r", tl.getTime());
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
	double tmp0 = (vv[ii][0] - comv[0]) * sin(kk * (xx[ii][2] - hz));
	sum += tmp0;
	record.push_back (tmp0);
	double tmp = 0.;
	tmp += (vv[ii][0] - comv[0]) * (vv[ii][0] - comv[0]);
	tmp += (vv[ii][1] - comv[1]) * (vv[ii][1] - comv[1]);
	tmp += (vv[ii][2] - comv[2]) * (vv[ii][2] - comv[2]);
	sum2 += tmp * sin(kk * (xx[ii][2] - hz));
      }
    }
  }
  printf ("\n");

  sum /= double (posi1.size() * nframes);
  sum *= sum;
  sum2 /= double (posi1.size() * nframes);

  double integral = 0.5 * (hh - ll) - 1./(4. * kk) * sin(2. * kk * (hh - ll));
  integral *= 2;
  cout << "# integral is " << integral << endl;

  BlockAverage ba;
  ba.processData (record, nblock);
  double avg = ba.getAvg();
  double avgError = ba.getAvgError();
  
  // cout << "velo mode^2 is " << sum * integral * box[0] * box[1] << endl;
  printf ("%f  %e %e   %e  %e\n", kk,
	  avg * avg * integral * box[0] * box[1],
	  (avg * avgError + avgError * avgError) * integral * box[0] * box[1],
	  sum * integral * box[0] * box[1],
	  sum2 * integral * box[0] * box[1]);
  
  return 0;
}


