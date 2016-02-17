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
  std::string ofile;
  double hbin;
  int nbin;
  double begin, end;
  unsigned nDataInBlock;
  unsigned every;
  po::options_description desc ("Allow options");
  
  desc.add_options()
      ("help,h", "print this message")
      ("bin-size,s", po::value<double > (&hbin)->default_value (0.1, "the size of the bin"))
      ("begin,b", po::value<double > (&begin)->default_value (10.0, "the starting time"))
      ("n-data-block,n", po::value<unsigned > (&nDataInBlock)->default_value (10, "number of data in each block"))
      ("end,e", po::value<double > (&end)->default_value (0.0, "the ending time"))
      ("every", po::value<unsigned > (&every)->default_value (1, "skip frame"))
      ("conf,c",  po::value<std::string > (&cfile)->default_value (std::string("conf.gro"), "input conf file name"))
      ("input,f",  po::value<std::string > (&ifile)->default_value (std::string("traj.trr"), "input traj file name"))
      ("output,o",  po::value<std::string > (&ofile)->default_value (std::string("dist.vx2.out"), "output of vx2 distribution"));

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
  int nCountFrame = 0;
  
  if (! tl.reinit (ifile.c_str())){
    cerr << "cannot init reading " << ifile << endl;
    return 1;
  }
  
  box = tl.getBox ();
  double hz = box[2] * 0.5;

  nbin = int ((box[2] + hbin * 0.5) /hbin);
  vector<BlockAverage_acc > avgs (nbin);
  vector<BlockAverage_acc > avgsw (nbin);
  BlockAverage_acc avgvx, avgvy, avgvz;
  avgvx.reinit (nDataInBlock);
  avgvy.reinit (nDataInBlock);
  avgvz.reinit (nDataInBlock);
  for (unsigned ii = 0; ii < avgs.size(); ++ii){
    avgs[ii].reinit (nDataInBlock);
    avgsw[ii].reinit (nDataInBlock);
  }

  double dt;
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
    if ((nCountFrame-1) % 1000 == 0){
      printf ("# dt: %f load frame at time %f       \r", dt, tl.getTime());
      fflush (stdout);
    }
    
    nframes ++;

    tl.getFrame (xx, vv, ff);
    if (xx.size () != atomname1.size()){
      cerr << "inconsistent vec sizes: "
	   << "xx: " << xx.size()
	   << "conf: " << atomname1.size() << endl;
      return 1;
    }

    // int counttp03 = 0;
    // for (unsigned ii = 0; ii < xx.size(); ++ii){
    //   if (atomname1[ii] == string("tp03")){
    // 	counttp03 ++;
    //   }
    // }
    double comx[3];
    comx[0] = 0.;
    comx[1] = 0.;
    comx[2] = 0.;
    double countcomv = 0.;
    for (unsigned ii = 0; ii < vv.size(); ++ii){
      if (atomname1[ii] == string("tp03")){
	comx[0] += xx[ii][0];
	comx[1] += xx[ii][1];
	comx[2] += xx[ii][2];
	countcomv += 1.;
      }
    }
    comx[0] /= double(countcomv);
    comx[1] /= double(countcomv);
    comx[2] /= double(countcomv);
    double shift = hz - comx[2];
    
    for (unsigned ii = 0; ii < xx.size(); ++ii){
      if (atomname1[ii] == string("tp03")){
	int idx = (xx[ii][2] + shift) / hbin;
	if (idx < 0) idx += nbin;
	else if (idx >= nbin) idx -= nbin;
	avgs[idx].deposite (vv[ii][0] * vv[ii][0]);
      }
      // if (atomname1[ii] != string("tp03")){
      else {
	int idx = (xx[ii][2] + shift) / hbin;
	if (idx < 0) idx += nbin;
	else if (idx >= nbin) idx -= nbin;
	avgsw[idx].deposite (vv[ii][0] * vv[ii][0]);
      }
      avgvx.deposite (vv[ii][0] * vv[ii][0]);
      avgvy.deposite (vv[ii][1] * vv[ii][1]);
      avgvz.deposite (vv[ii][2] * vv[ii][2]);
    }
  }
  printf ("\n");

  avgvx.calculate ();
  avgvy.calculate ();
  avgvz.calculate ();
  for (int ii = 0; ii < nbin; ++ii){
    avgs[ii].calculate ();
    avgsw[ii].calculate ();
  }

  FILE * fpo = fopen (ofile.c_str(), "w");
  for (int ii = 0; ii < nbin; ++ii){
    double xx = (0.5 + ii) * hbin - hz;
    fprintf (fpo, "%f %e %e\n", xx, avgs[ii].getAvg(), avgsw[ii].getAvg());
  }
  fclose (fpo);

  printf ("avgvx2: %f  avgvy2: %f  avgvz2: %f\n",
	  avgvx.getAvg(),
	  avgvy.getAvg(),
	  avgvz.getAvg());
  
  return 0;
}


