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
  unsigned nDataInBlock;
  unsigned every;
  po::options_description desc ("Allow options");
  
  desc.add_options()
    ("help,h", "print this message")
    ("mode-start", po::value<double > (&kkStart)->default_value (0.1, "the start of mode k"))
    ("mode-step", po::value<double > (&kkStep)->default_value (0.1, "the step of mode k"))
    ("mode-end", po::value<double > (&kkEnd)->default_value (1.0, "the end of mode k"))
    ("layer,l", po::value<double > (&ll)->default_value (0.0, "the layer width l"))
    ("begin,b", po::value<double > (&begin)->default_value (10.0, "the starting time"))
    ("n-data-block,n", po::value<unsigned > (&nDataInBlock)->default_value (10, "number of block"))
    ("end,e", po::value<double > (&end)->default_value (0.0, "the ending time"))
    ("every", po::value<unsigned > (&every)->default_value (1, "skip frame"))
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
  int nCountFrame = 0;
  // double sum = 0.;
  // double sum2 = 0.;
  // vector<vector<double > > record(nkk);
  vector<BlockAverage_acc > bas1 (nkk);
  vector<BlockAverage_acc > bas2 (nkk);
  vector<BlockAverage_acc > bas3 (nkk);
  vector<BlockAverage_acc > bas4 (nkk);
  for (unsigned ii = 0; ii < nkk; ++ii){
    bas1[ii].reinit (nDataInBlock);
    bas2[ii].reinit (nDataInBlock);
    bas3[ii].reinit (nDataInBlock);
    bas4[ii].reinit (nDataInBlock);
  }

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
    if (nCountFrame++ % every != 0) continue;
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
	for (unsigned mm = 0; mm < kks.size(); ++mm){	
	  double tmp1 = (vv[ii][0] - comv[0]) * sin(kks[mm] * (xx[ii][2] - hz));
	  // sum += tmp0;
	  // record[mm].push_back (tmp0);
	  bas1[mm].deposite (tmp1);
	  double tmp = 0.;
	  tmp += (vv[ii][0] - comv[0]) * (vv[ii][0] - comv[0]);
	  tmp += (vv[ii][1] - comv[1]) * (vv[ii][1] - comv[1]);
	  tmp += (vv[ii][2] - comv[2]) * (vv[ii][2] - comv[2]);
	  // sum2 += tmp * sin(kk * (xx[ii][2] - hz));
	  double tmp2 = (tmp) * sin(kks[mm] * (xx[ii][2] - hz));
	  bas2[mm].deposite (tmp2);
	  double tmp4 = sqrt(tmp) * fabs(sin(kks[mm] * (xx[ii][2] - hz)));
	  bas4[mm].deposite (tmp4);
	}
      }
    }
    for (unsigned mm = 0; mm < kks.size(); ++mm){	
      double tmpsum = 0.;
      int tmpcount = 0;
      for (unsigned ii = 0; ii < xx.size(); ++ii){
	if (atomname1[ii] == string("tp03")){
	  tmpsum += (vv[ii][0] - comv[0]) * sin(kks[mm] * (xx[ii][2] - hz));
	  tmpcount ++;
	}
      }
      tmpsum /= double (tmpcount);
      bas3[mm].deposite (fabs(tmpsum));
    }
  }
  printf ("\n");

  // sum /= double (posi1.size() * nframes);
  // sum *= sum;
  // sum2 /= double (posi1.size() * nframes);

  vector<double > integrals (nkk);
  for (unsigned ii = 0; ii < nkk; ++ii){
    double integral = 0.5 * (hh - ll) - 1./(4. * kks[ii]) * sin(2. * kks[ii] * (hh - ll));
    integrals[ii] = (integral * 2);
  }
  // double integral = 0.5 * (hh - ll) - 1./(4. * kk) * sin(2. * kk * (hh - ll));
  // integral *= 2;
  // cout << "# integral is " << integral << endl;

  vector<double > avgs (nkk);
  vector<double > avgErrors (nkk);
  for (unsigned ii = 0; ii < nkk; ++ii){  
    bas1[ii].calculate ();
    bas2[ii].calculate ();
    bas3[ii].calculate ();
    bas4[ii].calculate ();
    avgs[ii] = (bas1[ii].getAvg());
    avgErrors[ii] = (bas1[ii].getAvgError());
  }
  
  // cout << "velo mode^2 is " << sum * integral * box[0] * box[1] << endl;
  printf ("# kk		energy (error)		mode (error)		mode2 (error)		<fabs(mode)> (error)	 <|v|.|sin(kz)|> (error)\n");
  for (unsigned ii = 0; ii < nkk; ++ii){  
    printf ("%f  %e %e    %e %e    %e %e    %e %e    %e %e\n",
	    kks[ii],
	    avgs[ii] * avgs[ii] * integrals[ii] * box[0] * box[1],
	    (avgs[ii] * avgErrors[ii] + avgErrors[ii] * avgErrors[ii]) * integrals[ii] * box[0] * box[1],
	    avgs[ii], avgErrors[ii],
	    bas2[ii].getAvg(), bas2[ii].getAvgError(),
	    bas3[ii].getAvg(), bas3[ii].getAvgError(),
	    bas4[ii].getAvg(), bas4[ii].getAvgError()
	    );
  }
  
  return 0;
}


