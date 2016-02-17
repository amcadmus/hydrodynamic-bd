#ifndef __LammpsInput_wanghan_h__
#define __LammpsInput_wanghan_h__

#include "StringSplit.h"

using namespace std;

class LammpsInputData 
{
public:
  double xlo, xhi;
  double ylo, yhi;
  double zlo, zhi;
  int natoms;
  vector<double > masses;
  vector<int > type;
  vector<vector<double > > xx;
  void getBlockMasses (ifstream & file);
  void getBlockAtoms  (ifstream & file);
public:
  void read (const string & fname);
}
    ;


#endif
