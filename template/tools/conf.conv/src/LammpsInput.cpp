#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <string>

#include "LammpsInput.h"

#define MAX_LINE_LENGTH 2048

static inline void
die_reach_end (const istream & input,
		    const string & file,
		    const int & line)
{
  if (!input){
    cerr << "cannot getline at file "
	 << file << " : "
	 << line << endl;
    exit (1);
  }
}

static inline void
die_wrong_format (const string & file,
		  const int & line)
{
  cerr << "wrong format error happing at file "
       << file << " : "
       << line << endl;
  exit (1);
}

void LammpsInputData::
getBlockMasses (ifstream & file)
{
  char line[MAX_LINE_LENGTH];
  vector<string > words;  
  die_reach_end (file.getline (line, MAX_LINE_LENGTH), __FILE__, __LINE__);
  for (unsigned ii = 0; ii < masses.size(); ++ii){
    die_reach_end (file.getline (line, MAX_LINE_LENGTH), __FILE__, __LINE__);
    StringOperation::split (line, words);
    if (words.size() < 2) die_wrong_format (__FILE__, __LINE__);
    masses[ii] = atof(words[1].c_str());
  }
  
}

void LammpsInputData::
getBlockAtoms (ifstream & file)
{
  char line[MAX_LINE_LENGTH];
  vector<string > words;  
  die_reach_end (file.getline (line, MAX_LINE_LENGTH), __FILE__, __LINE__);
  for (unsigned ii = 0; ii < type.size(); ++ii){
    die_reach_end (file.getline (line, MAX_LINE_LENGTH), __FILE__, __LINE__);
    StringOperation::split (line, words);
    if (words.size() < 5) die_wrong_format (__FILE__, __LINE__);
    type[ii] = atoi(words[1].c_str());
    xx[ii][0] = atof(words[2].c_str());
    xx[ii][1] = atof(words[3].c_str());
    xx[ii][2] = atof(words[4].c_str());
  }
  
}

void LammpsInputData::
read (const string & fname)
{
  ifstream file (fname.c_str());
  if (! file.is_open()){
    cerr << "cannot open file " << fname << endl;
    exit (1);
  }

  char line[MAX_LINE_LENGTH];
  vector<string > words;

  die_reach_end (file.getline (line, MAX_LINE_LENGTH), __FILE__, __LINE__);
  die_reach_end (file.getline (line, MAX_LINE_LENGTH), __FILE__, __LINE__);
  die_reach_end (file.getline (line, MAX_LINE_LENGTH), __FILE__, __LINE__);
  StringOperation::split (string(line), words);
  if (words.size() < 2) die_wrong_format (__FILE__, __LINE__);
  if ((natoms = atoi(words[0].c_str())) < 0){
      cerr << "number of atoms smaller than 0" << endl;
      exit (2);
  }
  xx.resize (natoms, vector<double > (3, 0.0));
  type.resize (natoms, 0);

  while (! file.eof()){
    file.getline (line, MAX_LINE_LENGTH);
    StringOperation::split (string(line), words);
    cout << line << endl;
    if (words.size () >= 4 &&
	words[2] == string("xlo") && words[3] == string("xhi")){
      xlo = atof(words[0].c_str());
      xhi = atof(words[1].c_str());
    }
    else if (words.size () >= 4 &&
	     words[2] == string("ylo") && words[3] == string("yhi")){
      ylo = atof(words[0].c_str());
      yhi = atof(words[1].c_str());
    }
    else if (words.size () >= 4 &&
	     words[2] == string("zlo") && words[3] == string("zhi")){
      zlo = atof(words[0].c_str());
      zhi = atof(words[1].c_str());
    }
    else if (words.size () >= 3 &&
	     words[1] == string("atom") && words[2] == string("types")){
      masses.resize (atof(words[0].c_str()));
    }
    else if (words.size () >= 1 &&
	     words[0] == "Masses"){
      getBlockMasses (file);
    }
    else if (words.size () >= 1 &&
	     words[0] == "Atoms"){
      getBlockAtoms (file);
    }
  }
}

