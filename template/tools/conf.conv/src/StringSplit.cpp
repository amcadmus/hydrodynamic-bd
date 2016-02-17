#include "StringSplit.h"

void StringOperation::
split (const std::string & in,
       std::vector<std::string > & out)
{
  std::istringstream iss(in);
  out.clear();
  
  do {
    std::string sub;
    iss >> sub;
    out.push_back (sub);
  // std::vector<std::string > tokens;
  // tokens.push_back (" ");
  // tokens.push_back ("\t");
  // std::copy(std::istream_iterator<std::string>(iss),
  // 	    std::istream_iterator<std::string>(),
  // 	    std::back_inserter<std::vector<std::string> >(tokens));
  } while (iss);

  out.pop_back();
}


