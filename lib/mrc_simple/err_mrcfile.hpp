#ifndef _MRCFILE_ERR_HPP
#define _MRCFILE_ERR_HPP

#include <string>
using namespace std;

class MrcfileErr : public std::exception {
  string msg;
public:
  MrcfileErr(const char *description):msg(description) {}
  MrcfileErr(string description):msg(description) {}
  virtual const char *what() const throw() { return msg.c_str(); }
  virtual ~MrcfileErr() throw (){}
};

#endif //#ifndef _MRCFILE_ERR_HPP
