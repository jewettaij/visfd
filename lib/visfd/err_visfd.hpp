#ifndef _VISFD_ERR_HPP
#define _VISFD_ERR_HPP

#include <string>
using namespace std;

class VisfdErr : public std::exception {
  string msg;
public:
  VisfdErr(const char *description):msg(description) {}
  VisfdErr(string description):msg(description) {}
  virtual const char *what() const throw() { return msg.c_str(); }
  virtual ~VisfdErr() throw (){}
};

#endif //#ifndef _VISFD_ERR_HPP
