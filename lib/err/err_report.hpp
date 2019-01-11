#ifndef _ERR_REPORT_HPP   //This check insures that C-preprocessor only includes
#define _ERR_REPORT_HPP   //the following text once (not multiple times)

#include <string>
using namespace std;

class InputErr {
  string msg;
public:
  InputErr(const char *description):msg(description) {}
  InputErr(string description):msg(description) {}
  virtual const char *what() const throw() { return msg.c_str(); }
};

#endif //#ifndef _ERR_REPORT_HPP
