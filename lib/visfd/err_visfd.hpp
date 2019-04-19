/// @file err_visfd.hpp
/// @brief
/// A simple class used for reporting errors to the caller
/// which are specific to the VISFD library.

#ifndef _VISFD_ERR_HPP
#define _VISFD_ERR_HPP

#include <string>


namespace visfd {


class VisfdErr : public std::exception {
  std::string msg;
public:
  VisfdErr(const char *description):msg(description) {}
  VisfdErr(std::string description):msg(description) {}
  virtual const char *what() const throw() { return msg.c_str(); }
  virtual ~VisfdErr() throw (){}
};


} //namespace visfd


#endif //#ifndef _VISFD_ERR_HPP
