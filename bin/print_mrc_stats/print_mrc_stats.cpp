#include <iostream>
#include <fstream>
using namespace std;
#include "err.hpp"
#include "mrc_simple.hpp"

// (Note: For gcc version 4.8.3, you must compile using: g++ -std=c++11)



int main(int argc, char **argv) {
  if (argc <= 1) {
    cerr << "This program expects at an argument: input_file [output_file]\n"
         << "No tomogram file specified.  Exiting.\n" << endl;
    exit(0);
  }

  try {
    // Read the file
    string in_file_name(argv[1]);
    cerr << "Reading tomogram \""<<in_file_name<<"\"" << endl;
    fstream mrc_file;
    mrc_file.open(in_file_name.c_str(), ios::binary | ios::in);
    if (! mrc_file) 
      throw InputErr("Error: unable to open \"" + in_file_name + "\" for reading.\n");
    MrcHeader header;
    header.Read(mrc_file);
    header.PrintStats(cout);
    mrc_file.close();
  }
  catch (const std::exception& e) {
    cerr << "\n" << e.what() << endl;
    exit(1);
  }
}

