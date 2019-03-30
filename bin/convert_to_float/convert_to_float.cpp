#include <iostream>
using namespace std;
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
    MrcSimple tomo;
    tomo.Read(in_file_name);//You can also use "tomo.Read(cin);" or "cin>>tomo;"
    tomo.PrintStats(cerr);  //Optional (display the tomogram size & format)
    WarnMRCSignedBytes(tomo, in_file_name, cerr);

    // Write the file
    if (argc > 2) {
      string out_file_name(argv[2]);
      cerr << "writing tomogram (in float mode)" << endl;
      tomo.Write(out_file_name); // You can also use "tomo.Write(cout);" 
                                 // or "cout<<tomo;"
    }
  }
  catch (string s) {
    cerr << s << endl; // In case of file format error, print message and exit
    exit(1);
  }
}

