#include <iostream>
using namespace std;
#include "mrc_simple.hpp"

// (Note: For gcc version 4.8.3, you must compile using: g++ -std=c++11)



int main(int argc, char **argv) {
  if (argc <= 2) {
    cerr << "\n"
         << "ERROR: This program expects two arguments:  input_file  output_file\n"
         << endl;
    exit(1);
  }

  try {
    // Read the file
    string in_file_name(argv[1]);
    string out_file_name(argv[2]);
    cerr << "Reading tomogram \""<<in_file_name<<"\"" << endl;
    MrcSimple tomo;
    tomo.Read(in_file_name);//You can also use "tomo.Read(cin);" or "cin>>tomo;"
    tomo.PrintStats(cerr);  //Optional (display the tomogram size & format)

    if (tomo.header.use_signed_bytes &&
      tomo.header.mode == MrcHeader::MRC_MODE_BYTE) {
      cerr <<
        "\n"
        "WARNING: File \""<< in_file_name <<"\"\n"
        "         uses signed bytes.\n"
        "         IMOD programs (such as 3dmod) will probably add 128 to each voxel's\n"
        "         brightnesses to insure that they are all positive.\n"
        "         This program does not do that.  If, for some reason, you\n"
        "         want to insure that all brightnesses are positive, use this instead:\n"
        " newstack -mode 2 -in "<<in_file_name<<" -ou "<<out_file_name<<"\n"
                          << endl;
    }

    cerr << "writing tomogram (in float mode)" << endl;
    tomo.Write(out_file_name); // You can also use "tomo.Write(cout);" 
                               // or "cout<<tomo;"
  }
  catch (string s) {
    cerr << s << endl; // In case of file format error, print message and exit
    exit(1);
  }
}

