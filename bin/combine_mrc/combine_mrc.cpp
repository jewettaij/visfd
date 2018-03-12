#include <iostream>
#include <string>
#include <sstream>
#include <vector>
using namespace std;
#include <err_report.h>
#include <mrc_simple.h>
#include <threshold.h>
#include "settings.h"


// (Note: For gcc version 4.8.3, you must compile using: g++ -std=c++11)



int main(int argc, char **argv) {
  try {
    Settings settings; // parse the command-line argument list from the shell
    settings.ParseArgs(argc, argv);

    // Read the input tomograms
    cerr << "Reading tomogram1 \""<<settings.in1_file_name<<"\"" << endl;
    MrcSimple tomo1;
    tomo1.Read(settings.in1_file_name,
	       settings.in_rescale01 && !settings.in1_use_thresholds);
    tomo1.PrintStats(cerr);

    cerr << "\n"
	 << "Reading tomogram2 \""<<settings.in1_file_name<<"\"" << endl;
    MrcSimple tomo2;
    tomo2.Read(settings.in2_file_name,
	       settings.in_rescale01 && !settings.in2_use_thresholds);
    tomo2.PrintStats(cerr);

    if ((tomo1.header.nvoxels[0] != tomo2.header.nvoxels[0]) ||
	(tomo1.header.nvoxels[1] != tomo2.header.nvoxels[1]) ||
	(tomo1.header.nvoxels[2] != tomo2.header.nvoxels[2]))
      throw InputErr("Error: The size of the two input tomograms does not match.\n");

    // Optional: Now apply the threshold filters to each input tomogram
    if (settings.in1_use_thresholds)
      for (int iz=0; iz<tomo1.header.nvoxels[2]; iz++)
	for (int iy=0; iy<tomo1.header.nvoxels[1]; iy++)
	  for (int ix=0; ix<tomo1.header.nvoxels[0]; ix++)
	    tomo1.aaafI[iz][iy][ix] =
	      Threshold4(tomo1.aaafI[iz][iy][ix],
			 settings.in1_threshold_01_a,
			 settings.in1_threshold_01_b,
			 settings.in1_threshold_10_a,
			 settings.in1_threshold_10_b);

    if (settings.in2_use_thresholds)
      for (int iz=0; iz<tomo2.header.nvoxels[2]; iz++)
	for (int iy=0; iy<tomo2.header.nvoxels[1]; iy++)
	  for (int ix=0; ix<tomo2.header.nvoxels[0]; ix++)
	    tomo2.aaafI[iz][iy][ix] =
	      Threshold4(tomo2.aaafI[iz][iy][ix],
			 settings.in2_threshold_01_a,
			 settings.in2_threshold_01_b,
			 settings.in2_threshold_10_a,
			 settings.in2_threshold_10_b);

    // ---- make an array that will store the new tomogram we will create ----

    cerr << "allocating space for new tomogram..." << endl;
    MrcSimple out_tomo = tomo1; //this will take care of allocating the array


    // ---- Now merge the two tomograms ----
    if (settings.operator_char == '*')
      for (int iz=0; iz<tomo1.header.nvoxels[2]; iz++) 
	for (int iy=0; iy<tomo1.header.nvoxels[1]; iy++) 
	  for (int ix=0; ix<tomo1.header.nvoxels[0]; ix++) 
	    out_tomo.aaafI[iz][iy][ix] = 
	      tomo1.aaafI[iz][iy][ix] * tomo2.aaafI[iz][iy][ix];
    else if (settings.operator_char == '+')
      for (int iz=0; iz<tomo1.header.nvoxels[2]; iz++) 
	for (int iy=0; iy<tomo1.header.nvoxels[1]; iy++) 
	  for (int ix=0; ix<tomo1.header.nvoxels[0]; ix++) 
	    out_tomo.aaafI[iz][iy][ix] = 
	      tomo1.aaafI[iz][iy][ix] + tomo2.aaafI[iz][iy][ix];
    else if (settings.operator_char == '/')
      for (int iz=0; iz<tomo1.header.nvoxels[2]; iz++) 
	for (int iy=0; iy<tomo1.header.nvoxels[1]; iy++) 
	  for (int ix=0; ix<tomo1.header.nvoxels[0]; ix++) 
	    out_tomo.aaafI[iz][iy][ix] = 
	      tomo1.aaafI[iz][iy][ix] / tomo2.aaafI[iz][iy][ix];
    else if (settings.operator_char == '-')
      for (int iz=0; iz<tomo1.header.nvoxels[2]; iz++) 
	for (int iy=0; iy<tomo1.header.nvoxels[1]; iy++) 
	  for (int ix=0; ix<tomo1.header.nvoxels[0]; ix++)
	    out_tomo.aaafI[iz][iy][ix] = 
	      tomo1.aaafI[iz][iy][ix] - tomo2.aaafI[iz][iy][ix];
    else
      throw (string("Error: Unrecognized binary operation: \"") +
             string(1,settings.operator_char) + string("\"\n") +
             string("       Must be one of: \"+\", \"*\", \"-\", \"/\"\n"));

    // ---- thresholding and masking: ----
    
    cerr << "thresholding the output tomogram...\n"
	 << "  (Keeping the brightness of the voxels in the final tomogram between 0 1)\n" << endl;

    for (int iz=0; iz<out_tomo.header.nvoxels[2]; iz++)
      for (int iy=0; iy<out_tomo.header.nvoxels[1]; iy++)
        for (int ix=0; ix<out_tomo.header.nvoxels[0]; ix++)
	  if (out_tomo.aaafI[iz][iy][ix] < 0)
	    out_tomo.aaafI[iz][iy][ix] = 0;
	  else if (out_tomo.aaafI[iz][iy][ix] > 1)
	    out_tomo.aaafI[iz][iy][ix] = 1;

    // Write the file
    if (settings.out_file_name != "") {
      cerr << "writing tomogram (in float mode)" << endl;
      out_tomo.Write(settings.out_file_name);
      //(You can also use "out_tomo.Write(cout);" or "cout<<out_tomo;")
    }

  } //try {
  catch (InputErr& e) {
    cerr << "\n" << e.what() << endl;
    exit(1);
  }
} //main()

