#include <cmath>
#include <string>
#include <iostream>
using namespace std;
#include <err_report.hpp>
#include <mrc_simple.hpp>
#include <threshold.hpp>
#include "settings.hpp"

// (Note: For gcc version 4.8.3, you must compile using: g++ -std=c++11)


int main(int argc, char **argv) {
  try {

    Settings settings; // parse the command-line argument list from the shell
    settings.ParseArgs(argc, argv);

    // Read the input tomogram
    cerr << "Reading tomogram \""<<settings.in_file_name<<"\"" << endl;
    MrcSimple tomo_in;
    tomo_in.Read(settings.in_file_name, false);
    // (Note: You can also use "tomo_in.Read(cin);" or "cin >> tomo;")
    tomo_in.PrintStats(cerr);      //Optional (display the tomogram size & format)

    // ---- mask ----

    // Optional: if there is a "mask", read that too
    MrcSimple mask;
    if (settings.mask_file_name != "") {
      cerr << "Reading mask \""<<settings.mask_file_name<<"\"" << endl;
      mask.Read(settings.mask_file_name, false);
      if ((mask.header.nvoxels[0] != tomo_in.header.nvoxels[0]) ||
          (mask.header.nvoxels[1] != tomo_in.header.nvoxels[1]) ||
          (mask.header.nvoxels[2] != tomo_in.header.nvoxels[2]))
        throw InputErr("Error: The size of the mask does not match the size of the tomogram.\n");
      // The mask should be 1 everywhere we want to consider, and 0 elsewhere.
      if (settings.use_mask_select) {
        for (int iz=0; iz<mask.header.nvoxels[2]; iz++)
          for (int iy=0; iy<mask.header.nvoxels[1]; iy++)
            for (int ix=0; ix<mask.header.nvoxels[0]; ix++)
              if (mask.aaafI[iz][iy][ix] == settings.mask_select)
                mask.aaafI[iz][iy][ix] = 1.0;
              else
                mask.aaafI[iz][iy][ix] = 0.0;
      }
    }

    if (settings.rescale01_in)
      tomo_in.Rescale01(mask.aaafI);



    // ---- Voxel size? ----


    float voxel_width[3] = {1.0, 1.0, 1.0};
    if (settings.voxel_width > 0.0) {
      // Did the user manually specify the width of each voxel?
      voxel_width[0] = settings.voxel_width;
      voxel_width[1] = settings.voxel_width;
      voxel_width[2] = settings.voxel_width;
    }
    else {
      // Otherwise, infer it from the header of the MRC file
      voxel_width[0] = tomo_in.header.cellA[0]/tomo_in.header.nvoxels[0];
      voxel_width[1] = tomo_in.header.cellA[1]/tomo_in.header.nvoxels[1];
      voxel_width[2] = tomo_in.header.cellA[2]/tomo_in.header.nvoxels[2];
      if (settings.voxel_width_divide_by_10) {
        voxel_width[0] *= 0.1;
        voxel_width[1] *= 0.1;
        voxel_width[2] *= 0.1;
      }
      cerr << "voxel width in physical units = ("
           << voxel_width[0] << ", "
           << voxel_width[1] << ", "
           << voxel_width[2] << ")\n";
    }

    if ((voxel_width[0] <= 0.0) ||
        (voxel_width[1] <= 0.0) ||
        (voxel_width[2] <= 0.0))
      throw InputErr("Error in tomogram header: Invalid voxel width(s).\n"
                     "Use the -w argument to specify the voxel width.");

    //if (abs((voxel_width[0] - voxel_width[1])
    //        /
    //        (0.5*(voxel_width[0] + voxel_width[1]))) > 0.0001)
    //  throw InputErr("Error in tomogram header: Unequal voxel widths in the x and y directions.\n"
    //                 "Use the -w argument to specify the voxel width.");

    if ((voxel_width[0] != voxel_width[1]) ||
        (voxel_width[1] != voxel_width[2]))
      throw InputErr("Error in tomogram header: Unequal voxel widths in the x, y and directions.\n"
                     "Use the -w argument to specify the voxel width.");


    float voxel_volume = voxel_width[0] * voxel_width[1] * voxel_width[2];

    float sum_multiplier = 1.0;
    if (settings.multiply_by_voxel_volume)
      sum_multiplier = voxel_volume;

    // ----- thresholding: -----
    
    if (settings.use_thresholds) {

      cerr << "thresholding..." << endl;

      for (int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
        for (int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
          for (int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
            if (! settings.use_dual_thresholds) {
              tomo_in.aaafI[iz][iy][ix] =
                Threshold2(tomo_in.aaafI[iz][iy][ix],
                           settings.in_threshold_01_a,
                           settings.in_threshold_01_b,
                           (settings.in_thresh2_use_clipping
                            ? settings.in_threshold_01_a
                            : settings.in_thresh_a_value),
                           (settings.in_thresh2_use_clipping
                            ? settings.in_threshold_01_b
                            : settings.in_thresh_b_value));
            }
            else
              tomo_in.aaafI[iz][iy][ix] =
                Threshold4(tomo_in.aaafI[iz][iy][ix],
                           settings.in_threshold_01_a,
                           settings.in_threshold_01_b,
                           settings.in_threshold_10_a,
                           settings.in_threshold_10_b,
                           settings.in_thresh_a_value,
                           settings.in_thresh_b_value);
          }
        }
      }
    }

    double sum = 0.0;
    //double denominator = 0.0;
    for (int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
      for (int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
        for (int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
          if (mask.aaafI) {
            // compute a weighted sum using the mask values for weights
            sum += tomo_in.aaafI[iz][iy][ix] * mask.aaafI[iz][iy][ix];
            //denominator += 1.0;
          }
          else
            sum += tomo_in.aaafI[iz][iy][ix];
        }
      }
    }

    //if (mask.aaafI)
    //  sum /= denominator; // divide by the sum of the weights (mask values)

    cout << sum * sum_multiplier << endl;

  } //try {
  catch (InputErr& e) {
    cerr << "\n" << e.what() << endl;
    exit(1);
  }
} // main()

