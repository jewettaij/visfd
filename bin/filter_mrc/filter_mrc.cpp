/// @brief
/// This program allows the user to run one filter operation on an image file.
/// This file contains main() as well as many functions named "Handle...()".
/// Depending on which filter was selected, a different "Handle()" function
/// will be invoked.
/// Each of these "Handle...()" functions will collect the parameters supplied
/// by the user and invoke an corresponding function from the visfd library.


#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <array>
//#include <fftw3.h>  not needed yet
using namespace std;

#ifndef DISABLE_OPENMP
#include <omp.h>       // (OpenMP-specific)
#endif

#include "file_io.hpp"
#include "handlers.hpp"
#include "handlers_unsupported.hpp"


string g_program_name("filter_mrc");
string g_version_string("0.19.4");
string g_date_string("2018-5-26");





int main(int argc, char **argv) {
  cerr << g_program_name << " v" << g_version_string << ", " 
       << g_date_string << "\n" << flush;
  try {
    Settings settings; // parse the command-line argument list from the shell
    settings.ParseArgs(argc, argv);


    #ifndef DISABLE_OPENMP
    #pragma omp parallel
    {
      int rank, nthr;
      rank = omp_get_thread_num();
      //cerr << "rank=" << rank << endl;
      if (rank == 0) {
        nthr = omp_get_num_threads();
        cerr << "  (Using up to "
             << nthr << " threads (cpu cores). You can change this using the \"-np n\"\n"
             << "   argument, or by setting the OMP_NUM_THREADS environment variable.)" << endl;
      }
    }
    #else
    cerr << " (Serial version)" << endl;
    #endif //#ifndef DISABLE_OPENMP


    MrcSimple tomo_in;
    if (settings.in_file_name != "") {
      // Read the input tomogram
      cerr << "Reading tomogram \""<<settings.in_file_name<<"\"" << endl;
      tomo_in.Read(settings.in_file_name, false);
      // (Note: You can also use "tomo_in.Read(cin);" or "cin >> tomo;")
      tomo_in.PrintStats(cerr);      //Optional (display the tomogram size & format)
      WarnMRCSignedBytes(tomo_in, settings.in_file_name, cerr);
    }
    else if ((settings.in_set_image_size[0] > 0) &&
             (settings.in_set_image_size[1] > 0) &&
             (settings.in_set_image_size[2] > 0))
    {
      tomo_in.Resize(settings.in_set_image_size);
    }


    int image_size[3];
    for (int d = 0; d < 3; d++)
      image_size[d] = tomo_in.header.nvoxels[d];

    float voxel_width[3] = {1.0, 1.0, 1.0};

    // ---- Voxel size? ----

    if (settings.voxel_width > 0.0) {
      // Did the user manually specify the width of each voxel?
      voxel_width[0] = settings.voxel_width;
      voxel_width[1] = settings.voxel_width;
      voxel_width[2] = settings.voxel_width;
    }
    else {
      // Otherwise, infer it from the header of the MRC file
      voxel_width[0] = tomo_in.header.cellA[0]/image_size[0];
      voxel_width[1] = tomo_in.header.cellA[1]/image_size[1];
      voxel_width[2] = tomo_in.header.cellA[2]/image_size[2];
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

    if ((abs((voxel_width[0] - voxel_width[1]) /
             (0.5*(voxel_width[0] + voxel_width[1]))) > 0.0001) ||
        (abs((voxel_width[0] - voxel_width[2]) /
             (0.5*(voxel_width[0] + voxel_width[2]))) > 0.0001))
    {
      stringstream err_msg;
      err_msg << "Error in tomogram header: Unequal voxel widths in the x, y and z directions:\n"
        "  voxel_width_x = " << voxel_width[0] << "\n"
        "  voxel_width_y = " << voxel_width[1] << "\n"
        "  voxel_width_z = " << voxel_width[2] << "\n"
        "Tomograms with non-cubic voxels are not supported.\n"
        "Use interpolation to stretch the image in the shorter\n"
        "directions so that all 3 voxel widths agree.\n"
        "Alternatively, use the \"-w WIDTH\" argument to specify the same voxel width for\n"
        "all 3 directions.  (For example \"-w "<<voxel_width[0]<<"\")\n";
      throw VisfdErr(err_msg.str().c_str());
    }

    settings.surface_tv_sigma /= voxel_width[0];
    for (int d=0; d<3; d++) {
      settings.width_a[d] /= voxel_width[d];
      settings.width_b[d] /= voxel_width[d];
      settings.dogsf_width[d] /= voxel_width[d];
      #ifndef DISABLE_BOOTSTRAPPING
      settings.f_bs_scramble_radius[d] /= voxel_width[d];
      settings.bs_scramble_radius[d] = ceil(settings.bs_scramble_radius[d]);
      #endif
      #ifndef DISABLE_TEMPLATE_MATCHING
      settings.template_background_radius[d] /= voxel_width[d];
      #endif
      //Commenting out:  The user cannot set filter_truncate_halfwidth directly anymore:
      //settings.filter_truncate_halfwidth[d] = floor(settings.filter_truncate_ratio * 
      //                                              MAX(settings.width_a[d],
      //                                                  settings.width_b[d]));
    }

    // At some point I was trying to be as general as possible and allowed
    // for the possibility that voxels need not be cubes (same width x,y,z)
    // Now, I realize that allowing for this possibility would slow down some
    // calculations considerably, so I just assume cube-shaped voxels:
    //assert((voxel_width[0] == voxel_width[1]) &&
    //       (voxel_width[1] == voxel_width[2]));
    for (size_t ir = 0; ir < settings.blob_diameters.size(); ir++)
      settings.blob_diameters[ir] /= voxel_width[0];

    settings.sphere_decals_diameter /= voxel_width[0];
    if (! settings.sphere_decals_shell_thickness_is_ratio)
      settings.sphere_decals_shell_thickness /= voxel_width[0];

    // now use the voxel_width (distance-to-voxel converter)
    // to read in coordinates from various files:
    if (! settings.is_training_data_pos_in_voxels)
      for (size_t i = 0; i < settings.training_data_pos_crds.size(); i++)
        for (int d = 0; d < 3; d++)
          settings.training_data_pos_crds[i][d] /= voxel_width[d];

    if (! settings.is_training_data_neg_in_voxels)
      for (size_t i = 0; i < settings.training_data_neg_crds.size(); i++)
        for (int d = 0; d < 3; d++)
          settings.training_data_neg_crds[i][d] /= voxel_width[d];

    if (! settings.is_must_link_constraints_in_voxels)
      for (size_t i = 0; i < settings.must_link_constraints.size(); i++)
        for (size_t j = 0; j < settings.must_link_constraints[i].size(); j++)
          for (int d = 0; d < 3; d++)
            settings.must_link_constraints[i][j][d] /= voxel_width[d];


    // ---- mask ----

    // Optional: if there is a "mask", read that too
    MrcSimple mask;
    if (settings.mask_file_name != "") {
      cerr << "Reading mask \""<<settings.mask_file_name<<"\"" << endl;
      mask.Read(settings.mask_file_name, false);
      WarnMRCSignedBytes(mask, settings.mask_file_name, cerr);
      if ((mask.header.nvoxels[0] != image_size[0]) ||
          (mask.header.nvoxels[1] != image_size[1]) ||
          (mask.header.nvoxels[2] != image_size[2])) {
        if ((image_size[0]!=0) && (image_size[1]!=0) && (image_size[2]!=0)) {
          image_size[0] = mask.header.nvoxels[0];
          image_size[1] = mask.header.nvoxels[1];
          image_size[2] = mask.header.nvoxels[2];
          tomo_in.Resize(image_size);
        }
        else
          throw VisfdErr("Error: The size of the mask image does not match the size of the input image.\n");
      }
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

    else if (settings.mask_rectangle_xmin <= settings.mask_rectangle_xmax) {

      mask = tomo_in; //allocate an array for an image the same size as tomo_in
      if (! settings.is_mask_rectangle_in_voxels) {
        settings.mask_rectangle_xmin /= voxel_width[0];
        settings.mask_rectangle_xmax /= voxel_width[0];
        settings.mask_rectangle_ymin /= voxel_width[1];
        settings.mask_rectangle_ymax /= voxel_width[1];
        settings.mask_rectangle_zmin /= voxel_width[2];
        settings.mask_rectangle_zmax /= voxel_width[2];
      }

      for (int iz=0; iz<mask.header.nvoxels[2]; iz++) {
        for (int iy=0; iy<mask.header.nvoxels[1]; iy++) {
          for (int ix=0; ix<mask.header.nvoxels[0]; ix++) {
            if ((floor(settings.mask_rectangle_xmin) <= ix) &&
                (ix <= ceil(settings.mask_rectangle_xmax)) && 
                (floor(settings.mask_rectangle_ymin) <= iy) &&
                (iy <= ceil(settings.mask_rectangle_ymax)) && 
                (floor(settings.mask_rectangle_zmin) <= iz) &&
                (iz <= ceil(settings.mask_rectangle_zmax)))
              mask.aaafI[iz][iy][ix] = 1.0;
            else
              mask.aaafI[iz][iy][ix] = 0.0;
          }
        }
      }

    } //else if (settings.mask_rectangle_xmin <= settings.mask_rectangle_xmax) {

    if (settings.rescale_min_max_in) {
      tomo_in.Rescale01(mask.aaafI,
                        settings.in_rescale_min,
                        settings.in_rescale_max);
      tomo_in.FindMinMaxMean();
    }

    // ---- make an array that will store the new tomogram we will create ----

    cerr << "allocating space for new tomogram..." << endl;
    MrcSimple tomo_out = tomo_in; //this will take care of allocating the array

    if ((voxel_width[0] <= 0.0) ||
        (voxel_width[1] <= 0.0) ||
        (voxel_width[2] <= 0.0))
      throw VisfdErr("Error in tomogram header: Invalid voxel width(s).\n"
                     "Use the -w argument to specify the voxel width.");


    // ---- filtering ----
    // perform an operation which generates a new image based on the old image


    if (settings.filter_type == Settings::NONE) {
      cerr << "filter_type = Intensity Map <No convolution filter specified>\n";
      // Not needed:
      //for (int iz = 0; iz < size[2]; iz++)
      //  for (int iy = 0; iy < size[1]; iy++)
      //    for (int ix = 0; ix < size[0]; ix++)
      //      tomo_out.aaafI[iz][iy][ix]=tomo_in.aaafI[iz][iy][ix];
      // (We have copied the contents from tomo_in into tomo_out already.)
    } 



    else if (settings.filter_type == Settings::GAUSS) {

      HandleGauss(settings, tomo_in, tomo_out, mask, voxel_width);
        
    } // if (settings.filter_type == Settings::GAUSS)



    else if (settings.filter_type == Settings::GGAUSS) {

      HandleGGauss(settings, tomo_in, tomo_out, mask, voxel_width);

    } //else if (settings.filter_type == Settings::GGAUSS)



    else if (settings.filter_type == Settings::DOG) {

      HandleDog(settings, tomo_in, tomo_out, mask, voxel_width);

    } //if (settings.filter_type == Settings::DOG)



    else if (settings.filter_type == Settings::DOGG) {

      HandleDogg(settings, tomo_in, tomo_out, mask, voxel_width);

    } //if (settings.filter_type == Settings::DOGG)



    #ifndef DISABLE_DOGGXY
    else if (settings.filter_type == Settings::DOGGXY) {

      HandleDoggXY(settings, tomo_in, tomo_out, mask, voxel_width);

    } //else if (settings.filter_type = Settings::DOGGXY)
    #endif


    else if (settings.filter_type == Settings::DOG_SCALE_FREE) {

      HandleDogScaleFree(settings, tomo_in, tomo_out, mask, voxel_width);

    } //if (settings.filter_type == Settings::DOG_SCALE_FREE)



    else if (settings.filter_type == Settings::BLOB) {

      HandleBlobDetector(settings, tomo_in, tomo_out, mask, voxel_width);

    } //if (settings.filter_type == Settings::DOG)



    else if (settings.filter_type == Settings::RIDGE_SURFACE) {

      // find surface ridges (ie membranes or wide tubes)
      HandleRidgeDetector(settings, tomo_in, tomo_out, mask, voxel_width);

    }


    else if (settings.filter_type == Settings::WATERSHED) {

      // perform watershed segmentation
      HandleWatershed(settings, tomo_in, tomo_out, mask, voxel_width);

    }


    else if (settings.filter_type == Settings::CLUSTER_CONNECTED) {

      // cluster adjacent nearby voxels into disconnected "islands"
      // (this is similar to watershed segmentation)
      HandleClusterConnected(settings, tomo_in, tomo_out, mask, voxel_width);

    }


    else if (settings.filter_type == Settings::LOCAL_FLUCTUATIONS) {

      HandleLocalFluctuations(settings, tomo_in, tomo_out, mask, voxel_width);

    } //else if (settings.filter_type == Settings::TEMPLATE_GGAUSS)




    // ----- template matching with error reporting (probably not useful) -----


    else if (settings.filter_type == Settings::TEMPLATE_GGAUSS) {

      #ifndef DISABLE_TEMPLATE_MATCHING
      HandleTemplateGGauss(settings, tomo_in, tomo_out, mask, voxel_width);
      #endif //#ifndef DISABLE_TEMPLATE_MATCHING

    } //else if (settings.filter_type == Settings::TEMPLATE_GGAUSS)



    else if (settings.filter_type == Settings::TEMPLATE_GAUSS) {

      #ifndef DISABLE_TEMPLATE_MATCHING
      HandleTemplateGauss(settings, tomo_in, tomo_out, mask, voxel_width);
      #endif //#ifndef DISABLE_TEMPLATE_MATCHING

    } //else if (settings.filter_type == Settings::TEMPLATE_GAUSS)


    else if (settings.filter_type == Settings::BLOB_RADIAL_INTENSITY) {

      #ifndef DISABLE_INTENSITY_PROFILES
      HandleBlobRadialIntensity(settings, tomo_in, tomo_out, mask, voxel_width);
      #endif //#ifndef DISABLE_TEMPLATE_MATCHING

    } //else if (settings.filter_type == Settings::TEMPLATE_GAUSS)




    else if (settings.filter_type == Settings::BOOTSTRAP_DOGG) {

      #ifdef DISABLE_BOOSTRAPPING
      HandleBootstrappDogg(settings, tomo_in, tomo_out, mask, voxel_width);
      #endif //#ifndef DISABLE_BOOTSTRAPPING

    } //else if (settings.filter_type == Settings::BOOTSTRAP_DOGG) {


    // ----- distance_filter -----

    else if (settings.filter_type == settings.DISTANCE_TO_POINTS) {

      HandleDistanceToPoints(settings, tomo_in, tomo_out, mask, voxel_width);

    }

    // ----- sphere_decals_filter -----

    else if (settings.filter_type == settings.SPHERE_DECALS) {

      HandleDrawSpheres(settings, tomo_in, tomo_out, mask, voxel_width);

    }

    else if (settings.filter_type == settings.SPHERE_NONMAX_SUPPRESSION) {

      vector<array<float,3> > crds;
      vector<float> diameters;
      vector<float> scores;

      HandleBlobsNonmaxSuppression(settings,
                                   mask,
                                   voxel_width,
                                   crds,
                                   diameters,
                                   scores);

    }

    else {
      assert(false);  //should be one of the choices above
    }





    // --- Exchange light voxels for dark voxels ? ---

    if (settings.invert_output)
      tomo_out.Invert(mask.aaafI);




    // ----- thresholding: -----
    
    if (settings.use_intensity_map) {

      HandleThresholds(settings, tomo_in, tomo_out, mask, voxel_width);

    }

    // ----- find local minima or maxima -----

    if (settings.find_minima || settings.find_maxima) {

      HandleExtrema(settings, tomo_in, tomo_out, mask, voxel_width);

    }





    // ----- masking: -----
    // Also, after thresholding, check again to see if the mask is zero
    // at this location.  If so, make sure that aaafI is 0 there
    // (or whatever the user has requested us to put there).
    // Do this after thresholding and inverting.
    if ((mask.aaafI) && (settings.use_mask_out))
      for (int iz=0; iz<mask.header.nvoxels[2]; iz++)
        for (int iy=0; iy<mask.header.nvoxels[1]; iy++)
          for (int ix=0; ix<mask.header.nvoxels[0]; ix++)
            if (mask.aaafI[iz][iy][ix] == 0.0)
              tomo_out.aaafI[iz][iy][ix] = settings.mask_out;

    // --- Rescale so that the lowest, highest voxels have density 0 and 1? ---

    if (settings.rescale_min_max_out) {
      tomo_out.Rescale01(mask.aaafI,
                         settings.out_rescale_min,
                         settings.out_rescale_max);
    }

    tomo_out.FindMinMaxMean();


    // ------ Write the file ------
    if (settings.out_file_name != "") {
      cerr << "writing tomogram (in 32-bit float mode)" << endl;
      tomo_out.Write(settings.out_file_name);
      // (You can also use "file_stream << tomo_out;")
    }

  } // try {

  catch (const std::exception& e) {
    cerr << "\n" << e.what() << endl;
    exit(1);
  }
} //main()


