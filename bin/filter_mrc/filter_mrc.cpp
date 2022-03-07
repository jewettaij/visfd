/// @brief
/// This program allows the user to run one filter operation on an image file.
/// This file contains main() as well as many functions named "Handle...()".
/// After reading the image files, depending on which filter was selected,
/// a different "Handle()" function will be invoked.
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
string g_version_string("0.29.21");
string g_date_string("2022-3-06");




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
      // By default, set the voxel_width to 1.
      // The voxel_width = ratio of cellA / nvoxels.  So we make sure it's 1:
      tomo_in.header.cellA[0] = tomo_in.header.nvoxels[0];
      tomo_in.header.cellA[1] = tomo_in.header.nvoxels[1];
      tomo_in.header.cellA[2] = tomo_in.header.nvoxels[2];
    }


    // ---- mask ----

    // Optional: if there is a "mask", read that too
    MrcSimple mask;
    if (settings.mask_file_name != "") {
      cerr << "Reading mask \""<<settings.mask_file_name<<"\"" << endl;
      mask.Read(settings.mask_file_name, false);
      WarnMRCSignedBytes(mask, settings.mask_file_name, cerr);
      if ((mask.header.nvoxels[0] != tomo_in.header.nvoxels[0]) ||
          (mask.header.nvoxels[1] != tomo_in.header.nvoxels[1]) ||
          (mask.header.nvoxels[2] != tomo_in.header.nvoxels[2])) {
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


    // ---- How many voxels are in the image we will be processing? ----
    int image_size[3]; // <-- store here
    for (int d = 0; d < 3; d++)
      image_size[d] = tomo_in.header.nvoxels[d];

    // ---- Voxel width? ----
    float voxel_width[3] = {-1.0 , -1.0, -1.0}; //initial impossible value
    DetermineVoxelWidth(settings, tomo_in, mask, voxel_width);


    // ---- Reduce the size of the input image (using binning)? ----

    for (int d=0; d < 3; d++) {
      // To be safe, back up the original image size and voxel width first
      settings.image_size_orig[d] = tomo_in.header.nvoxels[d]; //image size
      settings.cellA_orig[d] = tomo_in.header.cellA[d]; //voxel width
    }
    // Did the user explicitly ask us to resample the image using binning?
    if (settings.resize_with_binning > 1) {
      HandleBinning(settings, tomo_in, mask, voxel_width);
      // Since we reduced the resolution, we should increase voxel width
      DetermineVoxelWidth(settings, tomo_in, mask, voxel_width);
      // Update the "image_size" array (just a copy of tomo_in.header.nvoxels[])
      for (int d = 0; d < 3; d++)
        image_size[d] = tomo_in.header.nvoxels[d];
    } //if (settings.resize_with_binning > 0)

    // If the user did not ask us to bin the image, perhaps we should anyway?
    else if (settings.resize_with_binning == 0) {
      settings.resize_with_binning = 1;
      if (settings.tv_sigma > 0) {
        float feature_detection_sigma = settings.width_a[0];
        float voting_distance = settings.tv_sigma;
        if (feature_detection_sigma > 1.8*voxel_width[0])
        {
          settings.resize_with_binning = int(ceil(feature_detection_sigma /
                                                  (1.8*voxel_width[0])));
          cerr <<
            "------------------------------------------------------------------\n"
            "------------------------------------------------------------------\n"
            "------------------------------------------------------------------\n"
            "---  WARNING: Tensor-voting requested with a feature width of\n"
            "---           sigma = " << feature_detection_sigma/voxel_width[0]
               << ",\n"
            "---           and a voting distance of\n"
            "--- voting_distance = " << voting_distance/voxel_width[0]
               << "  IN VOXELS\n"
            "---           This will be very slow unless binning is used.\n"
            "--- BINNING THE IMAGE BY A FACTOR OF " <<settings.resize_with_binning
               << "\n"
            "---           This may reduce the precision of the detected features.\n"
            "---           To prevent this, use the \"-bin 1\" argument.\n"
            "------------------------------------------------------------------\n"
            "------------------------------------------------------------------\n"
            "------------------------------------------------------------------"
               << endl;
          // perhaps later, I'll change the criteria to something like this:
          // if (voting_distance > 10.0*voxel_width[0]) ...
          HandleBinning(settings, tomo_in, mask, voxel_width);
          // Since we reduced the resolution, we should recalculate voxel width
          DetermineVoxelWidth(settings, tomo_in, mask, voxel_width);
          // Update the "image_size" array
          for (int d = 0; d < 3; d++)
            image_size[d] = tomo_in.header.nvoxels[d];
        }
      } //if (settings.tv_sigma > 0)
      else if (settings.blob_diameters.size() > 0)
      {
        if (settings.blob_diameters[0] > 15.0*voxel_width[0])
        {
          settings.resize_with_binning = int(ceil(settings.blob_diameters[0] /
                                                  (15.0*voxel_width[0])));
          cerr <<
            "------------------------------------------------------------------\n"
            "------------------------------------------------------------------\n"
            "------------------------------------------------------------------\n"
            "---  WARNING: Blob detection requested with a minimum sigma of\n"
            "---           sigma = " <<settings.blob_diameters[0]/voxel_width[0]
               << "   IN VOXELS\n"
            "---           This will be very slow unless binning is used.\n"
            "--- BINNING THE IMAGE BY A FACTOR OF " <<settings.resize_with_binning
               << "\n"
            "---           This may reduce the precision of the detected blobs.\n"
            "---           To prevent this, use the \"-bin 1\" argument.\n"
            "------------------------------------------------------------------\n"
            "------------------------------------------------------------------\n"
            "------------------------------------------------------------------"
               << endl;
          // perhaps later, I'll change the criteria to something like this:
          // if (voting_distance > 10.0*voxel_width[0]) ...
          HandleBinning(settings, tomo_in, mask, voxel_width);
          // Since we reduced the resolution, we should recalculate voxel width
          DetermineVoxelWidth(settings, tomo_in, mask, voxel_width);
          // Update the "image_size" array
          for (int d = 0; d < 3; d++)
            image_size[d] = tomo_in.header.nvoxels[d];
        } //if (settings.blob_diameters[0] > 10.0*voxel_width[0])
      } //else if (settings.blob_diameters.size() > 0)
    } //else if (settings.resize_with_binning == 0)




    // --- Rescale parameter coordinates according to voxel_width ---


    // Now that we know voxel_width, rescale coordinates
    // of the regions that define the image mask (if applicable).
    if (settings.vMaskRegions.size() > 0) {

      // first, make sure the mask is allocated and initialized
      if (mask.aaafI == nullptr) {
        mask=tomo_in; //allocate an array for an image the same size as tomo_in
        for (int iz=0; iz<mask.header.nvoxels[2]; iz++)
          for (int iy=0; iy<mask.header.nvoxels[1]; iy++)
            for (int ix=0; ix<mask.header.nvoxels[0]; ix++)
              mask.aaafI[iz][iy][ix] = 0.0; //ignore each voxel by default
      }

      // now, we loop through the vMaskRegions array, rescaling their coords
      for (int i = 0; i < settings.vMaskRegions.size(); i++) {
        //determine how to scale the coordinates
        float scale = 1.0;
        if (settings.is_mask_crds_in_voxels) {
          // Did we use binning?  If so, divide coordinates by the bin size.
          scale = 1.0 / settings.resize_with_binning;
          // Note that we only do this if we are ignoring the voxel_width,
          // because the code that takes care of binning also updates
          // the voxel_width automatically.  So if we are dividing by the
          // voxel width, we can igore settings.resize_with_binning.
        }
        else
          // Divide coordinates by the voxel width.
          // assume voxel_width[0] == voxel_width[1] == voxel_width[2]
          // (I haven't yet implemented non-cubical bin-widths.-A 2021-7-04)
          scale = 1.0 / voxel_width[0];

        switch (settings.vMaskRegions[i].type) {
        case SimpleRegion<float>::RECT:
          {
            settings.vMaskRegions[i].data.rect.xmin *= scale;
            settings.vMaskRegions[i].data.rect.xmax *= scale;
            settings.vMaskRegions[i].data.rect.ymin *= scale;
            settings.vMaskRegions[i].data.rect.ymax *= scale;
            settings.vMaskRegions[i].data.rect.zmin *= scale;
            settings.vMaskRegions[i].data.rect.zmax *= scale;
          }
          break;
        case SimpleRegion<float>::SPHERE:
          {
            settings.vMaskRegions[i].data.sphere.r  *= scale;
            settings.vMaskRegions[i].data.sphere.x0 *= scale;
            settings.vMaskRegions[i].data.sphere.y0 *= scale;
            settings.vMaskRegions[i].data.sphere.z0 *= scale;
          }
          break;
        default:
          assert(false); //this line should not be reached
          break;
        } //switch (settings.vMaskRegions[i].type)
      } //for (int i = 0; i < settings.vMaskRegions.size(); i++)
      // Now use VISFD's "DrawRegions()" function to fill the
      // voxels in these regions of the mask with brightness=1,
      // which means we want to include these voxels in the mask.
      // The remaining voxels will have brightness=0,
      // (unless they were read from a file beforehand)
      // ...and will be excluded from the mask.

      DrawRegions(image_size,
                  mask.aaafI,
                  static_cast<const float* const* const*>(nullptr),
                  settings.vMaskRegions,
                  true);//<--allows us to add and subtract regions from the mask

    } //if (settings.vMaskRegions.size() > 0)


    // Now that we know the voxel_width, rescale morphological filter widths
    // At some point I was trying to be as general as possible and allowed
    // for the possibility that voxels need not be cubes (same width x,y,z)
    // Now, I realize that allowing for this possibility would slow down some
    // calculations considerably, so I just assume cube-shaped voxels:
    //assert((voxel_width[0] == voxel_width[1]) &&
    //       (voxel_width[1] == voxel_width[2]));

    settings.morphology_r /= voxel_width[0];
    settings.morphology_rmax /= voxel_width[0];


    // Note: If it is a positive number, max_distance_to_feature is specified
    // by the user in units of voxels (not physical distance).
    // But it is a negative number, we need to convert the distance units.
    if (settings.max_distance_to_feature < 0.0)
      settings.max_distance_to_feature /= -voxel_width[0];
    else
      settings.max_distance_to_feature /= settings.resize_with_binning;

    // Now that we know the voxel_width, rescale any tensor-voting
    // parameters that have units of physical distance.
    settings.tv_sigma /= voxel_width[0];
    for (int d=0; d<3; d++) {
      settings.width_a[d] /= voxel_width[d];
      settings.width_b[d] /= voxel_width[d];
      settings.log_width[d] /= voxel_width[d];
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

    // Now that we know the voxel_width, rescale the blob diameters(for drawing)
    // (Again, I assume the voxel width is the same in the x,y,z directions.)
    for (size_t ir = 0; ir < settings.blob_diameters.size(); ir++)
      settings.blob_diameters[ir] /= voxel_width[0];
    if (! settings.sphere_decals_shell_thickness_is_ratio)
      settings.sphere_decals_shell_thickness /= voxel_width[0];
    else
      settings.sphere_decals_shell_thickness /= settings.resize_with_binning;

    // Now that we know the voxel_width, rescale the coordinates
    // we have read from various files:
    for (size_t i = 0; i < settings.training_pos_crds.size(); i++)
      for (int d = 0; d < 3; d++)
        if (! settings.is_training_pos_in_voxels)
          settings.training_pos_crds[i][d] /= voxel_width[d];
        else
          settings.training_pos_crds[i][d] /= settings.resize_with_binning;

    for (size_t i = 0; i < settings.training_neg_crds.size(); i++)
      for (int d = 0; d < 3; d++)
        if (! settings.is_training_neg_in_voxels)
          settings.training_neg_crds[i][d] /= voxel_width[d];
        else
          settings.training_neg_crds[i][d] /= settings.resize_with_binning;

    for (int I=0; I < settings.multi_is_training_pos_in_voxels.size(); I++)
      for (size_t i = 0; i < settings.training_pos_crds.size(); i++)
        for (int d = 0; d < 3; d++)
          if (! settings.multi_is_training_pos_in_voxels[I])
            settings.multi_training_pos_crds[I][i][d] /= voxel_width[d];
          else
            settings.multi_training_pos_crds[I][i][d] /=
              settings.resize_with_binning;

    for (int I=0; I < settings.multi_is_training_neg_in_voxels.size(); I++)
      for (size_t i = 0; i < settings.training_neg_crds.size(); i++)
        for (int d = 0; d < 3; d++)
          if (! settings.multi_is_training_neg_in_voxels[I])
            settings.multi_training_neg_crds[I][i][d] /= voxel_width[d];
          else
            settings.multi_training_neg_crds[I][i][d] /=
              settings.resize_with_binning;

    for (size_t i = 0; i < settings.must_link_constraints.size(); i++)
      for (size_t j = 0; j < settings.must_link_constraints[i].size(); j++)
        for (int d = 0; d < 3; d++)
          if (! settings.is_must_link_constraints_in_voxels)
            settings.must_link_constraints[i][j][d] /= voxel_width[d];
          else
            settings.must_link_constraints[i][j][d] /=
              settings.resize_with_binning;




    // (rescale the brightness of the image so it lies between 0 and 1?)
    if (settings.rescale_min_max_in) {
      tomo_in.Rescale01(mask.aaafI,
                        settings.in_rescale_min,
                        settings.in_rescale_max);
      tomo_in.FindMinMaxMean();
    }




    // ---- make an array that will store the new tomogram we will create ----

    cerr << "allocating space for new 3D image..." << endl;
    MrcSimple tomo_out = tomo_in; //this will take care of allocating the array

    if (((voxel_width[0] <= 0.0) ||
         (voxel_width[1] <= 0.0) ||
         (voxel_width[2] <= 0.0))
        &&
        (!
         (settings.filter_type == Settings::BLOB_NONMAX_SUPPRESSION) ||
         (settings.filter_type == Settings::BLOB_NONMAX_SUPERVISED_MULTI)))
      throw VisfdErr("Error in tomogram header: Invalid voxel width(s).\n"
                     "Use the -w argument to specify the voxel width.");

    // If the user customized the voxel width (using the "-w" argument,
    // for example), then include this information in the header portion
    // of the output tomogram.
    if (voxel_width[0] != tomo_in.header.cellA[0]/image_size[0])
      tomo_out.header.cellA[0] = image_size[0] * voxel_width[0];
    if (voxel_width[1] != tomo_in.header.cellA[1]/image_size[1])
      tomo_out.header.cellA[1] = image_size[1] * voxel_width[1];
    if (voxel_width[2] != tomo_in.header.cellA[2]/image_size[2])
      tomo_out.header.cellA[2] = image_size[2] * voxel_width[2];




    // ---- Select the primary operation we will perform on the image ----



    if (settings.filter_type == Settings::NONE) {

      cerr << "filter_type = Intensity Map <No convolution filter specified>\n";

    } 



    else if (settings.filter_type == Settings::FIND_EXTREMA) {

      // Find the local minima or maxima in the image
      HandleExtrema(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::DILATION) {

      // Apply a greyscale dilation filter to the image.
      HandleDilation(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::EROSION) {

      // Apply a greyscale erosion filter to the image.
      HandleErosion(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::OPENING) {

      // Apply a greyscale opening filter to the image.
      HandleOpening(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::CLOSING) {

      // Apply a greyscale closing filter to the image.
      HandleClosing(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::TOP_HAT_WHITE) {

      // Apply a greyscale "white top-hat" filter to the image.
      HandleTopHatWhite(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::TOP_HAT_BLACK) {

      // Apply a greyscale "black top-hat" filter to the image.
      HandleTopHatBlack(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::MEDIAN) {

      #ifndef CXX17_UNSUPPORTED
      // Apply a median filter to the image.
      HandleMedian(settings, tomo_in, tomo_out, mask, voxel_width);
      #endif

    }



    else if (settings.filter_type == Settings::GAUSS) {

      // Apply a Gaussian blur filter to the image.
      HandleGauss(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::GGAUSS) {

      // Apply a generalized Gaussian filter
      HandleGGauss(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::DOG) {

      // Apply a Difference-of-Gaussians filter
      HandleDog(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::DOGG) {

      // Apply a generalized DoG filter.
      //   (Note: This kind of operation is probably not useful.
      //          I may remove this feature in the future. -Andrew 2021-7-11)
      HandleDogg(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::DOGGXY) {

      #ifndef DISABLE_DOGGXY
      // Apply a generalized DoG filter in the XY direction
      // and a Gaussian filter in the Z direction.
      HandleDoggXY(settings, tomo_in, tomo_out, mask, voxel_width);
      #endif

    }



    else if (settings.filter_type == Settings::LOG_DOG) {

      // Apply a LoG filter (implemented using the DoG approximation).
      HandleLoGDoG(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::BLOB) {

      HandleBlobDetector(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::BLOB_NONMAX_SUPPRESSION) {

      vector<array<float,3> > crds;
      vector<float> diameters;
      vector<float> scores;

      // Discard overlapping or poor scoring blobs from a (user-supplied) list.
      // (Note: In this case, we are not analyzing the image.
      //        The list of blobs is supplied by the user, most likely by
      //        running this program at an earlier time using the "-blob"
      //        argument.)
      HandleBlobsNonmaxSuppression(settings,
                                   mask,
                                   voxel_width,
                                   crds,
                                   diameters,
                                   scores);

    }


    else if (settings.filter_type == Settings::BLOB_NONMAX_SUPERVISED_MULTI) {

      // Use training data to choose the right threshold(s) for discarding blobs
      HandleBlobScoreSupervisedMulti(settings,
                                     voxel_width);
    }



    else if ((settings.filter_type == Settings::SURFACE_EDGE) ||
             (settings.filter_type == Settings::SURFACE_RIDGE) ||
             (settings.filter_type == Settings::CURVE)) {

      // Detect 2-D surfaces or 1-D curves
      HandleTV(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::WATERSHED) {

      // Perform watershed segmentation
      HandleWatershed(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::LABEL_CONNECTED) {

      // Cluster adjacent nearby voxels with high "saliency" into "islands"
      // neighboring voxels (this is similar to watershed segmentation).
      HandleLabelConnected(settings, tomo_in, tomo_out, mask, voxel_width);

    }



    else if (settings.filter_type == Settings::LOCAL_FLUCTUATIONS) {

      // Apply a filter that measures the amount of brightness fluctuations
      // in the local neighborhood.
      HandleLocalFluctuations(settings, tomo_in, tomo_out, mask, voxel_width);

    }


    // ----- distance calculations -----

    else if (settings.filter_type == Settings::DISTANCE_TO_POINTS) {

      // Apply a filter which replaces the voxel brightness with the distance
      // to the nearest poing in a (user-supplied) point cloud.
      HandleDistanceToPoints(settings, tomo_in, tomo_out, mask, voxel_width);

    }


    else if (settings.filter_type == Settings::DISTANCE_TO_VOXELS) {

      // Read the coordinates of multiple points in 3-D (a pointcloud) and
      // calculate the distance of each of these points to the nearest voxel
      // whose brightness lies in a range selected by the user.
      HandleDistancePointsToFeature(settings, tomo_in, tomo_out, mask, voxel_width);

    }


    // ----- random_spheres -----

    else if (settings.filter_type == Settings::RANDOM_SPHERES) {

      // Generate a random point cloud of thick spheres which fills a 3-D region
      // of space corresponding to a region in the image whose voxel brightness
      // lies in a range selected by the user.
      HandleRandomSpheres(settings, tomo_in, tomo_out, mask, voxel_width);

    }


    // ----- sphere_decals_filter -----

    else if (settings.filter_type == Settings::DRAW_SPHERES) {

      // Draw a new image or annotate (draw on top of) an existing image.
      HandleDrawSpheres(settings, tomo_in, tomo_out, mask, voxel_width);

    }




    // ------- DEPRECIATED FEATURES -------


    else if (settings.filter_type == Settings::TEMPLATE_GGAUSS) {

      #ifndef DISABLE_TEMPLATE_MATCHING
      // Spherical template matching with RMS error reporting.
      // (Not very useful. I may remove this feature in the future -A 2021-7-11)
      HandleTemplateGGauss(settings, tomo_in, tomo_out, mask, voxel_width);
      #endif //#ifndef DISABLE_TEMPLATE_MATCHING

    }


    else if (settings.filter_type == Settings::TEMPLATE_GAUSS) {

      #ifndef DISABLE_TEMPLATE_MATCHING
      // Spherical template matching with RMS error reporting.
      // (Not very useful. I may remove this feature in the future -A 2021-7-11)
      HandleTemplateGauss(settings, tomo_in, tomo_out, mask, voxel_width);
      #endif //#ifndef DISABLE_TEMPLATE_MATCHING

    }


    else if (settings.filter_type == Settings::BOOTSTRAP_DOGG) {

      #ifdef DISABLE_BOOSTRAPPING
      // (I may remove this feature in the future -Andrew 2021-7-11)
      HandleBootstrappDogg(settings, tomo_in, tomo_out, mask, voxel_width);
      #endif //#ifndef DISABLE_BOOTSTRAPPING

    }


    else if (settings.filter_type == Settings::BLOB_RADIAL_INTENSITY) {

      #ifndef DISABLE_INTENSITY_PROFILES
      // This is a feature that a user requested for a specific project.
      // It is probably not relevant to most users.
      // (I may remove this feature in the future -Andrew 2021-7-11)
      HandleBlobRadialIntensity(settings, tomo_in, tomo_out, mask, voxel_width);
      #endif //#ifndef DISABLE_TEMPLATE_MATCHING

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





    // ----- masking: -----
    // Also, after thresholding, check again to see if the mask is zero
    // at this location.  If so, make sure that aaafI is 0 there
    // (or whatever the user has requested us to put there).
    // Do this after thresholding and inverting.
    if ((mask.aaafI) && (settings.specify_masked_brightness))
      for (int iz=0; iz<mask.header.nvoxels[2]; iz++)
        for (int iy=0; iy<mask.header.nvoxels[1]; iy++)
          for (int ix=0; ix<mask.header.nvoxels[0]; ix++)
            if (mask.aaafI[iz][iy][ix] == 0.0)
              tomo_out.aaafI[iz][iy][ix] = settings.masked_voxel_brightness;

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


