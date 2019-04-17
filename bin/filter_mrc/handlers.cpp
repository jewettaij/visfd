#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <array>
//#include <fftw3.h>  not needed yet
using namespace std;

#ifndef DISABLE_OPENMP
#include <omp.h>       // (OpenMP-specific)
#endif

#include <visfd.hpp>
#include <threshold.hpp>
#include <mrc_simple.hpp>
#include <random_gen.h>
#include "settings.hpp"
#include "file_io.hpp"
#include "filter3d_variants.hpp"
#include "feature_variants.hpp"
#include "feature_unsupported.hpp"
#include "handlers.hpp"
#include "handlers_unsupported.hpp"


void
HandleGGauss(Settings settings,
             MrcSimple &tomo_in,
             MrcSimple &tomo_out,
             MrcSimple &mask,
             float voxel_width[3])
{
  Filter3D<float, int>
    filter = GenFilterGenGauss3D(settings.width_a,
                                 settings.m_exp,
                                 settings.filter_truncate_ratio,
                                 settings.filter_truncate_threshold,
                                 static_cast<float*>(NULL),
                                 &cerr);

  filter.Apply(tomo_in.header.nvoxels,
               tomo_in.aaafI,
               tomo_out.aaafI,
               mask.aaafI,
               true,
               &cerr);

} //HandleGGauss()




void
HandleGauss(Settings settings,
            MrcSimple &tomo_in,
            MrcSimple &tomo_out,
            MrcSimple &mask,
            float voxel_width[3])
{
  cerr << "filter_type = Gaussian\n";

  float A; // let the user know what A coefficient was used

  A = ApplyGauss(tomo_in.header.nvoxels,
                 tomo_in.aaafI,
                 tomo_out.aaafI,   // <-- save result here
                 mask.aaafI,
                 settings.width_a,
                 settings.filter_truncate_ratio,
                 settings.filter_truncate_threshold,
                 true,
                 &cerr);

      
  cerr << " Filter Used:\n"
    " h(x,y,z)   = A*exp(-0.5*((x/σ_x)^2 + (y/σ_y)^2 + (z/σ_z)^))\n"
    " ... where  A = " << A << "\n" 
    "   (σ_x, σ_y, σ_z) = "
       << "(" << settings.width_a[0]
       << " " << settings.width_a[1]
       << " " << settings.width_a[2] << ")\n";
  cerr << " You can plot a slice of this function\n"
       << "     in the X direction using:\n"
       << " draw_filter_1D.py -gauss "
       << A << " " << settings.width_a[0] << endl;
  cerr << " and in the Y direction using:\n"
       << " draw_filter_1D.py -gauss "
       << A << " " << settings.width_a[1] << endl;
  cerr << " and in the Z direction using:\n"
       << " draw_filter_1D.py -gauss "
       << A << " " << settings.width_a[2] << endl;

} //HandleGauss()




void
HandleDogg(Settings settings,
           MrcSimple &tomo_in,
           MrcSimple &tomo_out,
           MrcSimple &mask,
           float voxel_width[3])
{
  cerr << "filter_type = Difference-of-Generalized-Gaussians (DOGG)\n";

  Filter3D<float, int> filter;
  float A, B;       // let the user know what A B coefficients were used

  filter = GenFilterDogg3D(settings.width_a,//"a" parameter in formula
                           settings.width_b,//"b" parameter in formula
                           settings.m_exp,  //"m" parameter in formula
                           settings.n_exp,  //"n" parameter in formula
                           settings.filter_truncate_ratio,
                           settings.filter_truncate_threshold,
                           &A,
                           &B,
                           &cerr);

  filter.Apply(tomo_in.header.nvoxels,
               tomo_in.aaafI,
               tomo_out.aaafI,
               mask.aaafI,
               false,
               &cerr);

} // HandleDogg()




void
HandleDog(Settings settings,
          MrcSimple &tomo_in,
          MrcSimple &tomo_out,
          MrcSimple &mask,
          float voxel_width[3])
{
  cerr << "filter_type = Difference of Gaussians (DOG)\n";

  float A, B;

  ApplyDog(tomo_in.header.nvoxels,
           tomo_in.aaafI,
           tomo_out.aaafI,
           mask.aaafI,
           settings.width_a,
           settings.width_b,
           settings.filter_truncate_ratio,
           settings.filter_truncate_threshold,
           &A,
           &B,
           &cerr);

  cerr << " Filter Used:\n"
    " h(x,y,z)   = h_a(x,y,z) - h_b(x,y,z)\n"
    " h_a(x,y,z) = A*exp(-0.5*((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))\n"
    " h_b(x,y,z) = B*exp(-0.5*((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))\n"
    "  ... where      A = " << A << "\n"
    "                 B = " << B << "\n" 
    "   (a_x, a_y, a_z) = "
       << "(" << settings.width_a[0]
       << " " << settings.width_a[1]
       << " " << settings.width_a[2] << ")\n"
    "   (b_x, b_y, b_z) = "
       << "(" << settings.width_b[0]
       << " " << settings.width_b[1]
       << " " << settings.width_b[2] << ")\n";
  cerr << " You can plot a slice of this function\n"
       << "     in the X direction using:\n"
    " draw_filter_1D.py -dog " << A << " " << B
       << " " << settings.width_a[0] << " " << settings.width_b[0] << endl;
  cerr << " and in the Y direction using:\n"
    " draw_filter_1D.py -dog " << A << " " << B
       << " " << settings.width_a[1] << " " << settings.width_b[1] << endl;
  cerr << " and in the Z direction using:\n"
    " draw_filter_1D.py -dog " << A << " " << B
       << " " << settings.width_a[2] << " " << settings.width_b[2] << endl;

} //HandleDog()





void
HandleDogScaleFree(Settings settings,
                   MrcSimple &tomo_in,
                   MrcSimple &tomo_out,
                   MrcSimple &mask,
                   float voxel_width[3])
{
  cerr << "filter_type = (Fast) Laplacian of Gaussians (LoG)\n"
       << "  (This will be approximated as a Difference of Gaussians,\n"
       << "   as explained below.)\n";
  //cerr << "filter_type = Difference of Gaussians Scale Free (DOGSF)\n";

  float A, B;

  ApplyLog(tomo_in.header.nvoxels,
           tomo_in.aaafI,
           tomo_out.aaafI,
           mask.aaafI,
           settings.dogsf_width,
           settings.delta_sigma_over_sigma,
           settings.filter_truncate_ratio,
           settings.filter_truncate_threshold,
           &A,
           &B,
           &cerr);

  cerr << " Filter Used:\n"
    " h(x,y,z)   = h_a(x,y,z) - h_b(x,y,z)\n"
    " h_a(x,y,z) = A*exp(-0.5*((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))\n"
    " h_b(x,y,z) = B*exp(-0.5*((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))\n"
    "  ... where      A = " << A << "\n"
    "                 B = " << B << "\n" 
    "   (a_x, a_y, a_z) = "
       << "(" << settings.dogsf_width[0]*(1-0.5*settings.delta_sigma_over_sigma)
       << " " << settings.dogsf_width[1]*(1-0.5*settings.delta_sigma_over_sigma)
       << " " << settings.dogsf_width[2]*(1-0.5*settings.delta_sigma_over_sigma)
       << ")\n"
    "   (b_x, b_y, b_z) = "
       << "(" << settings.dogsf_width[0]*(1+0.5*settings.delta_sigma_over_sigma)
       << " " << settings.dogsf_width[1]*(1+0.5*settings.delta_sigma_over_sigma)
       << " " << settings.dogsf_width[2]*(1+0.5*settings.delta_sigma_over_sigma)
       << ")\n";
  cerr << " You can plot a slice of this function\n"
       << "     in the X direction using:\n"
    " draw_filter_1D.py -dog " << A << " " << B
       << " " << settings.dogsf_width[0]*(1-0.5*settings.delta_sigma_over_sigma)
       << " " << settings.dogsf_width[0]*(1+0.5*settings.delta_sigma_over_sigma)
       << endl;
  cerr << " and in the Y direction using:\n"
    " draw_filter_1D.py -dog " << A << " " << B
       << " " << settings.dogsf_width[1]*(1-0.5*settings.delta_sigma_over_sigma)
       << " " << settings.dogsf_width[1]*(1+0.5*settings.delta_sigma_over_sigma)
       << endl;
  cerr << " and in the Z direction using:\n"
    " draw_filter_1D.py -dog " << A << " " << B
       << " " << settings.dogsf_width[2]*(1-0.5*settings.delta_sigma_over_sigma)
       << " " << settings.dogsf_width[2]*(1+0.5*settings.delta_sigma_over_sigma)
       << endl;

} // HandleDogScaleFree()








void
HandleBlobsNonmaxSuppression(Settings settings,
                             MrcSimple &mask,
                             float voxel_width[3],
                             vector<array<float,3> >& crds,
                             vector<float>& diameters,
                             vector<float>& scores)
{
  // At some point I was trying to be as general as possible and allowed
  // for the possibility that voxels need not be cubes (same width x,y,z)
  // Now, I realize that allowing for this possibility would slow down some
  // calculations considerably, so I just assume cube-shaped voxels:
  float voxel_width_ = voxel_width[0];
  assert((voxel_width[0] == voxel_width[1]) &&
         (voxel_width[1] == voxel_width[2]));


  ReadBlobCoordsFile(settings.in_coords_file_name,
                     &crds,
                     &diameters,
                     &scores,
                     voxel_width[0],
                     settings.sphere_decals_diameter,
                     settings.sphere_decals_scale,
                     settings.sphere_diameters_lower_bound,
                     settings.sphere_diameters_upper_bound,
                     settings.sphere_decals_foreground,
                     settings.score_lower_bound,
                     settings.score_upper_bound);

  cerr << " --- discarding blobs in file \n"
       << " \"" << settings.in_coords_file_name << "\" ---\n"
       << "\n";

  if ((crds.size() > 0) && (mask.aaafI != NULL)) {
    cerr << "  discarding blobs outside the mask" << endl;
    DiscardMaskedBlobs(crds,
                       diameters,
                       scores,
                       mask.aaafI);
  }

  cerr << "  discarding overlapping blobs" << endl;
  DiscardOverlappingBlobs(crds,
                          diameters, 
                          scores,
                          DO_NOT_SORT,
                          //REMOVE THIS CRUFT
                          //settings.find_extrema_occlusion_ratio,
                          settings.nonmax_min_radial_separation_ratio,
                          settings.nonmax_max_volume_overlap_large,
                          settings.nonmax_max_volume_overlap_small,
                          &cerr);

  cerr << " " << crds.size() << " blobs remaining" << endl;

  if (settings.out_coords_file_name != "") {
    fstream out_coords_file;
    out_coords_file.open(settings.out_coords_file_name.c_str(), ios::out);
    for (size_t i=0; i < crds.size(); i++) {
      out_coords_file << crds[i][0]*voxel_width[0] << " "
                      << crds[i][1]*voxel_width[1] << " "
                      << crds[i][2]*voxel_width[2] << " " 
                      << diameters[i]*voxel_width[0] << " " 
                      << scores[i]
                      << endl;
    }
  }

} //HandleBlobsNonmaxSuppression()





void
HandleVisualizeBlobs(Settings settings,
                     MrcSimple &tomo_in,
                     MrcSimple &tomo_out,
                     MrcSimple &mask,
                     float voxel_width[3])
{

  vector<array<float,3> > crds;
  vector<float> diameters;
  vector<float> scores;

  HandleBlobsNonmaxSuppression(settings,
                               mask,
                               voxel_width,
                               crds,
                               diameters,
                               scores);

  assert(crds.size() == diameters.size());
  assert(crds.size() == scores.size());
  size_t n_blobs = diameters.size();

  // Sometimes the user wants to display the decal with a particular brightness
  // We store that brightness inside the "scores[]" array
  if (! settings.sphere_decals_foreground_use_score)
    for (size_t i=0; i < n_blobs; i++)
      scores[i] = settings.sphere_decals_foreground;

  vector<float> shell_thicknesses(n_blobs);
  for (size_t i=0; i < n_blobs; i++) {
    float diameter = diameters[i];
    float shell_thickness = settings.sphere_decals_shell_thickness;
    if (settings.sphere_decals_shell_thickness_is_ratio) {
      shell_thickness *= diameter;
      if (shell_thickness < settings.sphere_decals_shell_thickness_min)
        shell_thickness = settings.sphere_decals_shell_thickness_min;
    }
    shell_thickness = settings.sphere_decals_shell_thickness;
    if (settings.sphere_decals_shell_thickness_is_ratio) {
      shell_thickness *= diameter;
      // The spherical shells superimposed on the tomogram should be
      // at least 1 voxel wide in order to be visible to the user
      if (shell_thickness < settings.sphere_decals_shell_thickness_min)
        shell_thickness = 1.0;
    }
    shell_thicknesses[i] = shell_thickness;
  }

  reverse(crds.begin(), crds.end());
  reverse(diameters.begin(), diameters.end());
  reverse(shell_thicknesses.begin(), shell_thicknesses.end());
  reverse(scores.begin(), scores.end());

  VisualizeBlobs(tomo_out.header.nvoxels,
                 tomo_out.aaafI,
                 mask.aaafI,
                 crds,
                 &diameters,
                 &shell_thicknesses,
                 &scores,
                 settings.sphere_decals_background,
                 tomo_in.aaafI,
                 settings.sphere_decals_background_scale,
                 settings.sphere_decals_foreground_norm);

} //HandleVisualizeBlobs()






void
HandleBlobDetector(Settings settings,
                   MrcSimple &tomo_in,
                   MrcSimple &tomo_out,
                   MrcSimple &mask,
                   float voxel_width[3])
{
  vector<array<float,3> > minima_crds_voxels;
  vector<array<float,3> > maxima_crds_voxels;

  vector<float> minima_diameters;
  vector<float> maxima_diameters;

  vector<float> minima_scores;
  vector<float> maxima_scores;


  // Optional: Preallocate space for BlobDogNM()
  float*** aaaafI[3];
  float* aafI[3];

  Alloc3D(tomo_in.header.nvoxels,
          &(aafI[0]),
          &(aaaafI[0]));

  if (tomo_out.aaafI) {
    // Optional: Instead of allocating aaaafI[1] and aafI[1], borrow the memory
    // you've already allocated to tomo_out.aaafI and afI. (Goal: Save memory.)
    aafI[1] = tomo_out.afI;
    aaaafI[1] = tomo_out.aaafI;
  } else {
    Alloc3D(tomo_in.header.nvoxels,
            &(aafI[1]),
            &(aaaafI[1]));
  }

  Alloc3D(tomo_in.header.nvoxels,
          &(aafI[2]),
          &(aaaafI[2]));
  // This way we can save memory and also save the 
  // filtered image to a file which we can view using IMOD.

  _BlobDogNM(tomo_in.header.nvoxels,
             tomo_in.aaafI,
             mask.aaafI,
             settings.blob_diameters,  // try detecting blobs of these diameters
             minima_crds_voxels,  // store minima x,y,z coords here
             maxima_crds_voxels,  // store maxima x,y,z coords here
             minima_diameters, // corresponding diameter for that minima
             maxima_diameters, // corresponding diameter for that maxima
             minima_scores, // what was the blob's score?
             maxima_scores, // ("score" = intensity after filtering)
             settings.delta_sigma_over_sigma, //difference in Gauss widths parameter
             settings.filter_truncate_ratio,
             settings.filter_truncate_threshold,
             settings.score_upper_bound,
             settings.score_lower_bound,
             settings.score_bounds_are_ratios,
             settings.nonmax_min_radial_separation_ratio,
             settings.nonmax_max_volume_overlap_large,
             settings.nonmax_max_volume_overlap_small,
             &cerr,
             aaaafI,
             aafI);


  long n_minima = minima_crds_voxels.size();
  long n_maxima = maxima_crds_voxels.size();
  vector<array<float,3> > minima_crds(n_minima);
  vector<array<float,3> > maxima_crds(n_maxima);

  // The user expects results in units of physical distance, not voxels.
  // Compensate for that now:
  for (int i = 0; i < n_minima; i++) {
    minima_crds[i][0] = minima_crds_voxels[i][0] * voxel_width[0];
    minima_crds[i][1] = minima_crds_voxels[i][1] * voxel_width[1];
    minima_crds[i][2] = minima_crds_voxels[i][2] * voxel_width[2];
    minima_diameters[i] *= voxel_width[0];      //Gaussian width has units of length
  }
  for (int i = 0; i < n_maxima; i++) {
    maxima_crds[i][0] = maxima_crds_voxels[i][0] * voxel_width[0];
    maxima_crds[i][1] = maxima_crds_voxels[i][1] * voxel_width[1];
    maxima_crds[i][2] = maxima_crds_voxels[i][2] * voxel_width[2];
    maxima_diameters[i] *= voxel_width[0];      //Gaussian width has units of length
  }

  //string out_file_name_base = settings.out_file_name;
  //if ((EndsWith(settings.out_file_name, ".rec")) ||
  //    (EndsWith(settings.out_file_name, ".mrc")))
  //  out_file_name_base =
  //    settings.out_file_name.substr(0,
  //                                  settings.out_file_name.length()-4);

  if ((minima_crds_voxels.size() > 0) && (settings.blob_minima_file_name != ""))
  {

    SortBlobs(minima_crds,
              minima_diameters,
              minima_scores,
              false);

    fstream minima_file;
    minima_file.open(settings.blob_minima_file_name.c_str(), ios::out);
    if (! minima_file)
      throw VisfdErr("Error: unable to open \""+ settings.blob_minima_file_name +"\" for reading.\n");
    for (int i=0; i < minima_crds_voxels.size(); i++) {
      minima_file << minima_crds[i][0] << " "
                  << minima_crds[i][1] << " "
                  << minima_crds[i][2] << " "
                  << minima_diameters[i] << " "
                //<< minima_diameters[i] / settings.blob_width_multiplier << " "
        // See comment at the end of this function:
                  << minima_scores[i] << "\n";
        // Here we assume the voxel width is the same in x,y,z directions
        // and equals voxel_width[0].
    }
  }


  if ((maxima_crds_voxels.size() > 0) && (settings.blob_maxima_file_name != ""))
  {

    SortBlobs(maxima_crds,
              maxima_diameters,
              maxima_scores,
              true);

    fstream maxima_file;
    maxima_file.open(settings.blob_maxima_file_name.c_str(), ios::out);
    if (! maxima_file)
      throw VisfdErr("Error: unable to open \""+ settings.blob_maxima_file_name +"\" for reading.\n");
    for (int i=0; i < maxima_crds_voxels.size(); i++) {
      maxima_file << maxima_crds[i][0] << " "
                  << maxima_crds[i][1] << " "
                  << maxima_crds[i][2] << " "
                  << maxima_diameters[i] << " "
                //<< maxima_diameters[i] / settings.blob_width_multiplier << " "
        // See comment at the end of this function:
                  << maxima_scores[i] << "\n";
        // Here we assume the voxel width is the same in x,y,z directions
        // and equals voxel_width[0].
    }
  }

  // ---------------------------------------------------------------------
  // OPTIONAL: Now render the final image with all of the spherical shells
  //           superimposed on top of the original image:
  //
  // Does the user want us to create a new image showing where the blobs are?
  if (tomo_out.aaafI) {

    vector<array<float,3> > display_crds_voxels(minima_crds_voxels);
    display_crds_voxels.insert(display_crds_voxels.end(),
                               maxima_crds_voxels.rbegin(),
                               maxima_crds_voxels.rend());

    //concatinate "maxima_diameters" with "minima_diameters"
    vector<float> display_diameters(minima_diameters);
    display_diameters.insert(display_diameters.end(),
                             maxima_diameters.rbegin(),
                             maxima_diameters.rend());

    //concatinate "maxima_scores" with "minima_scores"
    vector<float> display_scores(minima_scores);
    display_scores.insert(display_scores.end(),
                          maxima_scores.rbegin(),
                          maxima_scores.rend());

    vector<float> display_shell_thicknesses(display_crds_voxels.size());
    for (int i = 0; i < display_crds_voxels.size(); i++) {
      // Choose the size of the hollow spheres around each object so they are
      // large enough that the original objects underneath are still visible.
      display_diameters[i] = display_diameters[i] / voxel_width[0];
      display_shell_thicknesses[i] = settings.sphere_decals_shell_thickness;
      if (settings.sphere_decals_shell_thickness_is_ratio)
        display_shell_thicknesses[i] *= display_diameters[i];
      display_diameters[i] *= settings.sphere_decals_scale;
      // The spherical shells superimposed on the tomogram should be
      // at least 1 voxel wide in order to be visible to the user
      if (display_shell_thicknesses[i] < settings.sphere_decals_shell_thickness_min)
        display_shell_thicknesses[i] = 1.0;
    }

    VisualizeBlobs(tomo_out.header.nvoxels,
                   tomo_out.aaafI,
                   mask.aaafI,
                   display_crds_voxels,
                   &display_diameters,
                   &display_shell_thicknesses,
                   &display_scores,
                   settings.sphere_decals_background,
                   tomo_in.aaafI,
                   settings.sphere_decals_background_scale,
                   false);

  } //if (tomo_out.aaafI)



  Dealloc3D(tomo_in.header.nvoxels,
            &aafI[0],
            &aaaafI[0]);

  // Optional: Instead of allocating aaaafI[1] and aafI[1], borrow the memory
  // you've already allocated to tomo_out.aaafI and afI. (Goal: Save memory.)
  if (! tomo_out.aaafI) {  // <-- Was tomo_out.aaafI available?
    // However, if we didn't allocate it this memory, then don't delete it:
    Dealloc3D(tomo_in.header.nvoxels,
              &aafI[1],
              &aaaafI[1]);
  }

  Dealloc3D(tomo_in.header.nvoxels,
            &aafI[2],
            &aaaafI[2]);


} //HandleBlobDetector()






void
HandleThresholds(Settings settings,
                 MrcSimple &tomo_in,
                 MrcSimple &tomo_out,
                 MrcSimple &mask,
                 float voxel_width[3])
{
  cerr << "Applying thresholds" << endl;

  if (settings.out_thresh2_use_clipping_sigma) {
    float stddev_intensity = StdDevArr(tomo_out.header.nvoxels,
                                       tomo_out.aaafI,
                                       mask.aaafI);
    float ave_intensity = AverageArr(tomo_out.header.nvoxels,
                                     tomo_out.aaafI,
                                     mask.aaafI);

    settings.out_threshold_01_a = ave_intensity +
      settings.out_threshold_01_a * stddev_intensity;
    settings.out_threshold_01_b = ave_intensity +
      settings.out_threshold_01_b * stddev_intensity;
    cerr << "ave="<< ave_intensity <<", stddev="<<stddev_intensity << endl;
    cerr << "  Clipping intensities between ["
         << settings.out_threshold_01_a << ", "
         << settings.out_threshold_01_b << "]" << endl;
  }
  for (int iz=0; iz<tomo_out.header.nvoxels[2]; iz++) {
    for (int iy=0; iy<tomo_out.header.nvoxels[1]; iy++) {
      for (int ix=0; ix<tomo_out.header.nvoxels[0]; ix++) {
        if (settings.use_rescale_multiply) {
          tomo_out.aaafI[iz][iy][ix] *= settings.out_rescale_multiply;
          tomo_out.aaafI[iz][iy][ix] += settings.out_rescale_offset;
        }
        else if (settings.use_gauss_thresholds)
          tomo_out.aaafI[iz][iy][ix] =
            SelectIntensityRangeGauss(tomo_out.aaafI[iz][iy][ix],
                                      settings.out_thresh_gauss_x0,
                                      settings.out_thresh_gauss_sigma,
                                      settings.out_thresh_a_value,
                                      settings.out_thresh_b_value);
        else if (! settings.use_dual_thresholds) {
          if (settings.out_threshold_01_a == settings.out_threshold_01_b)
            tomo_out.aaafI[iz][iy][ix] = ((tomo_out.aaafI[iz][iy][ix] >
                                           settings.out_threshold_01_a)
                                          ? 1.0
                                          : 0.0);
          else
            tomo_out.aaafI[iz][iy][ix] =
              Threshold2(tomo_out.aaafI[iz][iy][ix],
                         settings.out_threshold_01_a,
                         settings.out_threshold_01_b,
                         (settings.out_thresh2_use_clipping
                          ? settings.out_threshold_01_a
                          : settings.out_thresh_a_value),
                         (settings.out_thresh2_use_clipping
                          ? settings.out_threshold_01_b
                          : settings.out_thresh_b_value));
        }
        else
          tomo_out.aaafI[iz][iy][ix] =
            Threshold4(tomo_out.aaafI[iz][iy][ix],
                       settings.out_threshold_01_a,
                       settings.out_threshold_01_b,
                       settings.out_threshold_10_a,
                       settings.out_threshold_10_b,
                       settings.out_thresh_a_value,
                       settings.out_thresh_b_value);
      }
    }
  }
} //HandleThresholds()



void
HandleExtrema(Settings settings,
              MrcSimple &tomo_in,
              MrcSimple &tomo_out,
              MrcSimple &mask,
              float voxel_width[3])
{
  assert(tomo_in.aaafI);
  //hopefully the next line is unnecessary. (it should already have been done)
  tomo_in.FindMinMaxMean(); //update the min, max, mean header params

  // optional: zero out the output image in case the user wants us to
  //           create an image showing where the voxels are
  assert(tomo_out.header.nvoxels[0] = tomo_in.header.nvoxels[0]);
  assert(tomo_out.header.nvoxels[1] = tomo_in.header.nvoxels[1]);
  assert(tomo_out.header.nvoxels[2] = tomo_in.header.nvoxels[2]);
  assert(tomo_out.aaafI);
  for (int iz = 0; iz < tomo_out.header.nvoxels[2]; iz++)
    for (int iy = 0; iy < tomo_out.header.nvoxels[1]; iy++)
      for (int ix = 0; ix < tomo_out.header.nvoxels[0]; ix++)
        tomo_out.aaafI[iz][iy][ix] = 0.0;

  //default values disable the thresholds:
  // REMOVE THIS CRUFT:
  //float local_minima_threshold = settings.find_minima_threshold;
  //float local_maxima_threshold = settings.find_maxima_threshold;

  float minima_threshold = settings.score_upper_bound;
  float maxima_threshold = settings.score_lower_bound;
  
  vector<array<float, 3> > minima_crds_voxels;
  vector<array<float, 3> > maxima_crds_voxels;
  vector<float> minima_scores;
  vector<float> maxima_scores;
  vector<size_t> minima_nvoxels;
  vector<size_t> maxima_nvoxels;
  vector<array<float, 3> > *pv_minima_crds_voxels = NULL;
  vector<array<float, 3> > *pv_maxima_crds_voxels = NULL;
  vector<float> *pv_minima_scores = NULL;
  vector<float> *pv_maxima_scores = NULL;
  vector<size_t> *pv_minima_nvoxels = NULL;
  vector<size_t> *pv_maxima_nvoxels = NULL;


  if (settings.find_minima) {
    pv_minima_crds_voxels = &minima_crds_voxels;
    pv_minima_scores = &minima_scores;
    pv_minima_nvoxels = &minima_nvoxels;
  }

  if (settings.find_maxima) {
    pv_maxima_crds_voxels = &maxima_crds_voxels;
    pv_maxima_scores = &maxima_scores;
    pv_maxima_nvoxels = &maxima_nvoxels;
  }

  FindExtrema(tomo_in.header.nvoxels,
              tomo_in.aaafI,
              mask.aaafI,
              pv_minima_crds_voxels,
              pv_maxima_crds_voxels,
              pv_minima_scores,
              pv_maxima_scores,
              pv_minima_nvoxels,
              pv_maxima_nvoxels,
              minima_threshold,
              maxima_threshold,
              settings.neighbor_connectivity,
              settings.extrema_on_boundary,
              tomo_out.aaafI, //<--an image showing where the minima are?
              &cerr);

  // non-max suppression
  // discards minima or maxima which lie too close together.
  // This requires that we choose a size for each minima or maxima.
  // (If two minima/maxima lie within this distance, one of them is discarded.)
  float use_this_diameter = (settings.sphere_decals_diameter *
                             settings.nonmax_min_radial_separation_ratio);
  if (use_this_diameter <= 0.0)
    use_this_diameter = 0.0;
  vector<float> minima_diameters(minima_crds_voxels.size(), use_this_diameter);
  vector<float> maxima_diameters(maxima_crds_voxels.size(), use_this_diameter);

  if ((minima_crds_voxels.size() > 0) && (mask.aaafI != NULL))
    DiscardMaskedBlobs(minima_crds_voxels,
                       minima_diameters,
                       minima_scores,
                       mask.aaafI);

  if ((maxima_crds_voxels.size() > 0) && (mask.aaafI != NULL))
    DiscardMaskedBlobs(maxima_crds_voxels,
                       maxima_diameters,
                       maxima_scores,
                       mask.aaafI);

  if ((settings.nonmax_min_radial_separation_ratio > 0.0) ||
      (settings.nonmax_max_volume_overlap_large < 1.0) ||
      (settings.nonmax_max_volume_overlap_small < 1.0)) {
    if ((settings.sphere_decals_diameter > 0) &&
        //(settings.find_extrema_occlusion_ratio > 0.0)) {
        (settings.nonmax_min_radial_separation_ratio > 0.0)) {
      DiscardOverlappingBlobs(minima_crds_voxels,
                              minima_diameters, 
                              minima_scores,
                              PRIORITIZE_LOW_SCORES,
                              //settings.find_extrema_occlusion_ratio,
                              settings.nonmax_min_radial_separation_ratio,
                              settings.nonmax_max_volume_overlap_large,
                              settings.nonmax_max_volume_overlap_small,
                              &cerr);
    }

    if ((settings.sphere_decals_diameter > 0) &&
        //(settings.find_extrema_occlusion_ratio > 0.0)) {
        (settings.nonmax_min_radial_separation_ratio > 0.0)) {
      DiscardOverlappingBlobs(maxima_crds_voxels,
                              maxima_diameters, 
                              maxima_scores,
                              PRIORITIZE_HIGH_SCORES,
                              //REMOVE THIS CRUFT
                              //settings.find_extrema_occlusion_ratio,
                              settings.nonmax_min_radial_separation_ratio,
                              settings.nonmax_max_volume_overlap_large,
                              settings.nonmax_max_volume_overlap_small,
                              &cerr);
    }

  } //if (settings.nonmax_min_radial_separation_ratio > 0)

  //string out_file_name_base = settings.out_file_name;
  //if ((EndsWith(settings.out_file_name, ".rec")) ||
  //    (EndsWith(settings.out_file_name, ".mrc")))
  //  out_file_name_base =
  //    settings.out_file_name.substr(0,
  //                                  settings.out_file_name.length()-4);
  if ((minima_crds_voxels.size()) > 0 && settings.find_minima) {
    fstream minima_file;
    minima_file.open(settings.find_minima_file_name.c_str(), ios::out);
    if (! minima_file)
      throw VisfdErr("Error: unable to open \""+ settings.find_minima_file_name +"\" for reading.\n");
    for (int i=0; i < minima_crds_voxels.size(); i++)
      minima_file << minima_crds_voxels[i][0] * voxel_width[0] << " "
                  << minima_crds_voxels[i][1] * voxel_width[1] << " "
                  << minima_crds_voxels[i][2] * voxel_width[2] << " "
                  //<< minima_diameters[i] << " "
                  << minima_nvoxels[i] << " "
                  << minima_scores[i] << "\n";
  }
  if ((maxima_crds_voxels.size() > 0) && settings.find_maxima) {
    fstream coords_file;
    coords_file.open(settings.find_maxima_file_name.c_str(), ios::out);
    if (! coords_file)
      throw VisfdErr("Error: unable to open \""+ settings.find_maxima_file_name +"\" for reading.\n");
    for (int i=0; i < maxima_crds_voxels.size(); i++)
      coords_file << maxima_crds_voxels[i][0] * voxel_width[0] << " "
                  << maxima_crds_voxels[i][1] * voxel_width[1] << " "
                  << maxima_crds_voxels[i][2] * voxel_width[2] << " "
                  //<< maxima_diameters[i] << " "
                  << maxima_nvoxels[i] << " "
                  << maxima_scores[i] << "\n";
  }
} //HandleExtrema()







void
HandleLocalFluctuations(Settings settings,
                        MrcSimple &tomo_in,
                        MrcSimple &tomo_out,
                        MrcSimple &mask,
                        float voxel_width[3])
// Calculate the fluctuations of nearby voxel intensities
{
  LocalFluctuationsByRadius(tomo_in.header.nvoxels,
                            tomo_in.aaafI,
                            tomo_out.aaafI,
                            mask.aaafI,
                            settings.template_background_radius,
                            settings.template_background_exponent,
                            settings.filter_truncate_ratio,
                            settings.filter_truncate_threshold,
                            true,
                            &cerr);

  tomo_out.FindMinMaxMean();
}







void
HandleWatershed(Settings settings,
                MrcSimple &tomo_in,
                MrcSimple &tomo_out,
                MrcSimple &mask,
                float voxel_width[3])
{
  int image_size[3];
  for (int d = 0; d < 3; d++)
    image_size[d] = tomo_in.header.nvoxels[d];

  vector<array<int, 3> > extrema_crds;
  vector<float> extrema_scores;

  // Create a temporary array to store the basin membership for each voxel.
  // Because the number or clusters could (conceivably) exceed 10^6, we
  // should not make this a table of ints or floats.  Instead use "ptrdiff_t".
  ptrdiff_t *aiBasinId = NULL;
  ptrdiff_t ***aaaiBasinId = NULL;
  Alloc3D(tomo_in.header.nvoxels,
          &aiBasinId,
          &aaaiBasinId);
  // (Later on we will copy the contents of aaaiBasinId into tomo_out.aaafI
  //  which will be written to a file later.  This way the end-user can
  //  view the results.)

  Watershed(tomo_in.header.nvoxels,
            tomo_in.aaafI,
            aaaiBasinId,
            mask.aaafI,
            settings.watershed_threshold,
            (! settings.clusters_begin_at_maxima),
            settings.neighbor_connectivity,
            settings.watershed_show_boundaries,
            settings.watershed_boundary_label,
            &extrema_crds,
            &extrema_scores,
            &cerr);

  for (int iz = 0; iz < image_size[2]; ++iz)
    for (int iy = 0; iy < image_size[1]; ++iy)
      for (int ix = 0; ix < image_size[0]; ++ix)
        tomo_out.aaafI[iz][iy][ix] = aaaiBasinId[iz][iy][ix];

  Dealloc3D(tomo_in.header.nvoxels,
            &aiBasinId,
            &aaaiBasinId);

  // Did the user supply a mask?
  // Watershed() intentionally does not modify voxels which lie 
  // outside the mask.  These voxels will have random undefined values 
  // unless we assign them manually.  We do this below.
  float UNDEFINED = -1;
  for (int iz = 0; iz < image_size[2]; iz++)
    for (int iy = 0; iy < image_size[1]; iy++)
      for (int ix = 0; ix < image_size[0]; ix++)
        if (mask.aaafI && (mask.aaafI[iz][iy][ix] == 0.0))
          tomo_out.aaafI[iz][iy][ix] = UNDEFINED;
} //HandleWatershed()






void
HandleClusterConnected(Settings settings,
                       MrcSimple &tomo_in,
                       MrcSimple &tomo_out,
                       MrcSimple &mask,
                       float voxel_width[3])
{
  vector<vector<array<int, 3> > > *pMustLinkConstraints = NULL;

  if (settings.must_link_filename != "") {
    // Prepare the list of coordinates in the settings.must_link_constraints
    // array beforehand so that it will be ready by the time we have to
    // invoke ClusterConnected().
    ProcessLinkConstraints(settings.must_link_filename,
                           settings.must_link_constraints,
                           voxel_width);
    if (settings.must_link_constraints.size() > 0)
      pMustLinkConstraints = &settings.must_link_constraints;
  }

  int image_size[3];
  for (int d = 0; d < 3; d++)
    image_size[d] = tomo_in.header.nvoxels[d];

  vector<array<int, 3> > cluster_centers;
  vector<float> cluster_sizes;
  vector<float> cluster_saliencies;

  // Create a temporary array to store the cluster membership for each voxel.
  // Because the number or clusters could (conceivably) exceed 10^6, we
  // should not make this a table of ints or floats.  Instead use "ptrdiff_t".
  ptrdiff_t *aiClusterId = NULL;
  ptrdiff_t ***aaaiClusterId = NULL;
  Alloc3D(tomo_in.header.nvoxels,
          &aiClusterId,
          &aaaiClusterId);
  // (Later on we will copy the contents of aaaiClusterId into tomo_out.aaafI
  //  which will be written to a file later.  This way the end-user can
  //  view the results.)

  ClusterConnected(tomo_in.header.nvoxels, //image size
                   tomo_in.aaafI, //<-saliency
                   aaaiClusterId, //<-which cluster does each voxel belong to?  (results will be stored here)
                   mask.aaafI,
                   settings.connect_threshold_saliency,
                   static_cast<ptrdiff_t>(0), //this value is ignored, but it specifies the type of array we are using
                   true,  //(voxels not belonging to clusters are assigned the highest value = num_clusters+1)
                   static_cast<array<float, 3> ***>(NULL),
                   -std::numeric_limits<float>::infinity(),
                   -std::numeric_limits<float>::infinity(),
                   false, //normal default value for this (ignored) parameter
                   static_cast<float****>(NULL),
                   -std::numeric_limits<float>::infinity(),
                   -std::numeric_limits<float>::infinity(),
                   true,  //normal default value for this (ignored) parameter
                   1,
                   &cluster_centers,
                   &cluster_sizes,
                   &cluster_saliencies,
                   ClusterSortCriteria::SORT_BY_SIZE,
                   static_cast<float***>(NULL),
                   #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
                   static_cast<array<float, 3> ***>(NULL),
                   #endif
                   pMustLinkConstraints,
                   true, //(clusters begin at regions of high saliency)
                   &cerr);  //!< print progress to the user

  // Now, copy the contents of aaaClusterId into tomo_out.aaafI
  for (int iz = 0; iz < image_size[2]; ++iz)
    for (int iy = 0; iy < image_size[1]; ++iy)
      for (int ix = 0; ix < image_size[0]; ++ix)
        tomo_out.aaafI[iz][iy][ix] = aaaiClusterId[iz][iy][ix];

  Dealloc3D(tomo_in.header.nvoxels,
            &aiClusterId,
            &aaaiClusterId);

} //HandleClusterConnected()







void
HandleRidgeDetector(Settings settings,
                    MrcSimple &tomo_in,
                    MrcSimple &tomo_out,
                    MrcSimple &mask,
                    float voxel_width[3])
{
  cerr << "filter_type = surface ridge detector\n";

  vector<vector<array<int, 3> > > *pMustLinkConstraints = NULL;

  if (settings.must_link_filename != "") {
    // Prepare the list of coordinates in the settings.must_link_constraints
    // array beforehand so that it will be ready by the time we have to
    // invoke ClusterConnected().
    ProcessLinkConstraints(settings.must_link_filename,
                           settings.must_link_constraints,
                           voxel_width);
    if (settings.must_link_constraints.size() > 0)
      pMustLinkConstraints = &settings.must_link_constraints;
  }

  int image_size[3];
  for (int d = 0; d < 3; d++)
    image_size[d] = tomo_in.header.nvoxels[d];

  selfadjoint_eigen3::EigenOrderType eival_order;
  if (settings.ridges_are_maxima) {
    eival_order = selfadjoint_eigen3::INCREASING_EIVALS;
    //We want the first eigenvalue to be the smallest.
    //Since, we assume it is negative, equivalently we want the first eigenvalue
    //to be the the most negative eigenvalue (the largest in magnitude).
  }
  else {
    eival_order = selfadjoint_eigen3::DECREASING_EIVALS;
    //Otherwise, we want the first eigenvalue to be the largest.
    //(In this case, it is assumed to be positive.)
  }

  float sigma = settings.width_a[0];

  // Allocate space for the arrays where will will store the image intensity
  // gradient vectors and Hessian matrices during the calculation. This is messy
  // I apologize for not using simpler data structures, but I was trying
  // to reduce memory usage.

  // Allocate the image storing the gradient at every voxel.  This could be
  // implemented as "float ****aaaafGradient;", however instead we use:

  array<float, 3> ***aaaafGradient;

  // Alternatively, we could have defined it using:
  // float ****aaaafGradient;  (or  float ***aaaafGradient[3];)
  // however defining it as "array<float, 3> ***aaaafGradient;" makes it easier
  // to allocate memory for this array using "Alloc3d()" defined in "alloc3d.h".

  array<float, 3> *aafGradient;
  // The "Alloc3d()" function also needs an argument which is a pointer to
  // a 1-D array with the contents of aaaafGradient arranged consecutively.
  // That's what "aafGradient" is.

  // Now use Alloc3D() to allocate space for both aafGradient and aaaafGradient.

  Alloc3D(tomo_in.header.nvoxels,
          &(aafGradient),
          &(aaaafGradient));

  // The storage requirement for Hessians (6 floats) is large enough that
  // I decided to represent hessians using a CompactMultiChannelImage3D.
  // Internally this is a 4-dimensional array, however the last dimension
  // is only allocated (non-NULL) for voxels which were selected by the user
  // (ie voxels for which the mask is non-zero).  This can reduce memory usage
  // by a factor of up to 3 (for floats) for this array.
  CompactMultiChannelImage3D<float> tmp_tensor(6);
  tmp_tensor.Resize(tomo_in.header.nvoxels, mask.aaafI, &cerr);


  // How did the user specify how wide to make the filter window?
  if (settings.filter_truncate_ratio <= 0) {
    assert(settings.filter_truncate_threshold > 0.0);
    //    filter_truncate_threshold = exp(-(1/2)*filter_truncate_ratio^2);
    //    -> filter_truncate_ratio^2 = -2*log(filter_truncate_threshold)
    settings.filter_truncate_ratio = sqrt(-2*log(settings.filter_truncate_threshold));
  }

  MrcSimple tomo_background;
  bool subtract_background = (settings.width_b[0] > 0.0);
  if (subtract_background) {
    tomo_background = tomo_in;

    int truncate_halfwidth = floor(settings.width_b[0] *
                                   settings.filter_truncate_ratio);

    ApplyGauss(tomo_in.header.nvoxels,
               tomo_in.aaafI,
               tomo_background.aaafI,
               mask.aaafI,
               settings.width_b[0],
               truncate_halfwidth,
               true,
               &cerr);

    truncate_halfwidth = floor(settings.width_a[0] *
                               settings.filter_truncate_ratio);

    ApplyGauss(tomo_in.header.nvoxels,
               tomo_in.aaafI,
               tomo_out.aaafI,
               mask.aaafI,
               settings.width_a[0],
               truncate_halfwidth,
               true,
               &cerr);
  }

  // Using this blurred image, calculate the 2nd derivative matrix everywhere:

  cerr << "Applying a Gaussian blur of width sigma="
       << sigma*voxel_width[0] << " voxels" << endl;

  CalcHessian(tomo_in.header.nvoxels,
              tomo_in.aaafI,
              aaaafGradient,
              tmp_tensor.aaaafI,
              mask.aaafI,
              sigma,
              settings.filter_truncate_ratio,
              &cerr);

  cerr << "Diagonalizing the Hessians" << endl;


  // DELETE THIS DEBUGGING CRUFT
  //for(int iz=0; iz < image_size[2]; iz++)
  //  for(int iy=0; iy < image_size[1]; iy++)
  //    for(int ix=0; ix < image_size[0]; ix++)
  //      tomo_out.aaafI[iz][iy][ix] =
  //        FrobeniusNormSqdSym3(tmp_tensor.aaaafI[iz][iy][ix]);
  //return;
  // DELETE THIS DEBUGGING CRUFT


  // We need to store the direction of the most important eigenvector
  // somewhere.  To save space, why not store it in the aaaafGradient
  // array?  (At this point, we are no longer using it).

  array<float, 3> ***aaaafDirection = aaaafGradient;

  // This just effectively changes the name
  // from "aaaafGradient" to "aaaafDirection".
  // The purpose of the name change is to make it easier to read the code later.

  // Optional: store the saliency (score) of each voxel in tomo_out.aaafI
  for(int iz=0; iz < image_size[2]; iz++)
    for(int iy=0; iy < image_size[1]; iy++)
      for(int ix=0; ix < image_size[0]; ix++)
        tomo_out.aaafI[iz][iy][ix] = 0.0;

  for(int iz=0; iz < image_size[2]; iz++) {
    #pragma omp parallel for collapse(2)
    for(int iy=0; iy < image_size[1]; iy++) {
      for(int ix=0; ix < image_size[0]; ix++) {

        if (! tmp_tensor.aaaafI[iz][iy][ix]) //ignore voxels that are either
          continue;                         //in the mask or on the boundary

        float eivals[3];
        float eivects[3][3];

        ConvertFlatSym2Evects3(tmp_tensor.aaaafI[iz][iy][ix],
                               eivals,
                               eivects,
                               eival_order);

        float score;

        // DEBUG: REMOVE THE NEXT IF STATMENT AFTER DEBUGGING IS FINISHED
        #ifndef NDEBUG
        if ((ix==image_size[0]/2) && //if ((ix==78) &&
            (iy == image_size[1] / 2) &&
            (iz == image_size[2] / 2))
        {
          cerr << "[iz][iy][ix]=["<<iz<<"]["<<iy<<"]["<<ix<<"]\n"
               << "eivals = "<<eivals[0]<<","<<eivals[1]<<","<<eivals[2]<<"\n"
               << "eivects = \n"
               << "    "<<eivects[0][0]<<","<<eivects[0][1]<<","<<eivects[0][2]<<"\n"
               << "    "<<eivects[1][0]<<","<<eivects[1][1]<<","<<eivects[1][2]<<"\n"
               << "    "<<eivects[2][0]<<","<<eivects[2][1]<<","<<eivects[2][2]<<"\n"
               << endl;
        }
        #endif  //#ifndef NDEBUG


        score = ScoreHessianPlanar(eivals,
                                   aaaafGradient[iz][iy][ix]);

        float peak_height = 1.0;
        if (tomo_background.aaafI)
          peak_height = (tomo_in.aaafI[iz][iy][ix] -
                         tomo_background.aaafI[iz][iy][ix]);
        score *= peak_height;

        tomo_out.aaafI[iz][iy][ix] = score;

        #if 0
        //I will use this code eventually, but not yet
        float gradient_along_v1 = DotProduct3(grad, eivects[0]);
        float distance_to_ridge;
        if (lambda1 != 0)
          distance_to_ridge = abs(gradient_along_v1 / lambda1);
        else
          distance_to_ridge = std::numeric_limits<float>::infinity();

        bool ridge_located_in_same_voxel = true;
        for (int d=0; d<3; d++) {
          float ridge_voxel_location = eivects[0][d] * distance_to_ridge;
          if (abs(ridge_voxel_location) > 0.5)
            ridge_located_in_same_voxel = false;
        }
        //if ((distance_to_ridge < settings.ridge_detector_search_width*0.5) &&

        if (! ridge_located_in_same_voxel)
          out_tomo.aaafI[iz][iy][ix] = 0.0;  //then ignore discard this voxel
        #endif


        aaaafDirection[iz][iy][ix][0] = eivects[0][0];
        aaaafDirection[iz][iy][ix][1] = eivects[0][1];
        aaaafDirection[iz][iy][ix][2] = eivects[0][2];


      } //for(int ix=0; ix < image_size[0]; ix++)
    } //for(int iy=0; iy < image_size[1]; iy++)
  } //for(int iz=0; iz < image_size[2]; iz++)




  { // Use thresholding to reduce the number of voxels that we have to consider
    
    float hessian_score_threshold = settings.surface_hessian_score_threshold;
    if (settings.surface_hessian_score_threshold_is_a_fraction) {
      float hessian_score_fraction = settings.surface_hessian_score_threshold;
      size_t n_voxels = 0;
      for (int iz=0; iz < image_size[2]; iz++) {
        for (int iy=0; iy < image_size[1]; iy++) {
          for (int ix=0; ix < image_size[0]; ix++) {
            if (mask.aaafI && (mask.aaafI[iz][iy][ix] == 0))
              continue;
            n_voxels++;
          }
        }
      }
      cerr << " -- sorting ALL voxels by ridge saliency --\n" << endl;
      vector<float> saliencies(n_voxels);
      size_t i = 0;
      for (int iz=0; iz < image_size[2]; iz++) {
        for (int iy=0; iy < image_size[1]; iy++) {
          for (int ix=0; ix < image_size[0]; ix++) {
            if (mask.aaafI && (mask.aaafI[iz][iy][ix] == 0))
              continue;
            saliencies[i] = tomo_out.aaafI[iz][iy][ix];
            i++;
          }
        }
      }
      sort(saliencies.rbegin(),
           saliencies.rend());
      i = floor(n_voxels * hessian_score_fraction);
      hessian_score_threshold = saliencies[i];
    }

    // Now that we know what the threshold is (for not discarding)
    // discard voxels whose saliencies fall below this number.
    // (Do this by setting their saliency to 0.
    //  This will cause them to be ignored in later stages of the calculation.)
    for (int iz=0; iz < image_size[2]; iz++) {
      for (int iy=0; iy < image_size[1]; iy++) {
        for (int ix=0; ix < image_size[0]; ix++) {
          if (tomo_out.aaafI[iz][iy][ix] < hessian_score_threshold)
            tomo_out.aaafI[iz][iy][ix] = 0.0;
        }
      }
    }
  } // Use thresholding to reduce the number of voxels that we have to consider



  if (settings.surface_tv_sigma > 0.0) {
    assert(settings.filter_truncate_ratio > 0);

    TV3D<float, int, array<float,3>, float* >
      tv(settings.surface_tv_sigma,
         settings.surface_tv_exponent,
         settings.surface_tv_truncate_ratio);

    tv.TVDenseStick(tomo_in.header.nvoxels,
                    tomo_out.aaafI,
                    aaaafDirection,
                    tmp_tensor.aaaafI,
                    mask.aaafI,
                    mask.aaafI,
                    false,  // (we want to detect surfaces not curves)
                    //settings.surface_hessian_score_threshold,
                    true,   // (do normalize near rectangular image bounaries)
                    false,  // (diagonalize each tensor afterwards?)
                    &cerr);

    for(int iz=0; iz < image_size[2]; iz++) {
      for(int iy=0; iy < image_size[1]; iy++) {
        for(int ix=0; ix < image_size[0]; ix++) {
          float diagonalized_hessian[6];
          DiagonalizeFlatSym3(tmp_tensor.aaaafI[iz][iy][ix],
                              diagonalized_hessian,
                              eival_order);
          float score = ScoreTensorPlanar(diagonalized_hessian);
          float peak_height = 1.0;
          if (tomo_background.aaafI)
            peak_height = (tomo_in.aaafI[iz][iy][ix] -
                           tomo_background.aaafI[iz][iy][ix]);
          score *= peak_height;
          tomo_out.aaafI[iz][iy][ix] = score;
        }
      }
    }
  } // if (settings.surface_tv_sigma > 0.0)




  if (settings.cluster_connected_voxels)
  {
    tomo_in = tomo_out; // (horrible hack.  I should not modify tomo_in.)
                        //  allocate a new 3D array to store the saliency)

    // Copy the principal eigenvector of tmp_tensor into aaaafDirection
    for(int iz=0; iz < image_size[2]; iz++) {
      for(int iy=0; iy < image_size[1]; iy++) {
        for(int ix=0; ix < image_size[0]; ix++) {
          float eivals[3];
          float eivects[3][3];
          ConvertFlatSym2Evects3(tmp_tensor.aaaafI[iz][iy][ix],
                                 eivals,
                                 eivects,
                                 eival_order);
          for (int d=0; d<3; d++)
            aaaafDirection[iz][iy][ix][d] = eivects[0][d];
        }
      }
    }

    vector<array<int, 3> > cluster_centers;
    vector<float> cluster_sizes;
    vector<float> cluster_saliencies;

    // Create a temporary array to store the cluster membership for each voxel.
    // Because the number or clusters could (conceivably) exceed 10^6, we
    // should not make this a table of ints or floats.  Instead use "ptrdiff_t".
    ptrdiff_t *aiClusterId = NULL;
    ptrdiff_t ***aaaiClusterId = NULL;
    Alloc3D(tomo_in.header.nvoxels,
            &aiClusterId,
            &aaaiClusterId);
    // (Later on we will copy the contents of aaaiClusterId into tomo_out.aaafI
    //  which will be written to a file later.  This way the end-user can
    //  view the results.)

    ClusterConnected(tomo_in.header.nvoxels, //image size
                     tomo_in.aaafI, //<-saliency
                     aaaiClusterId, //<-which cluster does each voxel belong to?  (results will be stored here)
                     mask.aaafI,
                     settings.connect_threshold_saliency,
                     static_cast<ptrdiff_t>(0), //this value is ignored, but it specifies the type of array we are using
                     true,  //(voxels not belonging to clusters are assigned the highest value = num_clusters+1)
                     aaaafDirection,
                     settings.connect_threshold_vector_saliency,
                     settings.connect_threshold_vector_neighbor,
                     false, //eigenvector signs are arbitrary so ignore them
                     tmp_tensor.aaaafI,
                     settings.connect_threshold_tensor_saliency,
                     settings.connect_threshold_tensor_neighbor,
                     true,  //the tensor should be positive definite near the target
                     1,
                     &cluster_centers,
                     &cluster_sizes,
                     &cluster_saliencies,
                     ClusterSortCriteria::SORT_BY_SIZE,
                     static_cast<float***>(NULL),
                     #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
                     aaaafDirection,
                     #endif
                     pMustLinkConstraints,
                     true, //(clusters begin at regions of high saliency)
                     &cerr);  //!< print progress to the user

    // Now, copy the contents of aaaClusterId into tomo_out.aaafI
    for (int iz = 0; iz < image_size[2]; ++iz)
      for (int iy = 0; iy < image_size[1]; ++iy)
        for (int ix = 0; ix < image_size[0]; ++ix)
          tomo_out.aaafI[iz][iy][ix] = aaaiClusterId[iz][iy][ix];

    Dealloc3D(tomo_in.header.nvoxels,
              &aiClusterId,
              &aaaiClusterId);

  } // if (settings.cluster_connected_voxels)



  //Did the user ask us to generate output files containing surface orientation?
  if (settings.out_normals_fname != "")
  {

    // define the set of voxels we will use
    float ***aaafSelectedVoxels;

    // where is the lookup table to indicate which cluster a voxel belongs to?
    float ***aaafVoxel2Cluster = NULL;
    if (settings.cluster_connected_voxels)
      aaafVoxel2Cluster = tomo_out.aaafI; //tomo_out set by ClusterConnected()

    // (horrible hack.  I should not modify tomo_in.
    aaafSelectedVoxels = tomo_in.aaafI;
    //  Instead I should allocate a new 3D array to store the selected voxels.

    for (int iz = 0; iz < image_size[2]; ++iz) {
      for (int iy = 0; iy < image_size[1]; ++iy) {
        for (int ix = 0; ix < image_size[0]; ++ix) {
          aaafSelectedVoxels[iz][iy][ix] = 0.0;
          if (mask.aaafI && (mask.aaafI[iz][iy][ix] == 0.0))
            continue;
          if (! settings.cluster_connected_voxels)
            aaafSelectedVoxels[iz][iy][ix] = 1.0;
          else {
            assert(aaafVoxel2Cluster);
            if (settings.select_cluster == aaafVoxel2Cluster[iz][iy][ix])
              aaafSelectedVoxels[iz][iy][ix] = 1.0;
          }
        }
      }
    }

    WriteOrientedPointCloud(settings.out_normals_fname,
                            image_size,
                            aaaafDirection,
                            aaafSelectedVoxels,
                            voxel_width);

  } //if (settings.out_normals_fname != "")


  ////Did the user ask us to generate any output files?
  //if ((settings.out_normals_fname == "") &&
  //    (settings.out_file_name != ""))
  //  settings.out_normals_fname =
  //    settings.out_file_name + string(".bnpts");


  Dealloc3D(tomo_in.header.nvoxels,
            &(aafGradient),
            &(aaaafGradient));

} //HandleRidgeDetector()



