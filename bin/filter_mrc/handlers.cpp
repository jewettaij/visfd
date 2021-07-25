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
using namespace visfd;

#include <threshold.hpp>
#include <mrc_simple.hpp>
#include <random_gen.h>
#include <feature.hpp>
#include <segmentation.hpp>
#include <clustering.hpp>
#include <morphology.hpp>
#include <draw.hpp>
#include "err.hpp"
#include "settings.hpp"
#include "file_io.hpp"
#include "filter3d_variants.hpp"
#include "feature_variants.hpp"
#include "feature_unsupported.hpp"
#include "handlers.hpp"
#include "handlers_unsupported.hpp"


void
HandleDilation(const Settings &settings,
               MrcSimple &tomo_in,
               MrcSimple &tomo_out,
               MrcSimple &mask,
               float voxel_width[3])
{
  DilateSphere(settings.morphology_r,
               tomo_in.header.nvoxels,
               tomo_in.aaafI,
               tomo_out.aaafI,
               mask.aaafI,
               settings.morphology_rmax,
               settings.morphology_bmax,
               &cerr);
}


void
HandleErosion(const Settings &settings,
              MrcSimple &tomo_in,
              MrcSimple &tomo_out,
              MrcSimple &mask,
              float voxel_width[3])
{
  ErodeSphere(settings.morphology_r,
              tomo_in.header.nvoxels,
              tomo_in.aaafI,
              tomo_out.aaafI,
              mask.aaafI,
              settings.morphology_rmax,
              settings.morphology_bmax,
              &cerr);
}


void
HandleOpening(const Settings &settings,
              MrcSimple &tomo_in,
              MrcSimple &tomo_out,
              MrcSimple &mask,
              float voxel_width[3])
{
  OpenSphere(settings.morphology_r,
             tomo_in.header.nvoxels,
             tomo_in.aaafI,
             tomo_out.aaafI,
             mask.aaafI,
             settings.morphology_rmax,
             settings.morphology_bmax,
             &cerr);
}


void
HandleClosing(const Settings &settings,
              MrcSimple &tomo_in,
              MrcSimple &tomo_out,
              MrcSimple &mask,
              float voxel_width[3])
{
  CloseSphere(settings.morphology_r,
              tomo_in.header.nvoxels,
              tomo_in.aaafI,
              tomo_out.aaafI,
              mask.aaafI,
              settings.morphology_rmax,
              settings.morphology_bmax,
              &cerr);
}


void
HandleTopHatWhite(const Settings &settings,
                  MrcSimple &tomo_in,
                  MrcSimple &tomo_out,
                  MrcSimple &mask,
                  float voxel_width[3])
{
  WhiteTopHatSphere(settings.morphology_r,
                    tomo_in.header.nvoxels,
                    tomo_in.aaafI,
                    tomo_out.aaafI,
                    mask.aaafI,
                    settings.morphology_rmax,
                    settings.morphology_bmax,
                    &cerr);
}


void
HandleTopHatBlack(const Settings &settings,
                  MrcSimple &tomo_in,
                  MrcSimple &tomo_out,
                  MrcSimple &mask,
                  float voxel_width[3])
{
  BlackTopHatSphere(settings.morphology_r,
                    tomo_in.header.nvoxels,
                    tomo_in.aaafI,
                    tomo_out.aaafI,
                    mask.aaafI,
                    settings.morphology_rmax,
                    settings.morphology_bmax,
                    &cerr);
}



#ifndef CXX17_UNSUPPORTED
void
HandleMedian(const Settings &settings,
             MrcSimple &tomo_in,
             MrcSimple &tomo_out,
             MrcSimple &mask,
             float voxel_width[3])
{
  MedianSphere(settings.median_radius,
               tomo_in.header.nvoxels,
               tomo_in.aaafI,
               tomo_out.aaafI,
               mask.aaafI,
               &cerr);
}
#endif



void
HandleGGauss(const Settings &settings,
             MrcSimple &tomo_in,
             MrcSimple &tomo_out,
             MrcSimple &mask,
             float voxel_width[3])
{
  float A; // scaling coefficient (computed by normalization)

  Filter3D<float, int>
    filter = GenFilterGenGauss3D(settings.width_a,
                                 settings.m_exp,
                                 settings.filter_truncate_ratio,
                                 settings.filter_truncate_threshold,
                                 &A);

  filter.Apply(tomo_in.header.nvoxels,
               tomo_in.aaafI,
               tomo_out.aaafI,
               mask.aaafI,
               settings.normalize_near_boundaries,
               &cerr);

  cerr << " Filter Used:\n"
      " h(x,y,z)   = A*exp(-((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2)^(m/2))\n"
      "  ... where      A = " << A << "\n"
      "                 m = " << settings.m_exp << "\n" 
      "   (a_x, a_y, a_z) = "
       << "(" << settings.width_a[0]
       << " " << settings.width_a[1]
       << " " << settings.width_a[2] << ")\n";
    cerr << " You can plot a slice of this function\n"
         << "     in the X direction using:\n"
      " draw_filter_1D.py -ggauss " << A
         << " " << settings.width_a[0]
         << " " << settings.m_exp << endl;
    if ((settings.width_a[1] != settings.width_a[0]) ||
        (settings.width_a[2] != settings.width_a[0])) {
      cerr << " and in the Y direction using:\n"
        " draw_filter_1D.py -ggauss " << A
           << " " << settings.width_a[1]
           << " " << settings.m_exp << endl;
      cerr << " and in the Z direction using:\n"
        " draw_filter_1D.py -ggauss " << A
           << " " << settings.width_a[2]
           << " " << settings.m_exp << endl;
    }
} //HandleGGauss()




void
HandleGauss(const Settings &settings,
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
                 settings.normalize_near_boundaries,
                 &cerr);

  cerr << " Filter Used: A discrete Gaussian kernel, approximately equal to\n"
    " h(x,y,z)   ≈ A*exp(-0.5*((x/σ_x)^2 + (y/σ_y)^2 + (z/σ_z)^2))\n"
    " ... where  A = " << A << "\n"
    "   (σ_x, σ_y, σ_z) = "
       << "(" << settings.width_a[0]
       << " " << settings.width_a[1]
       << " " << settings.width_a[2] << ")\n";
  cerr << " You can plot a slice of this function\n"
       << "     in the X direction using:\n"
       << " draw_filter_1D.py -gauss "
       << A << " " << settings.width_a[0] << endl;
  if ((settings.width_a[1] != settings.width_a[0]) ||
      (settings.width_a[2] != settings.width_a[0])) {
    cerr << " and in the Y direction using:\n"
         << " draw_filter_1D.py -gauss "
         << A << " " << settings.width_a[1] << endl;
    cerr << " and in the Z direction using:\n"
         << " draw_filter_1D.py -gauss "
         << A << " " << settings.width_a[2] << endl;
  }
} //HandleGauss()




void
HandleDogg(const Settings &settings,
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
HandleDog(const Settings &settings,
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

  cerr << " Filter Used: The difference of two discrete Gaussian kernels ≈\n"
    " h(x,y,z)   = h_a(x,y,z) - h_b(x,y,z)\n"
    " h_a(x,y,z) ≈ A*exp(-0.5*((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))  (approximately)\n"
    " h_b(x,y,z) ≈ B*exp(-0.5*((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))  (approximately)\n"
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

  if ((settings.width_a[1] != settings.width_a[0]) ||
      (settings.width_a[2] != settings.width_a[0]) ||
      (settings.width_b[1] != settings.width_b[0]) ||
      (settings.width_b[2] != settings.width_b[0])) {
    cerr << " and in the Y direction using:\n"
      " draw_filter_1D.py -dog " << A << " " << B
         << " " << settings.width_a[1] << " " << settings.width_b[1] << endl;
    cerr << " and in the Z direction using:\n"
      " draw_filter_1D.py -dog " << A << " " << B
         << " " << settings.width_a[2] << " " << settings.width_b[2] << endl;
  }

} //HandleDog()





void
HandleLoGDoG(const Settings &settings,
             MrcSimple &tomo_in,
             MrcSimple &tomo_out,
             MrcSimple &mask,
             float voxel_width[3])
{
  cerr << "filter_type = Laplacian of Gaussians (LoG)\n"
       << "  (This will be approximated as a Difference of Gaussians,\n"
       << "   as explained below.)\n";

  float A, B;

  ApplyLog(tomo_in.header.nvoxels,
           tomo_in.aaafI,
           tomo_out.aaafI,
           mask.aaafI,
           settings.log_width,
           settings.delta_sigma_over_sigma,
           settings.filter_truncate_ratio,
           settings.filter_truncate_threshold,
           &A,
           &B,
           &cerr);

  cerr << " Filter Used: The difference of two discrete Gaussian kernels ≈\n"
    " h(x,y,z)   = h_a(x,y,z) - h_b(x,y,z)\n"
    " h_a(x,y,z) ≈ A*exp(-0.5*((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))  (approximately)\n"
    " h_b(x,y,z) ≈ B*exp(-0.5*((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))  (approximately)\n"
    "  ... where      A = " << A << "\n"
    "                 B = " << B << "\n" 
    "   (a_x, a_y, a_z) = "
       << "(" << settings.log_width[0]*(1-0.5*settings.delta_sigma_over_sigma)
       << " " << settings.log_width[1]*(1-0.5*settings.delta_sigma_over_sigma)
       << " " << settings.log_width[2]*(1-0.5*settings.delta_sigma_over_sigma)
       << ")\n"
    "   (b_x, b_y, b_z) = "
       << "(" << settings.log_width[0]*(1+0.5*settings.delta_sigma_over_sigma)
       << " " << settings.log_width[1]*(1+0.5*settings.delta_sigma_over_sigma)
       << " " << settings.log_width[2]*(1+0.5*settings.delta_sigma_over_sigma)
       << ")\n";
  cerr << " You can plot a slice of this function\n"
       << "     in the X direction using:\n"
    " draw_filter_1D.py -dog " << A << " " << B
       << " " << settings.log_width[0]*(1-0.5*settings.delta_sigma_over_sigma)
       << " " << settings.log_width[0]*(1+0.5*settings.delta_sigma_over_sigma)
       << endl;
  if ((settings.log_width[1] != settings.log_width[0]) ||
      (settings.log_width[2] != settings.log_width[0])) {
    cerr << " and in the Y direction using:\n"
      " draw_filter_1D.py -dog " << A << " " << B
         << " " << settings.log_width[1]*(1-0.5*settings.delta_sigma_over_sigma)
         << " " << settings.log_width[1]*(1+0.5*settings.delta_sigma_over_sigma)
         << endl;
    cerr << " and in the Z direction using:\n"
      " draw_filter_1D.py -dog " << A << " " << B
         << " " << settings.log_width[2]*(1-0.5*settings.delta_sigma_over_sigma)
         << " " << settings.log_width[2]*(1+0.5*settings.delta_sigma_over_sigma)
         << endl;
  }
} // HandleLoGDoG()








void
HandleBlobsNonmaxSuppression(const Settings &settings,
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
                     settings.sphere_decals_foreground,
                     settings.sphere_decals_scale);



  cerr << " --- discarding blobs in file \n"
       << "     \"" << settings.in_coords_file_name << "\" ---\n"
       << "\n";


  // Discard blobs based on score or size?

  if (((settings.score_lower_bound !=
        -std::numeric_limits<float>::infinity()) ||
       (settings.score_upper_bound !=
        std::numeric_limits<float>::infinity()))
      ||
      ((settings.sphere_diameters_lower_bound !=
        -std::numeric_limits<float>::infinity()) ||
       (settings.sphere_diameters_upper_bound !=
        std::numeric_limits<float>::infinity())))
  {
    // report to the user
    if ((settings.score_lower_bound !=
         -std::numeric_limits<float>::infinity()) &&
        (settings.score_upper_bound !=
         std::numeric_limits<float>::infinity()))
      cerr << "  discarding blobs based on score" << endl;
    if ((settings.sphere_diameters_lower_bound !=
         -std::numeric_limits<float>::infinity()) &&
        (settings.sphere_diameters_upper_bound !=
         std::numeric_limits<float>::infinity()))
      cerr << "  discarding blobs based on size (diameter)" << endl;

    // If the user supplied explicit minimum and maximum threshold scores or
    // diameters for blobs, discard the blobs that lie outside this range.

    vector<array<float,3> > crds_cpy;
    vector<float> diameters_cpy;
    vector<float> scores_cpy;
    for (int i = 0; i < crds.size(); i++) {
      if ((scores[i] >= settings.score_lower_bound) &&
          (scores[i] <= settings.score_upper_bound) &&
          (diameters[i] >= settings.sphere_diameters_lower_bound) &&
          (diameters[i] <= settings.sphere_diameters_upper_bound))
      {
        crds_cpy.push_back(crds[i]);
        diameters_cpy.push_back(diameters[i]);
        scores_cpy.push_back(scores[i]);
      }
    }
    crds = crds_cpy;
    diameters = diameters_cpy;
    scores = scores_cpy;

  } // else clause for "if (settings.auto_thresh_score && ..."


  // Discard blobs outside the mask?
  if ((crds.size() > 0) && (mask.aaafI != nullptr)) {
    cerr << "  discarding blobs outside the mask" << endl;
    DiscardMaskedBlobs(crds,
                       diameters,
                       scores,
                       mask.aaafI);
  }

  // Discard overlapping blobs?
  cerr << "  discarding overlapping blobs" << endl;
  DiscardOverlappingBlobs(crds,
                          diameters, 
                          scores,
                          settings.nonmax_min_radial_separation_ratio,
                          settings.nonmax_max_volume_overlap_large,
                          settings.nonmax_max_volume_overlap_small,
                          SORT_DECREASING_MAGNITUDE,
                          &cerr);

  cerr << " " << crds.size() << " blobs remaining" << endl;

  float score_lower_bound = settings.score_lower_bound;
  float score_upper_bound = settings.score_upper_bound;

  // Finally, use supervised learning to discard the remaining blobs?
  // (Always do this AFTER discarding overlapping blobs.)
  if (settings.auto_thresh_score &&
      (settings.training_pos_crds.size() > 0) &&
      (settings.training_neg_crds.size() > 0))
  {

    cerr << "  discarding blobs based on score using training data" << endl;

    
    DiscardBlobsByScoreSupervised(crds,
                                  diameters,
                                  scores,
                                  settings.training_pos_crds,
                                  settings.training_neg_crds,
                                  SORT_DECREASING_MAGNITUDE,
                                  &score_lower_bound,
                                  &score_upper_bound,
                                  &cerr);

    cerr << " " << crds.size() << " blobs remaining" << endl;
  }

  // Write out the remaining blobs to a file:
  if (settings.out_coords_file_name != "") {
    fstream out_coords_file;
    out_coords_file.open(settings.out_coords_file_name, ios::out);
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
HandleBlobScoreSupervisedMulti(const Settings &settings,
                               float voxel_width[3])
{
  float threshold_lower_bound;
  float threshold_upper_bound;

  int Nsets = settings.multi_in_coords_file_names.size();

  vector<vector<array<float,3> > > blob_crds_multi(Nsets);
  vector<vector<float> > blob_diameters_multi(Nsets);
  vector<vector<float> > blob_scores_multi(Nsets);

  for (int I = 0; I < settings.multi_in_coords_file_names.size(); ++I) {
    ReadBlobCoordsFile(settings.multi_in_coords_file_names[I],
                       &blob_crds_multi[I],
                       &blob_diameters_multi[I],
                       &blob_scores_multi[I],
                       voxel_width[0],
                       settings.sphere_decals_diameter,
                       settings.sphere_decals_foreground,
                       settings.sphere_decals_scale);
  }

  assert(Nsets == blob_crds_multi.size());
  assert(Nsets == blob_diameters_multi.size());
  assert(Nsets == blob_scores_multi.size());
  assert(Nsets == settings.multi_training_pos_crds.size());
  assert(Nsets == settings.multi_training_neg_crds.size());

  ChooseBlobScoreThresholdsMulti(blob_crds_multi,
                                 blob_diameters_multi,
                                 blob_scores_multi,
                                 settings.multi_training_pos_crds,
                                 settings.multi_training_neg_crds,
                                 &threshold_lower_bound,
                                 &threshold_upper_bound,
                                 SORT_DECREASING_MAGNITUDE,
                                 &cerr);

} //HandleBlobScoreSupervisedMulti()


  

void
HandleDrawSpheres(const Settings &settings,
                  MrcSimple &tomo_in,
                  MrcSimple &tomo_out,
                  MrcSimple &mask,
                  float voxel_width[3])
{

  vector<array<float,3> > crds;
  vector<float> diameters;
  vector<float> scores;

  MrcSimple unallocated_image;
  HandleBlobsNonmaxSuppression(settings,
                               //mask, <-- commenting out.  keep masked blobs
                               unallocated_image,//use a NULL image instead
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

  DrawSpheres(tomo_out.header.nvoxels,
              tomo_out.aaafI,
              mask.aaafI,
              crds,
              &diameters,
              &shell_thicknesses,
              &scores,
              tomo_in.aaafI,
              settings.sphere_decals_background,
              settings.sphere_decals_background_scale,
              settings.sphere_decals_background_norm,
              settings.sphere_decals_foreground_norm);

} //HandleDrawSpheres()






void
HandleBlobDetector(const Settings &settings,
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

  aaaafI[0] = Alloc3D<float>(tomo_in.header.nvoxels);

  if (tomo_out.aaafI) {
    // Optional: Instead of allocating aaaafI[1] and aafI[1], borrow the memory
    // you've already allocated to tomo_out.aaafI and afI. (Goal: Save memory.)
    aaaafI[1] = tomo_out.aaafI;
  }
  else
    aaaafI[1] = Alloc3D<float>(tomo_in.header.nvoxels);

  aaaafI[2] = Alloc3D<float>(tomo_in.header.nvoxels);

  // This way we can save memory and also save the 
  // filtered image to a file which we can view using IMOD.

  _BlobDogNM(tomo_in.header.nvoxels,
             tomo_in.aaafI,
             mask.aaafI,
             settings.blob_diameters,  // try detecting blobs of these diameters
             &minima_crds_voxels,  // store minima x,y,z coords here
             &maxima_crds_voxels,  // store maxima x,y,z coords here
             &minima_diameters, // corresponding diameter for that minima
             &maxima_diameters, // corresponding diameter for that maxima
             &minima_scores, // what was the blob's score?
             &maxima_scores, // ("score" = intensity after filtering)
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
             aaaafI);


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
              true,
              false,
              nullptr);

    fstream minima_file;
    minima_file.open(settings.blob_minima_file_name, ios::out);
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
              false,
              false,
              nullptr);

    fstream maxima_file;
    maxima_file.open(settings.blob_maxima_file_name, ios::out);
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

    DrawSpheres(tomo_out.header.nvoxels,
                tomo_out.aaafI,
                mask.aaafI,
                display_crds_voxels,
                &display_diameters,
                &display_shell_thicknesses,
                &display_scores,
                tomo_in.aaafI,
                settings.sphere_decals_background,
                settings.sphere_decals_background_scale,
                settings.sphere_decals_background_norm,
                false);


  } //if (tomo_out.aaafI)



  // Clean up.  Deallocate temporary arrays.
  Dealloc3D(aaaafI[0]);

  // Optional: Instead of allocating aaaafI[1] and aafI[1], borrow the memory
  // you've already allocated to tomo_out.aaafI and afI. (Goal: Save memory.)
  if (! tomo_out.aaafI)  // <-- Was tomo_out.aaafI available?
    // However, if we didn't allocate it this memory, then don't delete it:
    Dealloc3D(aaaafI[1]);

  Dealloc3D(aaaafI[2]);

} //HandleBlobDetector()






void
HandleThresholds(const Settings &settings,
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

    float in_threshold_01_a = ave_intensity +
      settings.in_threshold_01_a * stddev_intensity;
    float in_threshold_01_b = ave_intensity +
      settings.in_threshold_01_b * stddev_intensity;

    cerr << "ave="<< ave_intensity <<", stddev="<<stddev_intensity << endl;
    cerr << "  Clipping intensities between ["
         << in_threshold_01_a << ", "
         << in_threshold_01_b << "]" << endl;
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
          if (settings.in_threshold_01_a == settings.in_threshold_01_b)
            tomo_out.aaafI[iz][iy][ix] = ((tomo_out.aaafI[iz][iy][ix] >
                                           settings.in_threshold_01_a)
                                          ? settings.out_thresh_b_value
                                          : settings.out_thresh_a_value);
          else
            tomo_out.aaafI[iz][iy][ix] =
              Threshold2(tomo_out.aaafI[iz][iy][ix],
                         settings.in_threshold_01_a,
                         settings.in_threshold_01_b,
                         (settings.out_thresh2_use_clipping
                          ? settings.in_threshold_01_a
                          : settings.out_thresh_a_value),
                         (settings.out_thresh2_use_clipping
                          ? settings.in_threshold_01_b
                          : settings.out_thresh_b_value));
        }
        else
          tomo_out.aaafI[iz][iy][ix] =
            Threshold4(tomo_out.aaafI[iz][iy][ix],
                       settings.in_threshold_01_a,
                       settings.in_threshold_01_b,
                       settings.in_threshold_10_a,
                       settings.in_threshold_10_b,
                       settings.out_thresh_a_value,
                       settings.out_thresh_b_value);
      }
    }
  }
} //HandleThresholds()



void
HandleExtrema(const Settings &settings,
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
  vector<array<float, 3> > *pv_minima_crds_voxels = nullptr;
  vector<array<float, 3> > *pv_maxima_crds_voxels = nullptr;
  vector<float> *pv_minima_scores = nullptr;
  vector<float> *pv_maxima_scores = nullptr;
  vector<size_t> *pv_minima_nvoxels = nullptr;
  vector<size_t> *pv_maxima_nvoxels = nullptr;


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

  if ((minima_crds_voxels.size() > 0) && (mask.aaafI != nullptr))
    DiscardMaskedBlobs(minima_crds_voxels,
                       minima_diameters,
                       minima_scores,
                       mask.aaafI);

  if ((maxima_crds_voxels.size() > 0) && (mask.aaafI != nullptr))
    DiscardMaskedBlobs(maxima_crds_voxels,
                       maxima_diameters,
                       maxima_scores,
                       mask.aaafI);

  if ((settings.nonmax_min_radial_separation_ratio > 0.0) ||
      (settings.nonmax_max_volume_overlap_large < 1.0) ||
      (settings.nonmax_max_volume_overlap_small < 1.0)) {
    if ((settings.sphere_decals_diameter > 0) &&
        (settings.nonmax_min_radial_separation_ratio > 0.0)) {
      DiscardOverlappingBlobs(minima_crds_voxels,
                              minima_diameters, 
                              minima_scores,
                              settings.nonmax_min_radial_separation_ratio,
                              settings.nonmax_max_volume_overlap_large,
                              settings.nonmax_max_volume_overlap_small,
                              SORT_INCREASING,
                              &cerr);
    }

    if ((settings.sphere_decals_diameter > 0) &&
        (settings.nonmax_min_radial_separation_ratio > 0.0)) {
      DiscardOverlappingBlobs(maxima_crds_voxels,
                              maxima_diameters, 
                              maxima_scores,
                              settings.nonmax_min_radial_separation_ratio,
                              settings.nonmax_max_volume_overlap_large,
                              settings.nonmax_max_volume_overlap_small,
                              SORT_DECREASING,
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
    minima_file.open(settings.find_minima_file_name, ios::out);
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
    coords_file.open(settings.find_maxima_file_name, ios::out);
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
HandleLocalFluctuations(const Settings &settings,
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
                            settings.normalize_near_boundaries,
                            &cerr);
}







void
HandleWatershed(const Settings &settings,
                MrcSimple &tomo_in,
                MrcSimple &tomo_out,
                MrcSimple &mask,
                float voxel_width[3])
{
  int image_size[3];
  for (int d = 0; d < 3; d++)
    image_size[d] = tomo_in.header.nvoxels[d];

  vector<array<float, 3> > extrema_crds;
  vector<float> extrema_scores;

  // Create a temporary array to store the basin membership for each voxel.
  // Because the number or clusters could (conceivably) exceed 10^6, we
  // should not make this a table of ints or floats.  Instead use "ptrdiff_t".
  ptrdiff_t ***aaaiBasinId = Alloc3D<ptrdiff_t>(tomo_in.header.nvoxels);

  // (Later on we will copy the contents of aaaiBasinId into tomo_out.aaafI
  //  which will be written to a file later.  This way the end-user can
  //  view the results.)

  
  ptrdiff_t undefined_voxel_brightness = settings.undefined_voxel_brightness;
  if (settings.undefined_voxels_are_max)
    undefined_voxel_brightness = -1;


  //** DEBUGGING THE aaaiMarkers FEATURE.
  //** REMOVE THIS CODE EVENTUALLY.
  //** ptrdiff_t ***aaaiDebugMarkers = Alloc3D<ptrdiff_t>(image_size);
  //** for (int iz = 0; iz < image_size[2]; iz++)
  //**   for (int iy = 0; iy < image_size[1]; iy++)
  //**     for (int ix = 0; ix < image_size[0]; ix++)
  //**       aaaiDebugMarkers[iz][iy][ix] = 0;
  //** aaaiDebugMarkers[40][58][70] = 7;
  //** aaaiDebugMarkers[40][57][70] = 18;
  //** aaaiDebugMarkers[40][42][48] = 13;
  //** END OF: DEBUGGING aaaiDebugMarkers FEATURE


  size_t num_basins =
    Watershed(tomo_in.header.nvoxels,
              tomo_in.aaafI,
              aaaiBasinId,
              mask.aaafI,
              static_cast<ptrdiff_t ***>(nullptr),
              settings.watershed_threshold,
              (! settings.clusters_begin_at_maxima),
              settings.neighbor_connectivity,
              settings.watershed_show_boundaries,
              static_cast<ptrdiff_t>(settings.watershed_boundary_label),
              static_cast<ptrdiff_t>(-1), // an impossible value
              &extrema_crds,
              &extrema_scores,
              &cerr);

  ptrdiff_t max_label = aaaiBasinId[0][0][0];
  for (int iz = 0; iz < image_size[2]; ++iz)
    for (int iy = 0; iy < image_size[1]; ++iy)
      for (int ix = 0; ix < image_size[0]; ++ix)
        max_label = std::max(max_label, aaaiBasinId[iz][iy][ix]);

  // Now, copy the contents of aaaClusterId into tomo_out.aaafI
  for (int iz = 0; iz < image_size[2]; ++iz) {
    for (int iy = 0; iy < image_size[1]; ++iy) {
      for (int ix = 0; ix < image_size[0]; ++ix) {
        tomo_out.aaafI[iz][iy][ix] = aaaiBasinId[iz][iy][ix];
        // Special case: Sometimes we want undefined voxels (ie voxels whose
        // saliency exceeds settings.connect_threshold_saliency)
        // to have a brightness which is equal to the maximum brightness
        // in the image (which is num_clusters), plus one.
        // This makes it easier to visualize the clusters with low ID numbers.
        // Only do this if settings.undefined_voxels_are_max == true.
        if (aaaiBasinId[iz][iy][ix] == -1) {
          if (settings.undefined_voxels_are_max)
            tomo_out.aaafI[iz][iy][ix] = max_label + 1;
          else
            tomo_out.aaafI[iz][iy][ix] = settings.undefined_voxel_brightness;
        }
      }
    }
  }

  Dealloc3D(aaaiBasinId);

  // Did the user supply a mask?
  // Watershed() intentionally does not modify voxels which lie 
  // outside the mask.  These voxels will have random undefined values 
  // unless we assign them manually.  We do this below.
  for (int iz = 0; iz < image_size[2]; iz++)
    for (int iy = 0; iy < image_size[1]; iy++)
      for (int ix = 0; ix < image_size[0]; ix++)
        if (mask.aaafI && (mask.aaafI[iz][iy][ix] == 0.0))
          tomo_out.aaafI[iz][iy][ix] = settings.undefined_voxel_brightness;
} //HandleWatershed()






void
HandleClusterConnected(const Settings &settings,
                       MrcSimple &tomo_in,
                       MrcSimple &tomo_out,
                       MrcSimple &mask,
                       float voxel_width[3])
{
  const vector<vector<array<float, 3> > > *pMustLinkConstraints = nullptr;

  if (settings.must_link_constraints.size() > 0)
    pMustLinkConstraints = &settings.must_link_constraints;

  int image_size[3];
  for (int d = 0; d < 3; d++)
    image_size[d] = tomo_in.header.nvoxels[d];

  vector<array<float, 3> > cluster_centers;
  vector<float> cluster_sizes;
  vector<float> cluster_saliencies;

  // Create a temporary array to store the cluster membership for each voxel.
  // Because the number or clusters could (conceivably) exceed 10^6, we
  // should not make this a table of ints or floats.  Instead use "ptrdiff_t".

  ptrdiff_t ***aaaiClusterId = Alloc3D<ptrdiff_t>(tomo_in.header.nvoxels);

  // (Later on we will copy the contents of aaaiClusterId into tomo_out.aaafI
  //  which will be written to a file later.  This way the end-user can
  //  view the results.)

  ptrdiff_t undefined_voxel_brightness = settings.undefined_voxel_brightness;
  if (settings.undefined_voxels_are_max)
    undefined_voxel_brightness = -1;

  size_t num_clusters =
    ClusterConnected(tomo_in.header.nvoxels, //image size
                     tomo_in.aaafI, //<-saliency
                     aaaiClusterId, //<-which cluster does each voxel belong to?  (results will be stored here)
                     mask.aaafI,
                     settings.connect_threshold_saliency,
                     static_cast<array<float, 3> ***>(nullptr),
                     -std::numeric_limits<float>::infinity(),
                     -std::numeric_limits<float>::infinity(),
                     false, //normal default value for this (ignored) parameter
                     static_cast<float****>(nullptr),
                     -std::numeric_limits<float>::infinity(),
                     -std::numeric_limits<float>::infinity(),
                     true,  //normal default value for this (ignored) parameter
                     1,
                     static_cast<ptrdiff_t>(-1), // an impossible value
                     &cluster_centers,
                     &cluster_sizes,
                     &cluster_saliencies,
                     ClusterSortCriteria::SORT_BY_SIZE,
                     static_cast<float***>(nullptr),
                     #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
                     static_cast<array<float, 3> ***>(nullptr),
                     #endif
                     pMustLinkConstraints,
                     true, //(clusters begin at regions of high saliency)
                     &cerr);  //!< print progress to the user

  ptrdiff_t max_label = aaaiClusterId[0][0][0];
  for (int iz = 0; iz < image_size[2]; ++iz)
    for (int iy = 0; iy < image_size[1]; ++iy)
      for (int ix = 0; ix < image_size[0]; ++ix)
        max_label = std::max(max_label, aaaiClusterId[iz][iy][ix]);

  // Now, copy the contents of aaaClusterId into tomo_out.aaafI
  for (int iz = 0; iz < image_size[2]; ++iz) {
    for (int iy = 0; iy < image_size[1]; ++iy) {
      for (int ix = 0; ix < image_size[0]; ++ix) {
        tomo_out.aaafI[iz][iy][ix] = aaaiClusterId[iz][iy][ix];
        // Special case: Sometimes we want undefined voxels (ie voxels whose
        // saliency exceeds settings.connect_threshold_saliency)
        // to have a brightness which is equal to the maximum brightness
        // in the image (which is num_clusters), plus one.
        // This makes it easier to visualize the clusters with low ID numbers.
        // Only do this if settings.undefined_voxels_are_max == true.
        if (aaaiClusterId[iz][iy][ix] == -1) {
          if (settings.undefined_voxels_are_max)
            tomo_out.aaafI[iz][iy][ix] = max_label + 1;
          else
            tomo_out.aaafI[iz][iy][ix] = settings.undefined_voxel_brightness;
        }
      }
    }
  }

  Dealloc3D(aaaiClusterId);

} //HandleClusterConnected()





void
HandleTV(const Settings &settings,
         MrcSimple &tomo_in,
         MrcSimple &tomo_out,
         MrcSimple &mask,
         float voxel_width[3])
{
  cerr << "filter_type = surface ridge detector\n";

  const vector<vector<array<float, 3> > > *pMustLinkConstraints = nullptr;
  
  if (settings.must_link_constraints.size() > 0)
    pMustLinkConstraints = &settings.must_link_constraints;

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

  // Now use Alloc3D() to allocate space for both aafGradient and aaaafGradient.

  aaaafGradient = Alloc3D<array<float, 3> >(image_size);

  // The storage requirement for Hessians (6 floats) is large enough that
  // I decided to represent hessians using a CompactMultiChannelImage3D.
  // Internally this is a 4-dimensional array, however the last dimension
  // is only allocated (non-null) for voxels which were selected by the user
  // (ie voxels for which the mask is non-zero).  This can reduce memory usage
  // by a factor of up to 3 (assuming floats) for this array.
  CompactMultiChannelImage3D<float> hessian_tensor(6);
  hessian_tensor.Resize(image_size, mask.aaafI, &cerr);

  float filter_truncate_ratio = settings.filter_truncate_ratio;

  // How did the user specify how wide to make the filter window?
  if (filter_truncate_ratio <= 0) {
    assert(settings.filter_truncate_threshold > 0.0);
    //    filter_truncate_threshold = exp(-(1/2)*filter_truncate_ratio^2);
    //    -> filter_truncate_ratio^2 = -2*log(filter_truncate_threshold)
    filter_truncate_ratio = sqrt(-2*log(settings.filter_truncate_threshold));
  }

  MrcSimple tomo_background;
  bool subtract_background = (settings.width_b[0] > 0.0);
  if (subtract_background) {
    tomo_background = tomo_in;

    int truncate_halfwidth = floor(settings.width_b[0] *
                                   filter_truncate_ratio);

    ApplyGauss(image_size,
               tomo_in.aaafI,
               tomo_background.aaafI,
               mask.aaafI,
               settings.width_b[0],
               truncate_halfwidth,
               settings.normalize_near_boundaries,
               &cerr);

    truncate_halfwidth = floor(settings.width_a[0] *
                               filter_truncate_ratio);

    ApplyGauss(image_size,
               tomo_in.aaafI,
               tomo_out.aaafI,
               mask.aaafI,
               settings.width_a[0],
               truncate_halfwidth,
               settings.normalize_near_boundaries,
               &cerr);
  }

  // Using this blurred image, calculate the 2nd derivative matrix everywhere:

  cerr << "Applying a Gaussian blur of width sigma="
       << sigma*voxel_width[0] << " (physical distance)" << endl;

  if ((image_size[0] < 3) ||
      (image_size[1] < 3) ||
      (image_size[2] < 3))
    throw VisfdErr("Error: Ridge-detection requires an image that is at least 3 voxels\n"
                   "       wide in the x,y,z directions.\n");

  CalcHessian(image_size,
              tomo_in.aaafI,
              aaaafGradient,
              hessian_tensor.aaaafI,
              mask.aaafI,
              sigma,
              filter_truncate_ratio,
              &cerr);

  cerr << "Diagonalizing the Hessians" << endl;

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

        if (! hessian_tensor.aaaafI[iz][iy][ix]) //ignore voxels that are either
          continue;                         //in the mask or on the boundary

        float eivals[3];
        float eivects[3][3];

        ConvertFlatSym2Evects3(hessian_tensor.aaaafI[iz][iy][ix],
                               eivals,
                               eivects,
                               eival_order);

        float score;

        // DEBUG: REMOVE THE NEXT IF STATMENT AFTER DEBUGGING IS FINISHED
        #ifndef NDEBUG
        if ((ix == image_size[0] / 2) &&
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



        if (settings.filter_type == Settings::CURVE) {
          score = ScoreHessianLinear(eivals,
                                     aaaafGradient[iz][iy][ix]);
        }
        else if (settings.filter_type == Settings::SURFACE_EDGE) {
          score = sqrt(SQR(aaaafGradient[iz][iy][ix][0]) +
                       SQR(aaaafGradient[iz][iy][ix][1]) +
                       SQR(aaaafGradient[iz][iy][ix][2]));
        }
        else if (settings.filter_type == Settings::SURFACE_RIDGE) {
          score = ScoreHessianPlanar(eivals,
                                     aaaafGradient[iz][iy][ix]);
        }
        else
          assert(0);  //the filter_type should be one of the options above


        float peak_height = 1.0;
        if (tomo_background.aaafI)
          peak_height = (tomo_in.aaafI[iz][iy][ix] -
                         tomo_background.aaafI[iz][iy][ix]);
        score *= peak_height;

        tomo_out.aaafI[iz][iy][ix] = score;


        #if 0
        //COMMENTING OUT
        //I might use this code eventually, but not before testing it first.
        if (settings.discard_voxels_not_on_ridge) {
          float gradient_along_v1 = DotProduct3(grad, eivects[0]);
          float distance_to_ridge;
          if (lambda1 != 0)
            distance_to_ridge = abs(gradient_along_v1 / lambda1);
          else
            distance_to_ridge = std::numeric_limits<float>::infinity();
          bool ridge_nearby = true;
          if ((settings.max_distance_to_feature > 0.0) &&
              (abs(distance_to_ridge) > settings.max_distance_to_feature))
            ridge_nearby = false;
          if (! ridge_nearby)
            out_tomo.aaafI[iz][iy][ix] = 0.0;  //then ignore discard this voxel
        } //if (settings.discard_voxels_not_on_ridge)
        #endif



        if (settings.filter_type != Settings::SURFACE_EDGE) {
          // Currently the aaaafDirection = aaaafGradient.
          // If the directional feature that we want to detect is
          // the edge of some sharp boundary between a light and dark region,
          // then leave aaaafDirection alone, because it is already pointing
          // in the direction of increasing brightness (along the gradient).
          //
          // However, if we want to detect ridges (thin membranes or curves)
          // then we should have aaaafDirection point in the direction
          // of the principal axis of the hessian.
          aaaafDirection[iz][iy][ix][0] = eivects[0][0];
          aaaafDirection[iz][iy][ix][1] = eivects[0][1];
          aaaafDirection[iz][iy][ix][2] = eivects[0][2];
        }


      } //for(int ix=0; ix < image_size[0]; ix++)
    } //for(int iy=0; iy < image_size[1]; iy++)
  } //for(int iz=0; iz < image_size[2]; iz++)




  { // Use thresholding to reduce the number of voxels that we have to consider
    
    float hessian_score_threshold = settings.hessian_score_threshold;
    if (settings.hessian_score_threshold_is_a_fraction) {
      float hessian_score_fraction = settings.hessian_score_threshold;
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




  // To save memory we will re-use the large array we used to store
  // the hessian that we calculated earlier.  But we will rename
  // this array "aaaafVoteTensor" to (hopefully) make it more clear in
  // the code that this tensor no longer necessarily stores the hessian.
  float ****aaaafVoteTensor = hessian_tensor.aaaafI;



  if (settings.tv_sigma > 0.0) {

    if (settings.load_intermediate_fname_base == "") {

      // Then the user is NOT loading (previously generated) files containing
      // the features we are looking for (like surfaces or curves).
      //
      // Hence we must calculate these features (using tensor-voting).

      assert(filter_truncate_ratio > 0);

      TV3D<float, int, array<float,3>, float* >
        tv(settings.tv_sigma,
           settings.tv_exponent,
           settings.tv_truncate_ratio);

      tv.TVDenseStick(image_size,
                      tomo_out.aaafI,
                      aaaafDirection,
                      aaaafVoteTensor,
                      mask.aaafI,
                      mask.aaafI,
                      (settings.filter_type == Settings::CURVE),
                      //settings.hessian_score_threshold,
                      false, // don't boost votes from voxels near mask boundary
                      false, // don't diagonalize each tensor afterwards
                      &cerr);

    } // if (settings.load_intermediate_fname_base == "")

    else {

      // Load the contents of the aaaafVoteTensor array
      // from a series of files (generated previously).
      //
      // First figure out the file name(s) we are going to use.
      // If the "settings.load_intermediate_fname_base" file name ends in
      // ".rec" or ".mrc", then remove these these suffixes.
      // First figure out the file name(s) we are going to use.
      string fname_base = settings.load_intermediate_fname_base;
      // Then load the files corresponding to the tensor
      for (int d = 0; d < 6; d++) {
        stringstream fname_ss;
        fname_ss << fname_base << "_tensor_" << d << ".rec";
        cerr << "loading \"" << fname_ss.str() << "\"" << endl;
        MrcSimple tomo_tmp;
        tomo_tmp.Read(fname_ss.str(), false);
        for(int iz=0; iz < image_size[2]; iz++)
          for(int iy=0; iy < image_size[1]; iy++)
            for(int ix=0; ix < image_size[0]; ix++)
              if (aaaafVoteTensor[iz][iy][ix])
                aaaafVoteTensor[iz][iy][ix][d] = tomo_tmp.aaafI[iz][iy][ix];
      }
    } // else clause for if (settings.load_intermediate_fname_base == "")




    // Now that the feature tensor has been calculated,
    // calculate the saliency of each voxel using this tensor.
    for(int iz=0; iz < image_size[2]; iz++) {
      for(int iy=0; iy < image_size[1]; iy++) {
        for(int ix=0; ix < image_size[0]; ix++) {
          if (aaaafVoteTensor[iz][iy][ix]) {
            float diagonalized_hessian[6];
            DiagonalizeFlatSym3(aaaafVoteTensor[iz][iy][ix],
                                diagonalized_hessian,
                                eival_order);
            float score;
            if (settings.filter_type == Settings::CURVE)
              score = ScoreTensorLinear(diagonalized_hessian);
            else
              score = ScoreTensorPlanar(diagonalized_hessian);
            float peak_height = 1.0;
            if (tomo_background.aaafI)
              peak_height = (tomo_in.aaafI[iz][iy][ix] -
                             tomo_background.aaafI[iz][iy][ix]);
            score *= peak_height;
            tomo_out.aaafI[iz][iy][ix] = score;
          }
        }
      }
    }
  } // if (settings.tv_sigma > 0.0)




  // Save the contents of the vote_tensor.aaaafI array for future use?
  if (settings.save_intermediate_fname_base != "") {

    // Then we will save the 4-dimensional vote_tensor.aaaafI array
    // in a series of 3-dimensional image files, with names ending in
    // "_tensor_0.rec", "_tensor_1.rec", ..., "_tensor_5.rec"
    //
    // First figure out the file name(s) we are going to use.
    // If the "settings.save_intermediate_fname_base" file name ends in
    // ".rec" or ".mrc", then remove these these suffixes.
    string fname_base = settings.save_intermediate_fname_base;
    // Then save the files corresponding to the tensor
    for (int d = 0; d < 6; d++) {
      stringstream fname_ss;
      fname_ss << fname_base << "_tensor_" << d << ".rec";
      cerr << "allocating space for \"" << fname_ss.str() << "\"" << endl;
      MrcSimple tomo_tmp = tomo_out;
      for(int iz=0; iz < image_size[2]; iz++)
        for(int iy=0; iy < image_size[1]; iy++)
          for(int ix=0; ix < image_size[0]; ix++)
            if (aaaafVoteTensor[iz][iy][ix])
              tomo_tmp.aaafI[iz][iy][ix] = aaaafVoteTensor[iz][iy][ix][d];
      cerr << "writing \"" << fname_ss.str() << "\"" << endl;
      tomo_tmp.Write(fname_ss.str());
    }

  } //if (settings.save_intermediate_fname_base != "")


  float ***aaafSaliency = nullptr;

  if (settings.cluster_connected_voxels)
  {
    // Later on, I will overwrite the contents of tomo_out.
    // So I need to copy it somewhere before this happens.
    // Unfortunately, memory is scarce due to all the arrays we have allocated.
    tomo_in = tomo_out; // (horrible hack.  I should not modify tomo_in.  I
                        // should allocate a new 3D array to store the saliency)
    aaafSaliency = tomo_in.aaafI;
    // Copy the principal eigenvector of vote_tensor into aaaafDirection
    for(int iz=0; iz < image_size[2]; iz++) {
      for(int iy=0; iy < image_size[1]; iy++) {
        for(int ix=0; ix < image_size[0]; ix++) {
          if (aaaafVoteTensor[iz][iy][ix]) {
            float eivals[3];
            float eivects[3][3];
            ConvertFlatSym2Evects3(aaaafVoteTensor[iz][iy][ix],
                                   eivals,
                                   eivects,
                                   eival_order);
            for (int d=0; d<3; d++)
              aaaafDirection[iz][iy][ix][d] = eivects[0][d];
          }
        }
      }
    }

    vector<array<float, 3> > cluster_centers;
    vector<float> cluster_sizes;
    vector<float> cluster_saliencies;

    // Create a temporary array to store the cluster membership for each voxel.
    // Because the number or clusters could (conceivably) exceed 10^6, we
    // should not make this a table of ints or floats.  Instead use "ptrdiff_t".

    ptrdiff_t ***aaaiClusterId = nullptr;
    aaaiClusterId = Alloc3D<ptrdiff_t>(image_size);

    // (Later on we will copy the contents of aaaiClusterId into tomo_out.aaafI
    //  which will be written to a file later.  This way the end-user can
    //  view the results.)

    size_t num_clusters =

      ClusterConnected(image_size, //image size
                       aaafSaliency,
                       aaaiClusterId, //<-which cluster does each voxel belong to?  (results will be stored here)
                       mask.aaafI,
                       settings.connect_threshold_saliency,
                       aaaafDirection,
                       settings.connect_threshold_vector_saliency,
                       settings.connect_threshold_vector_neighbor,
                       false, //eigenvector signs are arbitrary so ignore them
                       aaaafVoteTensor,
                       settings.connect_threshold_tensor_saliency,
                       settings.connect_threshold_tensor_neighbor,
                       true,  //the tensor should be positive definite near the target
                       1,
                       static_cast<ptrdiff_t>(-1), // an impossible value
                       &cluster_centers,
                       &cluster_sizes,
                       &cluster_saliencies,
                       ClusterSortCriteria::SORT_BY_SIZE,
                       static_cast<float***>(nullptr),
                       #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
                       aaaafDirection,
                       #endif
                       pMustLinkConstraints,
                       true, //(clusters begin at regions of high saliency)
                       &cerr);  //!< print progress to the user

    ptrdiff_t max_label = aaaiClusterId[0][0][0];
    for (int iz = 0; iz < image_size[2]; ++iz)
      for (int iy = 0; iy < image_size[1]; ++iy)
        for (int ix = 0; ix < image_size[0]; ++ix)
          max_label = std::max(max_label, aaaiClusterId[iz][iy][ix]);

    // Now, copy the contents of aaaClusterId into tomo_out.aaafI
    for (int iz = 0; iz < image_size[2]; ++iz) {
      for (int iy = 0; iy < image_size[1]; ++iy) {
        for (int ix = 0; ix < image_size[0]; ++ix) {
          tomo_out.aaafI[iz][iy][ix] = aaaiClusterId[iz][iy][ix];
          // Special case: Sometimes we want undefined voxels (ie voxels whose
          // saliency exceeds settings.connect_threshold_saliency)
          // to have a brightness which is equal to the maximum brightness
          // in the image (which is num_clusters), plus one.
          // This makes it easier to visualize the clusters with low ID numbers.
          // Only do this if settings.undefined_voxels_are_max == true.
          if (aaaiClusterId[iz][iy][ix] == -1) {
            if (settings.undefined_voxels_are_max)
              tomo_out.aaafI[iz][iy][ix] = max_label + 1;
            else
              tomo_out.aaafI[iz][iy][ix] = settings.undefined_voxel_brightness;
          }
        }
      }
    }

    Dealloc3D(aaaiClusterId);

  } // if (settings.cluster_connected_voxels)



  //Did the user ask us to generate output files containing surface orientation?
  if (settings.out_normals_fname != "")
  {
    vector<array<float,3> > coords;
    vector<array<float,3> > norms;
    // where is the lookup table to indicate which cluster a voxel belongs to?
    float ***aaafVoxel2Cluster = nullptr;
    if (settings.cluster_connected_voxels)
      aaafVoxel2Cluster = tomo_out.aaafI; //tomo_out set by ClusterConnected()

    for (int iz = 0; iz < image_size[2]; ++iz) {
      for (int iy = 0; iy < image_size[1]; ++iy) {
        for (int ix = 0; ix < image_size[0]; ++ix) {
          if (mask.aaafI && (mask.aaafI[iz][iy][ix] == 0.0))
            continue;
          if (! settings.cluster_connected_voxels) {
            // Then include all of the voxels in the image
            array<float,3> xyz;
            xyz[0] = ix * voxel_width[0];
            xyz[1] = iy * voxel_width[1];
            xyz[2] = iz * voxel_width[2];
            coords.push_back(xyz);
            array<float,3> normal;
            for (int d = 0; d < 3; d++)
              normal[d] = aaaafDirection[iz][iy][ix][d];
            norms.push_back(normal);
            continue;
          }
          else if (settings.select_cluster != aaafVoxel2Cluster[iz][iy][ix])
            continue;


          // If we are here then the user wanted us to cluster the voxels,
          // only select voxels from one of these clusters,
          // and this voxel is one of the voxels from that cluster.
          // In that case the cluster is supposed to consist of voxels which
          // lie on one of the surfaces detected in this image.

          // The surface is supposed to be flat and 2-dimensional.
          // But some of these selected voxels lie quite far from the surface.
          //
          // (Example:
          //  For lipid bilayers in cryo-EM tomogramd, after tensor-voting
          //  is finished, the resulting surfaces can be up to 15 voxels thick.
          //  It's because all of those voxels have saliencies above the cutoff,
          //  even though the actual membrane in the image is much thinner.)
          //
          // We need to figure out where the memrane surface actually is,
          // discard voxels that are far away from it, and place the positions
          // of the remaining voxels somewhere on the surface itself.
          // The points on the pointcloud should lie on this surface.

          array<float, 3> xyz;
          xyz[0] = ix;
          xyz[1] = iy;
          xyz[2] = iz;
          array<float,3> normal;
          float norm=length3<float,array<float,3> >(aaaafDirection[iz][iy][ix]);
          for (int d = 0; d < 3; d++) {
            normal[d] = aaaafDirection[iz][iy][ix][d] / norm;
            normal[d] *= aaafSaliency[iz][iy][ix];
          }

          if (settings.surface_normal_curve_ds > 0.0) {
            vector<float> vS;
            vector<float> vWeights;
            vector<array<float,3> > vxyz;
            float ds = settings.surface_normal_curve_ds;//=integration step-size
            // s = length of the curve so far
            float s;
            // r = 3D position on the curve
            array<float, 3> r = {float(ix), float(iy), float(iz)};
            // drds = tangent (direction) of the curve
            array<float, 3> drds;
            array<int, 3> ixyz = {ix, iy, iz};
            s = 0.0;
            while ((0 <= ixyz[0]) && (ixyz[0] < image_size[0]) &&
                   (0 <= ixyz[1]) && (ixyz[1] < image_size[1]) &&
                   (0 <= ixyz[2]) && (ixyz[2] < image_size[2]) &&
                   ((! mask.aaafI) ||
                    (mask.aaafI[ixyz[2]][ixyz[1]][ixyz[0]] != 0.0)) &&
                   //(std::abs(s) <= settings.max_distance_to_feature) &&
                   //(aaafSaliency[ixyz[2]][ixyz[1]][ixyz[0]] >=
                   // settings.connect_threshold_saliency))
                   (aaafVoxel2Cluster[ixyz[2]][ixyz[1]][ixyz[0]] ==
                    aaafVoxel2Cluster[iz][iy][ix]))
            {
              vS.push_back(s);
              vxyz.push_back(r);
              vWeights.push_back(aaafSaliency[ixyz[2]][ixyz[1]][ixyz[0]]);
              norm=length3<float,array<float,3> >(aaaafDirection[ixyz[2]][ixyz[1]][ixyz[0]]);
              for (int d = 0; d < 3; d++)
                drds[d] = aaaafDirection[ixyz[2]][ixyz[1]][ixyz[0]][d] / norm;
              s += ds;
              for (int d = 0; d < 3; d++) {
                r[d] += ds * drds[d];
                ixyz[d] = int(round(r[d]));
              }
            }
            // Now do the same in the opposite (-s) direction.
            vector<float> _vS;
            vector<float> _vWeights;
            vector<array<float,3> > _vxyz;
            r = {float(ix), float(iy), float(iz)};
            ixyz = {ix, iy, iz};
            s = 0.0;
            while (true)
            {
              norm=length3<float,array<float,3> >(aaaafDirection[ixyz[2]][ixyz[1]][ixyz[0]]);
              for (int d = 0; d < 3; d++)
                drds[d] = aaaafDirection[ixyz[2]][ixyz[1]][ixyz[0]][d] / norm;
              s -= ds;
              for (int d = 0; d < 3; d++) {
                r[d] -= ds * drds[d];
                ixyz[d] = int(round(r[d]));
              }
              if ((ixyz[0] < 0) || (image_size[0] <= ixyz[0]) ||
                  (ixyz[1] < 0) || (image_size[1] <= ixyz[0]) ||
                  (ixyz[2] < 0) || (image_size[2] <= ixyz[0]))
                break;
              if (mask.aaafI &&
                  (mask.aaafI[ixyz[2]][ixyz[1]][ixyz[0]] == 0.0))
                break;
              //if ((std::abs(s) > settings.max_distance_to_feature) ||
              //    (aaafSaliency[ixyz[2]][ixyz[1]][ixyz[0]] <
              //     settings.connect_threshold_saliency))
              //  break;
              if (aaafVoxel2Cluster[ixyz[2]][ixyz[1]][ixyz[0]] !=
                  aaafVoxel2Cluster[iz][iy][ix])
                break;
              _vS.push_back(s);
              _vxyz.push_back(r);
              _vWeights.push_back(aaafSaliency[ixyz[2]][ixyz[1]][ixyz[0]]);
            }
            // Now insert the contents from _vS at the beginning of vS.
            // making sure the the "s" values will be in increasing order.
            vS.insert(vS.begin(), _vS.rbegin(), _vS.rend());
            // Do the same for the other arrays.
            vxyz.insert(vxyz.begin(), _vxyz.rbegin(), _vxyz.rend());
            vWeights.insert(vWeights.begin(), _vWeights.rbegin(), _vWeights.rend());
            // Find the (weighteD) average S value and determine where it is.
            float sum_s = 0.0;
            float sum_w = 0.0;
            assert(vS.size() == vxyz.size());
            assert(vS.size() == vWeights.size());
            for (int i = 0; i < vS.size(); i++) {
              sum_s += vWeights[i] * vS[i];
              sum_w += vWeights[i];
            }
            float ave_s = sum_s / sum_w;
            int i = 0;
            while (i+1 < vS.size()) {
              i++;
              if ((vS[i-1] <= ave_s) && (ave_s <= vS[i]))
                break;
            }
            assert(i < vS.size());
            // i is now the index into the vxyz array indicating the
            // coordinates of the voxel which lies closest to the average
            // position along this curve.
            // Now store the coordinates of that voxel in xyz[]
            // so we can add it to the point cloud later.
            for (int d = 0; d < 3; d++)
              ixyz[d] = int(round(vxyz[i][d]));
            float norm = length3<float,array<float,3> >(aaaafDirection[ixyz[2]][ixyz[1]][ixyz[0]]);
            for (int d = 0; d < 3; d++)
              normal[d] = aaaafDirection[ixyz[2]][ixyz[1]][ixyz[0]][d] / norm;
            for (int d = 0; d < 3; d++) {
              if (i+1 < vS.size()) {
                xyz[d] = (vxyz[i][d] +
                          (vxyz[i+1][d]-vxyz[i][d])*((ave_s-vS[i]) /
                                                     (vS[i+1]-vS[i])));
              }
              else
                xyz[d] = vxyz[i][d];
              // We also want to use the normal direction at that location.
              // But scale the direction so that the normal magnitude is
              // proportional to the saliency at the --original-- voxel's
              // location (at ix, iy, iz).
              // (We will eventually include it with the point-cloud data also.)
              normal[d] *= aaafSaliency[iz][iy][ix];
            }
          } //if (settings.surface_normal_curve_ds > 0.0)


          // Do we want to futher refine the surface using the local Hessian?
          // (Comment: Most of the time, this makes the resulting surface
          //  -slightly- smoother.  However in some cases, at surface edges
          //  or other places where the ridge is not well defined, it
          //  causes more problems than it solves.  So I haven't decided if
          //  I should keep this feature or disable it.  -Andrew 2021-7-19)
          if (settings.surface_find_ridge) {
            int ix0 = int(round(xyz[0]));
            int iy0 = int(round(xyz[1]));
            int iz0 = int(round(xyz[2]));
            // For sub-voxel accuracy, try to find the closest ridge
            // to the current voxel location.
            // I assume it is located at the ridges of the saliency.
            // So figure out where those ridges are, discard voxels that
            // are too far away from them, and project the positions of
            // the remaining voxels onto the nearest point on the ridge surface.
            //
            // Calculate the gradient and hessian of the saliency
            float hessian[3][3];
            float gradient[3];
            assert(aaafSaliency);
            CalcHessianFiniteDifferences(aaafSaliency,
                                         ix0, iy0, iz0,
                                         hessian,
                                         image_size);
            CalcGradientFiniteDifferences(aaafSaliency,
                                          ix0, iy0, iz0,
                                          gradient,
                                          image_size);
            float eivects[3][3]; // eivect[0]=direction of largest 2nd derivative
            float eivals[3];     // eivals[0]=the 2nd derivative in that direction
            DiagonalizeSym3(hessian,
                            eivals,
                            eivects,
                            selfadjoint_eigen3::DECREASING_ABS_EIVALS);
            float gradient_along_v1 = DotProduct3(gradient, eivects[0]);
            // Now deal with sign ambiguity of the eigenvectors.
            // Choose the sign of the eigenvector that points along the gradient.
            if (gradient_along_v1 < 0.0) {
              // We want the dot product to be > 0
              gradient_along_v1 = -gradient_along_v1;
              for (int d = 0; d < 3; d++)
                eivects[0][d] = -eivects[0][d];
            }
            else if (gradient_along_v1 == 0.0)
              continue; // then ignore this voxel
            float distance_to_ridge;
            if (eivals[0] != 0)
              distance_to_ridge = gradient_along_v1 / eivals[0];
            else
              distance_to_ridge = std::numeric_limits<float>::infinity();
            // Only select voxels that lie near the saliency ridge.
            // If we are too far away from the ridge, discard this voxel
            if ((settings.max_distance_to_feature > 0.0) &&
                (abs(distance_to_ridge) > settings.max_distance_to_feature))
              continue;
            xyz[0] = ix0 - distance_to_ridge * eivects[0][0];
            xyz[1] = iy0 - distance_to_ridge * eivects[0][1];
            xyz[2] = iz0 - distance_to_ridge * eivects[0][2];
            // make sure the point lies within the image boundaries
            if ((xyz[0] < 0.0) || (image_size[0] < xyz[0]) ||
                (xyz[1] < 0.0) || (image_size[1] < xyz[1]) ||
                (xyz[2] < 0.0) || (image_size[2] < xyz[2]))
              continue; // if not, ignore this point
            for (int d = 0; d < 3; d++)
              // We want to export coordinates in units of physical distance
              xyz[d] *= voxel_width[d];
            // Note that eivects[0][] and aaaafDirection[] are both
            // vectors that store the surface normal.
            // I am guessing that aaaafDirection[] is presumably more accurate
            // choice because it comes directly from tensor-voting (although
            // the magnitude is hard to interpret, so you can't use it for
            // finding the ridge position).
            // So we set the surface normal to aaaafDirection[],not eivects[0][]
          } //if (settings.surface_find_ridge)

          coords.push_back(xyz); // add it to the point cloud
          norms.push_back(normal);
        } //for (int ix = 0; ix < image_size[0]; ++ix)
      } //for (int iy = 0; iy < image_size[1]; ++iy)
    } //for (int iz = 0; iz < image_size[2]; ++iz)

    WriteOrientedPointCloudPLY(settings.out_normals_fname,
                               coords,
                               norms);

  } //if (settings.out_normals_fname != "")


  ////Did the user ask us to generate any output files?
  //if ((settings.out_normals_fname == "") &&
  //    (settings.out_file_name != ""))
  //  settings.out_normals_fname =
  //    settings.out_file_name + string(".bnpts");


  Dealloc3D(aaaafGradient);

} //HandleTV()



void
HandleBinning(const Settings &settings,
              MrcSimple &tomo_in,
              MrcSimple &mask,
              float voxel_width[3])
{
  int nvoxels_resized[3];
  for (int d=0; d < 3; d++)
    nvoxels_resized[d] = (tomo_in.header.nvoxels[d] /
                          settings.resize_with_binning);
  double _voxel_width;
  if (settings.voxel_width > 0)
    _voxel_width = settings.voxel_width;
  else
    _voxel_width = tomo_in.header.cellA[0]/tomo_in.header.nvoxels[0];

  MrcSimple tomo_tmp = tomo_in;
  tomo_tmp.Resize(nvoxels_resized);
  // Now bin the original image, and save the result in tomo_tmp.aaafI
  BinArray3D(tomo_in.header.nvoxels,
             nvoxels_resized,
             tomo_in.aaafI,   // <-- source
             tomo_tmp.aaafI); // <-- dest

  _voxel_width *= settings.resize_with_binning;
  voxel_width[0] = _voxel_width;
  voxel_width[1] = _voxel_width;
  voxel_width[2] = _voxel_width;

  // Now update the voxel size related information
  // in the header section of "tomo_tmp".
  tomo_tmp.header.cellA[0] = _voxel_width * nvoxels_resized[0];
  tomo_tmp.header.cellA[1] = _voxel_width * nvoxels_resized[1];
  tomo_tmp.header.cellA[2] = _voxel_width * nvoxels_resized[2];

  // Now copy the contents of tomo_tmp into tomo_in.
  tomo_in.swap(tomo_tmp);
  // (Note: The old contents of "tomo_in.aaafI" will be freed
  //  when tomo_tmp is destroyed or modified with Resize().)

  // Did the user supply a mask?
  // If so, we should resize (bin) the mask as well.
  if (mask.header.nvoxels[0] > 0) {
    // Since we swapped tomo_in <--> tomo_tmp,
    // the tomo_tmp.aaafI array is no longer the correct size.
    // We can fix that by invoking Resize() again.
    tomo_tmp.Resize(nvoxels_resized);
    // Now bin the mask, and save the result in tomo_tmp.aaafI
    BinArray3D(mask.header.nvoxels,
               nvoxels_resized,
               mask.aaafI,      // <-- source
               tomo_tmp.aaafI); // <-- dest
    // Now update the voxel size related information
    // in the header section of "tomo_tmp".
    tomo_tmp.header.cellA[0] = _voxel_width * nvoxels_resized[0];
    tomo_tmp.header.cellA[1] = _voxel_width * nvoxels_resized[1];
    tomo_tmp.header.cellA[2] = _voxel_width * nvoxels_resized[2];

    // Copy the contents of tomo_tmp into mask.
    mask.swap(tomo_tmp);

    //(The old contents of "mask" will be freed when tomo_tmp is destroyed.)
  }
} // HandleBinning()



void
DetermineVoxelWidth(const Settings &settings,
                    MrcSimple &tomo_in,
                    MrcSimple &mask,
                    float voxel_width[3])
{
  voxel_width[0] = 1.0;
  voxel_width[1] = 1.0;
  voxel_width[2] = 1.0;

  int image_size[3];
  image_size[0] = tomo_in.header.nvoxels[0];
  image_size[1] = tomo_in.header.nvoxels[1];
  image_size[2] = tomo_in.header.nvoxels[2];

  if (settings.voxel_width > 0.0) {
    // Did the user manually specify the width of each voxel?
    voxel_width[0] = settings.voxel_width;
    voxel_width[1] = settings.voxel_width;
    voxel_width[2] = settings.voxel_width;
    // If the user altered the size of the voxel width by grouping
    // multiple voxels into single voxels (through binning),
    // then the effective size of the new voxels is larger than
    // the original voxels.  It is assumed that the voxel_width specified
    // by the user is the voxel width of the original image, not the binned
    // image.  So we must increase this voxel width to compensate.
    // (We don't have to do this if the voxel width is stored
    //  header of the MRC file, because that width is inferred from the
    //  header.cellA[] and image_size[] numbers, which are currently correct.)
    if (settings.resize_with_binning > 0) {
      for (int d=0; d < 3; d++)
        voxel_width[d] *= settings.resize_with_binning;
    }
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
    throw VisfdErr(err_msg.str());
  }
} //DetermineVoxelWidth()


