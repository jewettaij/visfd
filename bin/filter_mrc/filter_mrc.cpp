// (Note: For gcc version 4.8.3, you must compile using: g++ -std=c++11)


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

#include <err_report.hpp>
#include <alloc2d.hpp>
#include <alloc3d.hpp>
#include <filter1d.hpp>
#include <filter2d.hpp>
#include <filter3d.hpp>
#include <multichannel_image3d.hpp>
#include <lin3_utils.hpp>
#include <threshold.hpp>
#include <mrc_simple.hpp>
#include <random_gen.h>
#include "settings.hpp"
#include "filter3d_variants.hpp"
#include "unsupported.hpp"


string g_program_name("filter_mrc.cpp");
string g_version_string("0.12.6");
string g_date_string("2018-1-21");



static void
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




static void
HandleGauss(Settings settings,
            MrcSimple &tomo_in,
            MrcSimple &tomo_out,
            MrcSimple &mask,
            float voxel_width[3])
{
  cerr << "filter_type = Gaussian\n";

  float A; // let the user know what A coefficient was used

  A = ApplyGauss3D(tomo_in.header.nvoxels,
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




static void
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




static void
HandleDog(Settings settings,
          MrcSimple &tomo_in,
          MrcSimple &tomo_out,
          MrcSimple &mask,
          float voxel_width[3])
{
  cerr << "filter_type = Difference of Gaussians (DOG)\n";

  float A, B;

  ApplyDog3D(tomo_in.header.nvoxels,
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





static void
HandleDogScaleFree(Settings settings,
                   MrcSimple &tomo_in,
                   MrcSimple &tomo_out,
                   MrcSimple &mask,
                   float voxel_width[3])
{
  cerr << "filter_type = Fast Laplacian of Gaussians (LoG)\n"
       << "  (This will be approximated as a Difference of Gaussians,\n"
       << "   as explained below.)\n";
  //cerr << "filter_type = Difference of Gaussians Scale Free (DOGSF)\n";

  float A, B;

  ApplyLog3D(tomo_in.header.nvoxels,
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

} // HandleDogUseDelta()



/// @brief Find scale-invariant blobs in the image as a function of diameter.
///        This variant detects blobs and performs non-max suppression.
///
/// Algorithm described in:
///    Lindeberg,T., Int. J. Comput. Vision., 30(2):77-116, (1998)
/// This version refers to blobs by their diameter (instead of "sigma").
/// This version can discard blobs which overlap with existing blobs.
/// (This is sometimes called "non-max suppression".)
/// This function is not intended for public use.

template<class Scalar>
static void
BlobDogNM(int const image_size[3], //!<source image size
          Scalar const *const *const *aaafSource,   //!<source image
          Scalar const *const *const *aaafMask,     //!<ignore voxels where mask==0
          const vector<Scalar>& blob_diameters, //!< blob widths to try, ordered
          vector<array<Scalar,3> >& minima_crds, //!<store minima x,y,z coords here
          vector<array<Scalar,3> >& maxima_crds, //!<store maxima x,y,z coords here
          vector<Scalar>& minima_diameters, //!< corresponding width for that minima
          vector<Scalar>& maxima_diameters, //!< corresponding width for that maxima
          vector<Scalar>& minima_scores, //!< what was the blob's score?
          vector<Scalar>& maxima_scores, //!< (score = intensity after filtering)
          Scalar delta_sigma_over_sigma,//!< param for approximating LOG with DOG
          Scalar truncate_ratio,      //!< how many sigma before truncating?
          Scalar minima_threshold=0.5,  //!< discard blobs with unremarkable scores
          Scalar maxima_threshold=0.5,  //!< discard blobs with unremarkable scores
          bool    use_threshold_ratios=true, //!< threshold=ratio*best_score?
          Scalar sep_ratio_thresh=1.0,          //!< minimum radial separation between blobs
          Scalar nonmax_max_overlap_large=1.0,  //!< maximum volume overlap with larger blob
          Scalar nonmax_max_overlap_small=1.0,  //!< maximum volume overlap with smaller blob
          // optional arguments
          ostream *pReportProgress = NULL, //!< report progress to the user?
          Scalar ****aaaafI = NULL, //!<preallocated memory for filtered images
          Scalar **aafI = NULL     //!<preallocated memory for filtered images
          )
{

  BlobDogD(image_size,
           aaafSource,
           aaafMask,
           blob_diameters,
           minima_crds,
           maxima_crds,
           minima_diameters,
           maxima_diameters,
           minima_scores,
           maxima_scores,
           delta_sigma_over_sigma,
           truncate_ratio,
           minima_threshold,
           maxima_threshold,
           use_threshold_ratios,
           pReportProgress,
           aaaafI,
           aafI);

  if (pReportProgress)
    *pReportProgress << "----------- Removing overlapping blobs -----------\n" << endl;


  if (pReportProgress)
    *pReportProgress << "--- Discarding overlapping minima blobs ---\n";

  DiscardOverlappingBlobs(minima_crds,
                          minima_diameters,
                          minima_scores,
                          PRIORITIZE_LOW_SCORES,
                          sep_ratio_thresh,
                          nonmax_max_overlap_large,
                          nonmax_max_overlap_small,
                          pReportProgress);

  if (pReportProgress)
    *pReportProgress << "done --\n"
                     << "--- Discarding overlapping maxima blobs ---\n";

  DiscardOverlappingBlobs(maxima_crds,
                          maxima_diameters,
                          maxima_scores,
                          PRIORITIZE_HIGH_SCORES,
                          sep_ratio_thresh,
                          nonmax_max_overlap_large,
                          nonmax_max_overlap_small,
                          pReportProgress);

} //BlobDogNM()





/// @brief Find scale-invariant blobs in the image as a function of diameter.
///        In this minor variant, the user can specify the filter window width
///        either in units of sigma, or by specifying the decay threshold.
///        This function was not intended for public use.
///
/// Algorithm described in:
///    Lindeberg,T., Int. J. Comput. Vision., 30(2):77-116, (1998)
/// This version refers to blobs by their diameter (instead of "sigma").
/// This version can discard blobs which overlap with existing blobs.
/// (This is sometimes called "non-max suppression".)

template<class Scalar>
static
void
_BlobDogNM(int const image_size[3], //!<source image size
           Scalar const *const *const *aaafSource,   //!<source image
           Scalar const *const *const *aaafMask,     //!<ignore voxels where mask==0
           const vector<Scalar>& blob_diameters, //!<list of diameters to try (ordered)
           vector<array<Scalar,3> >& minima_crds, //!<store minima x,y,z coords here
           vector<array<Scalar,3> >& maxima_crds, //!<store maxima x,y,z coords here
           vector<Scalar>& minima_diameters, //!<corresponding radius for that minima
           vector<Scalar>& maxima_diameters, //!<corresponding radius for that maxima
           vector<Scalar>& minima_scores, //!<what was the blobs score?
           vector<Scalar>& maxima_scores, //!<(score = intensity after filtering)
           Scalar delta_sigma_over_sigma, //!<difference in Gauss widths parameter
           Scalar filter_truncate_ratio,     //!<how many sigma before truncating?
           Scalar filter_truncate_threshold, //!<decay in filter before truncating
           Scalar minima_threshold,       //!<discard unremarkable minima
           Scalar maxima_threshold,       //!<discard unremarkable maxima
           bool use_threshold_ratios=true, //!<threshold=ratio*best_score ?
           Scalar sep_ratio_thresh=1.0,          //!<minimum radial separation between blobs
           Scalar nonmax_max_overlap_large=1.0,  //!<maximum volume overlap with larger blob
           Scalar nonmax_max_overlap_small=1.0,  //!<maximum volume overlap with smaller blob
           ostream *pReportProgress = NULL,
           Scalar ****aaaafI = NULL, //!<preallocated memory for filtered images
           Scalar **aafI = NULL     //!<preallocated memory for filtered images
           )
{
  
  if (filter_truncate_ratio <= 0) {
    assert(filter_truncate_threshold > 0.0);
    //    filter_truncate_threshold = exp(-(1/2)*filter_truncate_ratio^2);
    //    -> filter_truncate_ratio^2 = -2*log(filter_truncate_threshold)
    filter_truncate_ratio = sqrt(-2*log(filter_truncate_threshold));
  }

  BlobDogNM(image_size,
            aaafSource,
            aaafMask,
            blob_diameters,
            minima_crds,
            maxima_crds,
            minima_diameters,
            maxima_diameters,
            minima_scores,
            maxima_scores,
            delta_sigma_over_sigma,
            filter_truncate_ratio,
            minima_threshold,
            maxima_threshold,
            use_threshold_ratios,
            sep_ratio_thresh,
            nonmax_max_overlap_large,
            nonmax_max_overlap_small,
            pReportProgress,
            aaaafI,
            aafI);

} //_BlobDogNM(...,filter_truncate_ratio,filter_truncate_threshold,...)




static void
HandleMinDistance(Settings settings,
                  MrcSimple &tomo_in,
                  MrcSimple &tomo_out,
                  MrcSimple &mask,
                  float voxel_width[3])
{
  fstream coords_file;
  coords_file.open(settings.in_coords_file_name.c_str(), ios::in);
  if (! coords_file)
    throw InputErr("Error: unable to open \""+
                   settings.in_coords_file_name +"\" for reading.\n");
  vector<array<int,3> > crds;
  while (coords_file) {
    double x, y, z;
    coords_file >> x;
    coords_file >> y;
    coords_file >> z;
    double ix, iy, iz;
    ix = floor((x / voxel_width[0]) + 0.5);
    iy = floor((y / voxel_width[1]) + 0.5);
    iz = floor((z / voxel_width[2]) + 0.5);
    array<int, 3> ixiyiz;
    ixiyiz[0] = ix;
    ixiyiz[1] = iy;
    ixiyiz[2] = iz;
    crds.push_back(ixiyiz);
  }

  cerr << " ------ calculating distance to points in "
       << settings.in_coords_file_name << " ------\n"
       << endl;

  // At some point I was trying to be as general as possible and allowed
  // for the possibility that voxels need not be cubes (same width x,y,z)
  // Now, I realize that allowing for this possibility would slow down some
  // calculations considerably, so I just assume cube-shaped voxels:
  float voxel_width_ = voxel_width[0];
  assert((voxel_width[0] == voxel_width[1]) &&
         (voxel_width[1] == voxel_width[2]));

  for (int iz=0; iz<tomo_out.header.nvoxels[2]; iz++) {
    cerr << "processing Z-plane " << iz+1 << " / "
         << tomo_out.header.nvoxels[2] << "\n";
    for (int iy=0; iy<tomo_out.header.nvoxels[1]; iy++) {
      for (int ix=0; ix<tomo_out.header.nvoxels[0]; ix++) {
        int rminsq_int = SQR(tomo_out.header.nvoxels[0] +
                             tomo_out.header.nvoxels[1] +
                             tomo_out.header.nvoxels[2]);

        if (mask.aaafI && (mask.aaafI[iz][iy][ix] == 0.0))
          continue;

        // Loop over all of the points and find the nearest one.
        // Calculate the distance to that point, and store that in the tomogram.
        // (Note: This is a terribly inefficient way to do this.
        //        A paint-bucket-like-tool using expanding spheres would
        //        be much faster and visit each voxel only once on average.)
        for (vector<array<int,3> >::const_iterator pxyz=crds.begin();
             pxyz != crds.end();
             pxyz++) {
          int rsqi=SQR(ix-(*pxyz)[0])+SQR(iy-(*pxyz)[1])+SQR(iz-(*pxyz)[2]);
          if (rsqi < rminsq_int)
            rminsq_int = rsqi;
        }
        float rmin = sqrt(rminsq_int * SQR(voxel_width_));
        tomo_out.aaafI[iz][iy][ix] = rmin;
      }
    }
  }
} //HandleMinDistance()














static void
HandleBlobsNonmaxSuppression(Settings settings,
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


  cerr << " ------ calculating distance and volume overlap in file: -------\n"
       << " " << settings.in_coords_file_name << "\n"
       << "\n";

  fstream coords_file;
  coords_file.open(settings.in_coords_file_name.c_str(), ios::in);
  if (! coords_file)
    throw InputErr("Error: unable to open \""+
                   settings.in_coords_file_name +"\" for reading.\n");

  bool custom_diameters = false;
  while (coords_file) {
    string strLine;
    getline(coords_file, strLine);
    if (strLine.size() == 0)
      continue;
    stringstream ssLine(strLine);
    double x, y, z;
    ssLine >> x;
    ssLine >> y;
    ssLine >> z;
    double ix, iy, iz;
    ix = floor((x / voxel_width[0]) + 0.5);
    iy = floor((y / voxel_width[1]) + 0.5);
    iz = floor((z / voxel_width[2]) + 0.5);
    array<float, 3> ixiyiz;
    ixiyiz[0] = ix;
    ixiyiz[1] = iy;
    ixiyiz[2] = iz;

    float diameter = -1.0;
    float score = settings.sphere_decals_foreground;
    if (ssLine) { // Does the file contain a 4th column? (the diameter)
      float _diameter;
      ssLine >> _diameter;
      // convert from physical distance to # of voxels:
      _diameter /= voxel_width[0];
      if (ssLine)
        diameter = _diameter;
      custom_diameters = true;
      if ((! ((settings.sphere_diameters_lower_bound <= diameter) &&
             (diameter <= settings.sphere_diameters_upper_bound)))
          &&
          (settings.sphere_diameters_lower_bound <=
           settings.sphere_diameters_upper_bound))
        // If the diameter lies outside of the desired range, ignore this blob
        continue; 
    }

    if (diameter < 0) { //If file does not contain a 4th column
      diameter = settings.sphere_decals_diameter;
      if (settings.sphere_decals_diameter < 0)
        diameter = 0.5;   // (sphere will be 1 voxel wide by default)
    }

    if (settings.sphere_decals_diameter >= 0) //override the diameter ?
      diameter = settings.sphere_decals_diameter;
    else
      diameter *= settings.sphere_decals_scale;

    if (ssLine) {
      float _score;
      ssLine >> _score;
      if (ssLine)
        score = _score;
    }

    if (! ((settings.score_lower_bound <= score) &&
           (score <= settings.score_upper_bound)))
      // If the score is not sufficiently high (or low), skip this blob
      continue; 

    crds.push_back(ixiyiz);
    diameters.push_back(diameter);
    scores.push_back(score);

  } //while (coords_file) {...


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





static void
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
                               voxel_width,
                               crds,
                               diameters,
                               scores);

  assert(crds.size() == diameters.size());
  assert(crds.size() == scores.size());

  // Sometimes the user wants to display the decal with a particular brightness
  // We store that brightness inside the "scores[]" array
  if (! settings.sphere_decals_foreground_use_score)
    for (size_t i=0; i < diameters.size(); i++)
      scores[i] = settings.sphere_decals_foreground;

  vector<float> shell_thicknesses(diameters.size());
  for (size_t i=0; i < diameters.size(); i++) {
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

  tomo_out = tomo_in; //copy the voxels from the original image to tomo_out

  VisualizeBlobs(tomo_out.header.nvoxels,
                 tomo_out.aaafI,
                 mask.aaafI,
                 crds,
                 diameters,
                 shell_thicknesses,
                 scores,
                 settings.sphere_decals_background,
                 settings.sphere_decals_background_scale,
                 settings.sphere_decals_foreground_norm);

} //HandleVisualizeBlobs()






static void
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
      throw InputErr("Error: unable to open \""+ settings.blob_minima_file_name +"\" for reading.\n");
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
      throw InputErr("Error: unable to open \""+ settings.blob_maxima_file_name +"\" for reading.\n");
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


    tomo_out = tomo_in; //copy the voxels from the original image to tomo_out

    VisualizeBlobs(tomo_out.header.nvoxels,
                   tomo_out.aaafI,
                   mask.aaafI,
                   display_crds_voxels,
                   display_diameters,
                   display_shell_thicknesses,
                   display_scores,
                   settings.sphere_decals_background,
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






static void
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
        if (! settings.use_dual_thresholds) {
          if (settings.out_threshold_01_a == settings.out_threshold_01_b)
            tomo_out.aaafI[iz][iy][ix] = ((tomo_out.aaafI[iz][iy][ix] >=
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



static void
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

  FindExtrema3D(tomo_in.header.nvoxels,
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
      throw InputErr("Error: unable to open \""+ settings.find_minima_file_name +"\" for reading.\n");
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
      throw InputErr("Error: unable to open \""+ settings.find_maxima_file_name +"\" for reading.\n");
    for (int i=0; i < maxima_crds_voxels.size(); i++)
      coords_file << maxima_crds_voxels[i][0] * voxel_width[0] << " "
                  << maxima_crds_voxels[i][1] * voxel_width[1] << " "
                  << maxima_crds_voxels[i][2] * voxel_width[2] << " "
                  //<< maxima_diameters[i] << " "
                  << maxima_nvoxels[i] << " "
                  << maxima_scores[i] << "\n";
  }
} //HandleExtrema()





static void
HandleLocalFluctuations(Settings settings,
                        MrcSimple &tomo_in,
                        MrcSimple &tomo_out,
                        MrcSimple &mask,
                        float voxel_width[3])
// Calculate the fluctuations of nearby voxel intensities
{

  // Filter weights, w_i:
  Filter3D<float, int>
    w = GenFilterGenGauss3D(settings.template_background_radius,
                            settings.template_background_exponent,
                            //template_profile.halfwidth);
                            settings.filter_truncate_ratio,
                            settings.filter_truncate_threshold,
                            static_cast<float*>(NULL),
                            static_cast<ostream*>(NULL));
  // GenFilterGenGauss3D() creates normalized gaussians with integral 1.
  // That's not what we want for the weights, w_i:
  // The weights w_i should be 1 in the viscinity we care about, and 0 outside
  // So, we can undo this normalization by dividing all the w_i weights
  // by their maximum value max(w_i)    (at the central peak)
  // This will mean the maximum w_i is 1, and the others decay to 0, as we want
  float wpeak = w.aaafH[0][0][0];
  w.MultiplyScalar(1.0 / wpeak);


  cerr << " ------ Calculating the average of nearby voxels: ------\n";
  // P = original image (after subtracting average nearby intensities):

  MrcSimple P = tomo_in; //this will take care of allocating the array

  float template_background_sigma[3];
  for (int d = 0; d < 3; d++)
    template_background_sigma[d] = (settings.template_background_radius[d]
                                    /
                                    sqrt(3.0));

  // First, let's calculate the weighted average voxel intensity in the 
  // source image
  if (settings.template_background_exponent == 2.0) {

    // then do it the fast way with seperable (ordinary) Gaussian filters
    ApplyGauss3D(tomo_in.header.nvoxels,
                 tomo_in.aaafI,
                 P.aaafI,    // <-- save result here
                 mask.aaafI,
                 template_background_sigma,
                 settings.filter_truncate_ratio,
                 settings.filter_truncate_threshold,
                 true,
                 &cerr);
  }
  else {
    w.Apply(tomo_in.header.nvoxels,
            tomo_in.aaafI,
            P.aaafI,    // <-- save result here
            mask.aaafI,
            true,
            &cerr);
  }

  // Subtract the average value from the image intensity, and store in P:
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
        P.aaafI[iz][iy][ix] = tomo_in.aaafI[iz][iy][ix]-P.aaafI[iz][iy][ix];

  // Calculate <P_, P_>

  // Because memory is so precious, we must reuse arrays whenever we can.
  // So we will now use "P" (which used to denote voxel intensity)
  // to store the square of intensity instead
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
        P.aaafI[iz][iy][ix] *= P.aaafI[iz][iy][ix];

  cerr << "\n"
    " ------ Calculating fluctuations around that average ------\n" << endl;

  //   Original code:  Commented out because it requires too much memory:
  //MrcSimple P_dot_P = tomo_in; //(this takes care of allocating the array)
  //   Memory efficient, ugly code:
  // We're done with tomo_in.  Re-use this array to store P_dot_P temporarilly,
  // but give it a new name, to make the code slightly less confusing to read:
  float ***P_dot_P_aaafI = tomo_in.aaafI;

  if (settings.template_background_exponent == 2.0) {
    // then do it the fast way with seperable (ordinary) Gaussian filters
    ApplyGauss3D(tomo_in.header.nvoxels,
                 P.aaafI,
                 //P_dot_P.aaafI,
                 P_dot_P_aaafI,    // <-- save result here
                 mask.aaafI,
                 template_background_sigma,
                 settings.filter_truncate_ratio,
                 settings.filter_truncate_threshold,
                 true,
                 &cerr);
  }
  else {
    w.Apply(tomo_in.header.nvoxels,
            P.aaafI,
            //P_dot_P.aaafI,
            P_dot_P_aaafI,  // <-- store result here
            mask.aaafI,
            false,
            &cerr);
  }

  // Now calculate "rms" (sqrt(variance))
  // Save the result in tomo_out

  // Rrite out RMS variance from the average of nearby voxels:
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
        //     variance = <P_, P_>

        //float variance = (P_dot_P.aaafI[iz][iy][ix]
        float variance = P_dot_P_aaafI[iz][iy][ix];

        //Optional:
        //Compensate for dividing by w.aaafH[][][] by "wpeak" earlier.
        //This enables us to interpret variance as RMSE, (a.k.a. root-
        //mean-squared-error.  The formula above only calculates the
        //"mean" if w_i are normalized, which they were before we divided
        //them all by wpeak.  (Also: This insures that this fast method is
        //equivalent to the "SlowDebug" method commented out above.)
        //(Note: It is important for other reasons that w_i (w.aaafH[][][])
        //       were not normalized until this point.)

        variance *= wpeak;

        if (variance < 0.0)
          variance = 0.0;

        tomo_out.aaafI[iz][iy][ix] = sqrt(variance);
      } //for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
    } //for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
  } //for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)

  tomo_out.FindMinMaxMean();
  
} //HandleLocalFluctuations()







template <typename T>
static int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}



template<class Scalar>
static void
WriteBNPTSfile(string filename,
               vector<array<Scalar,3> > coords,
               vector<array<Scalar,3> > norms)
{
  assert(coords.size() == norms.size());
  size_t n = coords.size();
  fstream bnpts_file;
  bnpts_file.open(filename.c_str(), ios::out | ios::binary);
  for (size_t i=0; i < n; i++) {
    float xyz[3];
    xyz[0] = coords[i][0];  //(convert from "Scalar" to float)
    xyz[1] = coords[i][1];
    xyz[2] = coords[i][2];
    float norm[3];
    norm[0] = norms[i][0];  //(convert from "Scalar" to float)
    norm[1] = norms[i][1];
    norm[2] = norms[i][2];
    bnpts_file.write((char*)&(xyz[0]), sizeof(float));
    bnpts_file.write((char*)&(xyz[1]), sizeof(float));
    bnpts_file.write((char*)&(xyz[2]), sizeof(float));
    bnpts_file.write((char*)&(norm[0]), sizeof(float));
    bnpts_file.write((char*)&(norm[1]), sizeof(float));
    bnpts_file.write((char*)&(norm[2]), sizeof(float));
  }
  bnpts_file.close();
}


template<class Scalar, class VectorContainer>
static void
WriteOrientedPointCloud(string pointcloud_file_name,
                        const int image_size[3],
                        Scalar const *const *const *aaafImage,
                        CompactMultiChannelImage3D<Scalar> &hessian,
                        Scalar threshold)
{
  assert(aaafImage);
  vector<array<Scalar,3> > coords;
  vector<array<Scalar,3> > norms;
  for (int iz=0; iz < image_size[2]; iz++) {
    for (int iy=0; iy < image_size[1]; iy++) {
      for (int ix=0; ix < image_size[0]; ix++) {
        Scalar metric = aaafImage[iz][iy][ix];
        if ((abs(metric) >= abs(threshold)) &&
            (metric*threshold > 0.0)) //same sign
        {
          // check if this is one of the voxels we want to consider
          if (! hessian.aaaafI[iz][iy][ix])
            // (if not, then this voxel lies outside the mask, so skip it)
            continue;

          array<Scalar,3> xyz;
          xyz[0] = ix;
          xyz[1] = iy;
          xyz[2] = iz;
          coords.push_back(xyz);

          array<Scalar,3> norm;
          Scalar quat[4];
          quat[0] = hessian.aaaafI[iz][iy][ix][3];
          quat[1] = hessian.aaaafI[iz][iy][ix][4];
          quat[2] = hessian.aaaafI[iz][iy][ix][5];
          quat[3] = hessian.aaaafI[iz][iy][ix][6];

          Scalar eigenvects[3][3];
          Quaternion2Matrix(quat, eigenvects);
          norm[0] = eigenvects[0][0];
          norm[1] = eigenvects[0][1];
          norm[2] = eigenvects[0][2];
          norms.push_back(norm);
        }
      }
    }
  }

  WriteBNPTSfile(pointcloud_file_name,
                 coords,
                 norms);

} //WriteOrientedPointCloud()




static void
HandleRidgeDetectorPlanar(Settings settings,
                          MrcSimple &tomo_in,
                          MrcSimple &tomo_out,
                          MrcSimple &mask,
                          float voxel_width[3])
{
  cerr << "filter_type = planar ridge detector\n";

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
  CompactMultiChannelImage3D<float> c_hessian(6);
  c_hessian.Resize(tomo_in.header.nvoxels, mask.aaafI, &cerr);
  

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
    ApplyGauss3D(tomo_in.header.nvoxels,
                 tomo_in.aaafI,
                 tomo_background.aaafI,
                 mask.aaafI,
                 settings.width_b[0],
                 truncate_halfwidth,
                 true,
                 &cerr);

    truncate_halfwidth = floor(settings.width_a[0] *
                               settings.filter_truncate_ratio);
    ApplyGauss3D(tomo_in.header.nvoxels,
                 tomo_in.aaafI,
                 tomo_out.aaafI,
                 mask.aaafI,
                 settings.width_a[0],
                 truncate_halfwidth,
                 true,
                 &cerr);
  }

  // Using this blurred image, calculate the 2nd derivative matrix everywhere:

  CalcHessian3D(tomo_in.header.nvoxels,
                tomo_in.aaafI,
                aaaafGradient,
                c_hessian.aaaafI,
                mask.aaafI,
                sigma,
                settings.filter_truncate_ratio,
                &cerr);

  // Now calculate the eigenvalues and eigenvectors of the hessian matrix
  // located at each voxel.
  // (This converts an array of 6 numbers representing the non-redundant
  //  entries from the symmetrix 3x3 matrix, into 3 eigenvalues,
  //  and 3 "Shoemake" coordinates, from which eigenvectors can be calculated)
  // To save space, I chose to store the result in the same array.

  DiagonalizeHessianImage3D(tomo_in.header.nvoxels,
                            c_hessian.aaaafI,
                            c_hessian.aaaafI,
                            mask.aaafI,
                            selfadjoint_eigen3::DECREASING_EIVALS,
                            &cerr);


  // We need to store the direction of the most important eigenvector
  // somewhere.  To save space, why not store it in the aaaafGradient
  // array?  (At this point, we are no longer using it).

  array<float, 3> ***aaaafStickDirection = aaaafGradient;

  // This just effectively changes the name
  // from "aaaafGradient" to "aaaafStickDirection".
  // The purpose of the name change is to make it easier to read the code later.


  // Optional: store the saliency (score) of each voxel in tomo_out.aaafI
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
        tomo_out.aaafI[iz][iy][ix] = 0.0;

  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {

        if (! c_hessian.aaaafI[iz][iy][ix]) //ignore voxels that are either
          continue;                         //in the mask or on the boundary

        float eivals[3];
        eivals[0] = c_hessian.aaaafI[iz][iy][ix][0]; //maximum eigenvalue
        eivals[1] = c_hessian.aaaafI[iz][iy][ix][1];
        eivals[2] = c_hessian.aaaafI[iz][iy][ix][2]; //minimum eigenvalue
        // To save space the eigenvectors were stored as a "Shoemake" 
        // coordinates (similar to quaternions) instead of a 3x3 matrix.
        // So we must unpack the eigenvectors.
        float shoemake[3]; // <- the eigenvectors stored in "Shoemake" format
        shoemake[0]       = c_hessian.aaaafI[iz][iy][ix][3];
        shoemake[1]       = c_hessian.aaaafI[iz][iy][ix][4];
        shoemake[2]       = c_hessian.aaaafI[iz][iy][ix][5];

        // Now lets extract the eigenvectors:
        float eivects[3][3];
        Shoemake2Matrix(shoemake, eivects); //convert to 3x3 matrix

        // REMOVE THIS CRUFT ?
        //float grad[3];
        //grad[0] = aaaafGradient[iz][iy][ix][0];
        //grad[1] = aaaafGradient[iz][iy][ix][1];
        //grad[2] = aaaafGradient[iz][iy][ix][2];

        float score;

        // DEBUG: REMOVE THE NEXT IF STATMENT AFTER DEBUGGING IS FINISHED
        #ifndef NDEBUG
        if ((ix==tomo_in.header.nvoxels[0]/2) &&
            (iy==tomo_in.header.nvoxels[1]/2) &&
            (iz==tomo_in.header.nvoxels[2]/2))
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

        score = ScoreHessianPlanar(c_hessian.aaaafI[iz][iy][ix],
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


        if (settings.planar_tv_sigma > 0.0) {
          aaaafStickDirection[iz][iy][ix][0] = eivects[0][0];
          aaaafStickDirection[iz][iy][ix][1] = eivects[0][1];
          aaaafStickDirection[iz][iy][ix][2] = eivects[0][2];
        }


      } //for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
    } //for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
  } //for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)





  if (settings.planar_tv_sigma > 0.0) {
    assert(settings.filter_truncate_ratio > 0);

    // Use thresholding to reduce the number of voxels that we have to consider
    // Later, voxels with 0 saliency will be ignored.
    for (int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
      for (int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
        for (int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
          if (tomo_out.aaafI[iz][iy][ix] < settings.planar_hessian_score_threshold)
            tomo_out.aaafI[iz][iy][ix] = 0.0;
        }
      }
    }

    TV3D<float, int, array<float,3>, float* >
      tv(settings.planar_tv_sigma,
         settings.planar_tv_exponent,
         settings.planar_tv_truncate_ratio);

    tv.TVDenseStick(tomo_in.header.nvoxels,
                    tomo_out.aaafI,
                    aaaafStickDirection,
                    c_hessian.aaaafI,
                    mask.aaafI,
                    mask.aaafI,
                    false,  // (we want to detect surfaces not curves)
                    //settings.planar_hessian_score_threshold,
                    true,   // (do normalize near rectangular image bounaries)
                    &cerr);

    for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
      for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
        for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
          float score = ScoreTensorPlanar(c_hessian.aaaafI[iz][iy][ix]);
          float peak_height = 1.0;
          if (tomo_background.aaafI)
            peak_height = (tomo_in.aaafI[iz][iy][ix] -
                           tomo_background.aaafI[iz][iy][ix]);
          score *= peak_height;
          tomo_out.aaafI[iz][iy][ix] = score;
        }
      }
    }
  } // if (settings.planar_tv_sigma > 0.0)



  ////Did the user ask us to generate any output files?
  //if ((settings.out_normals_fname == "") &&
  //    (settings.out_file_name != ""))
  //  settings.out_normals_fname =
  //    settings.out_file_name + string(".bnpts");


  #if 0
  //Did the user ask us to generate output files containing surface orientation?
  if (settings.out_normals_fname != "")
    WriteOrientedPointCloud(settings.out_normals_fname,
                            tomo_out.header.nvoxels,
                            tomo_out.aaafI,
                            c_hessian,
                            settings.planar_tv_score_threshold);
  #endif

  Dealloc3D(tomo_in.header.nvoxels,
            &(aafGradient),
            &(aaaafGradient));
} //HandleRidgeDetectorPlanar()




static void
HandleWatershed(Settings settings,
                MrcSimple &tomo_in,
                MrcSimple &tomo_out,
                MrcSimple &mask,
                float voxel_width[3])
{
  vector<array<int, 3> > extrema_crds;
  vector<float> extrema_scores;

  Watershed3D(tomo_in.header.nvoxels,
              tomo_in.aaafI,
              tomo_out.aaafI,
              mask.aaafI,
              settings.watershed_threshold,
              settings.watershed_use_minima,
              settings.neighbor_connectivity,
              settings.watershed_show_boundaries,
              &extrema_crds,
              &extrema_scores,
              &cerr);

  // Did the user supply a mask?
  // WatershedMeyers3D() intentionally does not modify voxels which lie 
  // outside the mask.  These voxels will have random undefined values 
  // unless we assign them manually.  We do this below.
  float UNDEFINED = -1;
  for (int iz = 0; iz < tomo_in.header.nvoxels[2]; iz++)
    for (int iy = 0; iy < tomo_in.header.nvoxels[1]; iy++)
      for (int ix = 0; ix < tomo_in.header.nvoxels[0]; ix++)
        if (mask.aaafI && (mask.aaafI[iz][iy][ix] == 0.0))
          tomo_out.aaafI[iz][iy][ix] = UNDEFINED;
} //HandleWatershed()






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
        cerr << "   (Using " << nthr << " threads (cpu cores).  You can change this using the \"-np n\"\n"
             << "    argument, or by setting the OMP_NUM_THREADS environment variable.)" << endl;
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
    }

    // ---- mask ----

    // Optional: if there is a "mask", read that too
    MrcSimple mask;
    if (settings.mask_file_name != "") {
      cerr << "Reading mask \""<<settings.mask_file_name<<"\"" << endl;
      mask.Read(settings.mask_file_name, false);
      if ((mask.header.nvoxels[0] != tomo_in.header.nvoxels[0]) ||
          (mask.header.nvoxels[1] != tomo_in.header.nvoxels[1]) ||
          (mask.header.nvoxels[2] != tomo_in.header.nvoxels[2]))
        throw InputErr("Error: The size of the mask image does not match the size of the input image.\n");
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

    // ---- make an array that will store the new tomogram we will create ----

    cerr << "allocating space for new tomogram..." << endl;
    MrcSimple tomo_out = tomo_in; //this will take care of allocating the array

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

    if ((voxel_width[0] != voxel_width[1]) ||
        (voxel_width[1] != voxel_width[2])) {
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
      throw InputErr(err_msg.str().c_str());
    }


    settings.planar_tv_sigma /= voxel_width[0];
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
    for (int ir = 0; ir < settings.blob_diameters.size(); ir++)
      settings.blob_diameters[ir] /= voxel_width[0];

    settings.sphere_decals_diameter /= voxel_width[0];
    if (! settings.sphere_decals_shell_thickness_is_ratio)
      settings.sphere_decals_shell_thickness /= voxel_width[0];



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



    else if (settings.filter_type == Settings::RIDGE_PLANAR) {

      // find planar ridges (ie membranes or wide tubes)
      HandleRidgeDetectorPlanar(settings, tomo_in, tomo_out, mask, voxel_width);

    }


    else if (settings.filter_type == Settings::WATERSHED) {

      // perform watershed segmentation
      HandleWatershed(settings, tomo_in, tomo_out, mask, voxel_width);

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




    else if (settings.filter_type == Settings::BOOTSTRAP_DOGG) {

      #ifdef DISABLE_BOOSTRAPPING
      HandleBootstrappDogg(settings, tomo_in, tomo_out, mask, voxel_width);
      #endif //#ifndef DISABLE_BOOTSTRAPPING

    } //else if (settings.filter_type == Settings::BOOTSTRAP_DOGG) {


    // ----- distance_filter -----

    else if (settings.filter_type == settings.MIN_DISTANCE) {

      HandleMinDistance(settings, tomo_in, tomo_out, mask, voxel_width);

    }

    // ----- sphere_decals_filter -----

    else if (settings.filter_type == settings.SPHERE_DECALS) {

      HandleVisualizeBlobs(settings, tomo_in, tomo_out, mask, voxel_width);

    }

    else if (settings.filter_type == settings.SPHERE_NONMAX_SUPPRESSION) {

      vector<array<float,3> > crds;
      vector<float> diameters;
      vector<float> scores;

      HandleBlobsNonmaxSuppression(settings,
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
    
    if (settings.use_thresholds) {

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

    tomo_out.FindMinMaxMean();


    // --- Rescale so that the lowest, highest voxels have density 0 and 1? ---

    if (settings.rescale01_out)
      tomo_out.Rescale01(mask.aaafI);








    // ------ Write the file ------
    if (settings.out_file_name != "") {
      cerr << "writing tomogram (in 32-bit float mode)" << endl;
      tomo_out.Write(settings.out_file_name);
      // (You can also use "file_stream << tomo_out;")
    }

  } // try {

  catch (InputErr& e) {
    cerr << "\n" << e.what() << endl;
    exit(1);
  }
} //main()


