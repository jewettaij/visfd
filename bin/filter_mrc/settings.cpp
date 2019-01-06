#include <sstream>
#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>
using namespace std;

#ifndef DISABLE_OPENMP
#include <omp.h>       // (OpenMP-specific)
#endif

#include <err_report.hpp>
#include "settings.hpp"


template<class RealNum>
inline RealNum MAX(RealNum x, RealNum y) { return ((x<y) ? y : x); }


Settings::Settings() {
  // Default settings
  in_file_name = "";
  invert_output = false;
  out_file_name = "";
  rescale01_in = false;
  rescale01_out = false;
  mask_file_name = "";
  mask_select = 1;
  use_mask_select = false;
  mask_out = 0.0;
  use_mask_out = true;
  voxel_width = 0.0;  //How many Angstroms per voxel? (if 0 then read from file)
  voxel_width_divide_by_10 = false; //Use nm instead of Angstroms?
  filter_type = NONE;

  // The code is a little bit confusing because I tried to cram as much
  // functionality into the smallest number of parameters possible:
  // The "difference of gaussians" (dog) filter used here requires
  //  8 parameters: a_x, a_y, a_z, b_x, b_y, b_z, m, n
  // For details, see the comment at the beginning of "settings.h".
  //
  // Default values of a_x, a_y, a_z, b_x, b_y, b_z, m, n are:

  width_a[0] = 70.0;  // a_x =70.0 gaussian width in x direction(physical units)
  width_a[1] = 70.0;  // a_y = 70.0  gaussian width in y direction
  width_a[2] = 110.0; // a_z = 110.0  gaussian width in z direction
  width_b[0] = 90.0;  // b_x = 90.0  2nd (outer) gaussian width in x direction
  width_b[1] = 90.0;  // b_y = 905.0  gaussian width in y direction
  width_b[2] = 130.0; // b_z = 130.0  gaussian width in z direction(<0 disables)
  m_exp = 2.0;       // exponent in generalized Gaussian formula (width a)
  n_exp = 2.0;       // exponent in generalized Gaussian formula (width b)

  //filter_truncate_halfwidth[0] = -1; // width of the filter used (in voxels)
  //filter_truncate_halfwidth[1] = -1; // (This is the size of the domain of the function
  //filter_truncate_halfwidth[2] = -1; //  which will be convolved with the image.
  //                          //  "-1" means unspecified.)

  filter_truncate_threshold=0.02;    //Filter intensity decay value before giving up
                            //When the filter strength is less than this value
                            //we ignore it. For difference-of-gaussian filters
                            //we choose the gaussian with the wider width. This
                            //parameter overrides other window-width settings.
                            //(Setting it to a number < 0 disables it.)

  filter_truncate_ratio = -1.0; //When averaging/filtering consider nearby
                            //voxels up to a distance of this many sigma away
                            //(Setting this to a number < 0 disables it.)

  dogsf_width[0] = 0.0;
  dogsf_width[1] = 0.0;
  dogsf_width[2] = 0.0;
  delta_sigma_over_sigma = 0.02;

  find_minima = false;
  find_maxima = false;
  find_minima_file_name = "";
  find_maxima_file_name = "";
  //REMOVE THIS CRUFT
  //find_extrema_occlusion_ratio = 1.0;
  in_coords_file_name = "";
  out_coords_file_name = "";
  sphere_decals_diameter = -1.0;
  sphere_decals_shell_thickness = -1.0;
  sphere_decals_foreground = 1.0;
  sphere_decals_background = 0.0;
  sphere_decals_background_scale = 0.25; //default dimming of original image to emphasize spheres
  sphere_decals_foreground_use_score = true;
  //sphere_decals_background_use_orig = true;
  sphere_decals_foreground_norm = false;
  sphere_decals_scale = 1.0;
  sphere_decals_shell_thickness_is_ratio = true;
  sphere_decals_shell_thickness = 0.05;
  sphere_decals_shell_thickness_min = 1.0;
  score_lower_bound = -std::numeric_limits<float>::infinity();
  score_upper_bound = std::numeric_limits<float>::infinity();
  score_bounds_are_ratios = true;
  sphere_diameters_lower_bound = -std::numeric_limits<float>::infinity();
  sphere_diameters_upper_bound = std::numeric_limits<float>::infinity();
  blob_width_multiplier = 1.0;
  nonmax_min_radial_separation_ratio = 0.0;
  nonmax_max_volume_overlap_small = 1.0;
  nonmax_max_volume_overlap_large = 1.0;

  watershed_use_minima = true;
  watershed_threshold = std::numeric_limits<float>::infinity();

  out_normals_fname = "";
  planar_hessian_score_threshold = 0.0;
  planar_tv_score_threshold = 0.0;
  planar_tv_sigma = 0.0;
  planar_tv_exponent = 4;
  planar_tv_num_iters = 0;
  planar_tv_truncate_ratio = sqrt(2.0);

  use_thresholds = false;
  use_dual_thresholds = false;
  out_threshold_01_a = 1.0;
  out_threshold_01_b = 1.0;
  out_threshold_10_a = 1.0;
  out_threshold_10_b = 1.0;
  out_thresh2_use_clipping = false;
  out_thresh_a_value = 0.0;
  out_thresh_b_value = 1.0;
  //missing_wedge_min[0]=-90.0; // By default, the "missing wedge" includes all
  //missing_wedge_max[0]=90.0;  // orientations around the Y axis which lie
  //missing_wedge_min[1]=-30.0; // between -30 and +30 degrees (relative to the
  //missing_wedge_max[1]=+30.0; // Z axis), independent of X-axis orientation.

  #ifndef DISABLE_BOOTSTRAPPING
  bs_ntests = 0;                   //disable
  bs_threshold = 0.0;              //disable
  float bs_threshold_sign = 0.0;   //disable
  bs_scramble_radius[0] = -1;      //impossible value
  bs_scramble_radius[1] = -1;      //impossible value
  bs_scramble_radius[2] = -1;      //impossible value
  f_bs_scramble_radius[0] = -1.0;  //impossible value
  f_bs_scramble_radius[1] = -1.0;  //impossible value
  f_bs_scramble_radius[2] = -1.0;  //impossible value
  bs_random_seed = 1;
  #endif //#ifndef DISABLE_BOOTSTRAPPING

  #ifndef DISABLE_TEMPLATE_MATCHING
  template_background_radius[0] -1.0;         //impossible value
  template_background_radius[1] -1.0;         //impossible value
  template_background_radius[2] -1.0;         //impossible value
  template_background_exponent = 2.0;         //default value
  //template_compare_exponent = 2.0;            //default value
  #endif //#ifndef DISABLE_TEMPLATE_MATCHING

} //Settings::Settings()



void Settings::ParseArgs(int argc, char **argv) {
  vector<string> vArgs;
  ConvertArgvToVectorOfStrings(argc, argv, vArgs);
  ParseArgs(vArgs);
}



void Settings::ConvertArgvToVectorOfStrings(int argc,
                                            char **argv,
                                            vector<string>& dest)
{
  dest.resize(argc);
  for (int i=0; i < argc; ++i)
    dest[i].assign(argv[i]);
}



void
Settings::ParseArgs(vector<string>& vArgs)
{
  bool exponents_set_by_user = false;
  bool delta_set_by_user = false;
  for (int i=1; i < vArgs.size(); ++i)
  {

    int num_arguments_deleted = 0;

    if ((vArgs[i] == "-in") || (vArgs[i] == "-i"))
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a file name.\n");
      in_file_name = vArgs[i+1];

      num_arguments_deleted = 2;

    } // if ((vArgs[i] == "-in") || (vArgs[i] == "-i"))


    else if ((vArgs[i] == "-out") || (vArgs[i] == "-o"))
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a file name.\n");
      out_file_name = vArgs[i+1];

      num_arguments_deleted = 2;

    } // if ((vArgs[i] == "-out") || (vArgs[i] == "-o"))


    else if (vArgs[i] == "-mask")
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a file name.\n");
      mask_file_name = vArgs[i+1];
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-mask")


    else if (vArgs[i] == "-mask-select")
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        use_mask_select = true;
        mask_select = stoi(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by an integer.\n");
      }
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-mask-select")


    else if (vArgs[i] == "-mask-out")
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        use_mask_out = true;
        mask_out = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number.\n");
      }
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-mask-out")


    else if (vArgs[i] == "-w")
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        voxel_width = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by voxel width.\n");
      }
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-mask")


    else if ((vArgs[i] == "-a2nm") || (vArgs[i] == "-ang-to-nm"))
    {
      voxel_width_divide_by_10 = true;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-a2nm")


    else if ((vArgs[i] == "-no-rescale") || (vArgs[i] == "-norescale"))
    {
      rescale01_in = false;
      rescale01_out = false;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-norescale")


    else if (vArgs[i] == "-rescale")
    {
      if (filter_type == NONE) {
        rescale01_in = true;
        out_threshold_01_a = 0.0;
        out_threshold_01_b = 1.0;
      }
      else {
        rescale01_out = true;
      }
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-rescale")


    else if (vArgs[i] == "-invert")
    {
      invert_output = true;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-invert")


    else if ((vArgs[i] == "-gauss-aniso") || (vArgs[i] == "-ggauss-aniso"))
    {
      try {
        if ((i+3 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-'))
          throw invalid_argument("");
        width_a[0] = stof(vArgs[i+1]);
        width_a[1] = stof(vArgs[i+2]);
        width_a[2] = stof(vArgs[i+3]);
        width_b[0] = -1.0;
        width_b[1] = -1.0;
        width_b[2] = -1.0;
        filter_type = GAUSS;
        if (vArgs[i] == "-ggauss-aniso")
          filter_type = GGAUSS;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 3 positive numbers:\n"
                       " s_x  s_y  s_z\n"
                       " the Gaussian widths in the X, Y, and Z direction.)\n");
      }
      num_arguments_deleted = 4;
    } //if (vArgs[i] == "-gauss-aniso")



    else if ((vArgs[i] == "-gauss") || (vArgs[i] == "-ggauss"))
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        width_a[0] = stof(vArgs[i+1]);
        width_a[1] = width_a[0];
        width_a[2] = width_a[0];
        width_b[0] = -1.0;
        width_b[1] = -1.0;
        width_b[2] = -1.0;
        filter_type = GAUSS;
        if (vArgs[i] == "-ggauss")
          filter_type = GGAUSS;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a positive number (\"s\"),\n"
                       " the Gaussian width\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-gauss")



    else if ((vArgs[i] == "-dog-aniso") ||
             (vArgs[i] == "-dogg-aniso"))
    {
      try {
        if ((i+6 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-') ||
            (vArgs[i+4] == "") || (vArgs[i+4][0] == '-') ||
            (vArgs[i+5] == "") || (vArgs[i+5][0] == '-') ||
            (vArgs[i+6] == "") || (vArgs[i+6][0] == '-'))
          throw invalid_argument("");
        width_a[0] = stof(vArgs[i+1]);
        width_a[1] = stof(vArgs[i+2]);
        width_a[2] = stof(vArgs[i+3]);
        width_b[0] = stof(vArgs[i+4]);
        width_b[1] = stof(vArgs[i+5]);
        width_b[2] = stof(vArgs[i+6]);
        //if (width_b[0] <= width_a[0])
        //  throw InputErr("Error: The two arguments following " + vArgs[i] + " must be\n"
        //                 " increasing.  (Ie., the 2nd argument must be > 1st argument.)\n");
        filter_type = DOG;
        if (vArgs[i] == "-dogg-aniso")
          filter_type = DOGG;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 6 positive numbers.\n");
      }
      num_arguments_deleted = 7;
    } //if (vArgs[i] == "-dog-aniso")



    else if ((vArgs[i] == "-dog") ||
             (vArgs[i] == "-dogg"))
    {
      try {
        if ((i+2 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-'))
          throw invalid_argument("");
        width_a[0] = stof(vArgs[i+1]);
        width_a[1] = width_a[0];
        width_a[2] = width_a[0];
        width_b[0] = stof(vArgs[i+2]);
        width_b[1] = width_b[0];
        width_b[2] = width_b[0];
        //if (width_b[0] <= width_a[0])
        //  throw InputErr("Error: The two arguments following " + vArgs[i] + " must be\n"
        //                 " increasing.  (Ie., the 2nd argument must be > 1st argument.)\n");
        filter_type = DOG;
        if (vArgs[i] == "-dogg")
          filter_type = DOGG;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 2 positive numbers.\n");
      }
      num_arguments_deleted = 3;
    } //if (vArgs[i] == "-dog")



    else if (vArgs[i] == "-doggxy-aniso")
    {
      try {
        if ((i+5 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-') ||
            (vArgs[i+4] == "") || (vArgs[i+4][0] == '-') ||
            (vArgs[i+5] == "") || (vArgs[i+5][0] == '-'))
          throw invalid_argument("");
        width_a[0] = stof(vArgs[i+1]);
        width_a[1] = stof(vArgs[i+2]);
        width_b[0] = stof(vArgs[i+3]);
        width_b[1] = stof(vArgs[i+4]);
        //The "-doggxy" filter is a Difference-of-Generalized-Gaussians 
        //in the X,Y  directions, multiplied by an ordinary Gaussian in 
        //the Z direction.
        width_a[2] = stof(vArgs[i+5]);
        width_b[2] = -1.0;   //(This disables the second (negative) Gaussian)
        //if (width_b[0] <= width_a[0])
        //  throw InputErr("Error: The two arguments following " + vArgs[i] + " must be\n"
        //                 " increasing.  (Ie., the 2nd argument must be > 1st argument.)\n");
        filter_type = DOGGXY;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 5 positive numbers:\n"
                       " a_x  a_y  b_x  b_y  a_z\n"
                       " (I.E., the \"A\" and \"B\" Gaussian widths in X and Y directions,\n"
                       "  followed by the Gaussian width in the Z direction.)\n");
      }
      num_arguments_deleted = 6;
    } //if (vArgs[i] == "-doggxy-aniso")


    else if (vArgs[i] == "-doggxy")
    {
      try {
        if ((i+5 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-'))
          throw invalid_argument("");
        width_a[0] = stof(vArgs[i+1]);
        width_a[1] = width_a[0];
        width_b[0] = stof(vArgs[i+2]);
        width_b[1] = width_b[0];
        //The "-doggxy" filter is a Difference-of-Generalized-Gaussians 
        //in the X,Y  directions, multiplied by an ordinary Gaussian in 
        //the Z direction.
        width_a[2] = stof(vArgs[i+3]);
        width_b[2] = -1.0;   //(This disables the second (negative) Gaussian)
        //if (width_b[0] <= width_a[0])
        //  throw InputErr("Error: The two arguments following " + vArgs[i] + " must be\n"
        //                 " increasing.  (Ie., the 2nd argument must be > 1st argument.)\n");
        filter_type = DOGGXY;
      }
      catch (invalid_argument& exc) {
          throw InputErr("Error: The " + vArgs[i] + 
                         " argument must be followed by 3 positive numbers:\n"
                         " a_xy  b_xy  a_z\n"
                         " (I.E., the \"A\" and \"-B\" Gaussian widths in the XY plane,\n"
                         "  followed by the Gaussian width in the Z direction.)\n");
      }
      num_arguments_deleted = 4;
    } //if (vArgs[i] == "-doggxy")



    else if ((vArgs[i] == "-blob-single") ||
             (vArgs[i] == "-blob1") ||
             (vArgs[i] == "-blob1-s") || (vArgs[i] == "-blob-s1") ||
             (vArgs[i] == "-blob1-r") || (vArgs[i] == "-blob-r1") ||
             (vArgs[i] == "-blob-1d") || (vArgs[i] == "-blob-d1"))
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        float blob_width_multiplier = 1.0;
        blob_width_multiplier = 1.0;
        if ((vArgs[i] == "-blob1-s") || (vArgs[i] == "-blob-s1"))
          // (for a solid uniform 3-D sphere)
          blob_width_multiplier = 2.0*sqrt(3.0);
          //REMOVE THIS CRUFT
          //blob_width_multiplier = 1.0/sqrt(3.0);
        if ((vArgs[i] == "-blob1-r") || (vArgs[i] == "-blob-r1"))
          // (for a solid uniform 3-D sphere)
          blob_width_multiplier = 2.0;
          //REMOVE THIS CRUFT
          //blob_width_multiplier = 1.0/sqrt(3.0);
        else if ((vArgs[i] == "-blob-single") ||
                 (vArgs[i] == "-blob1") ||
                 (vArgs[i] == "-blob1-d") || (vArgs[i] == "-blob-d1"))
          // (for a solid uniform 3-D sphere)
          blob_width_multiplier = 1.0;
          //REMOVE THIS CRUFT
          //blob_width_multiplier = 1.0/(2.0*sqrt(3.0));
        dogsf_width[0] = stof(vArgs[i+1]) * blob_width_multiplier;
        dogsf_width[1] = dogsf_width[0];
        dogsf_width[2] = dogsf_width[0];
        m_exp = 2.0;
        n_exp = 2.0;
        filter_type = DOG_SCALE_FREE;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a positive number:\n"
                       "       the (approximate) width of the objects of interest.\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-blob-single")


    else if (vArgs[i] == "-blob1-aniso")
    {
      try {
        if ((i+3 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-'))
          throw invalid_argument("");
        dogsf_width[0] = stof(vArgs[i+1]);
        dogsf_width[1] = stof(vArgs[i+2]);
        dogsf_width[2] = stof(vArgs[i+3]);
        m_exp = 2.0;
        n_exp = 2.0;
        filter_type = DOG_SCALE_FREE;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 3 positive numbers:\n"
                       "       the approximate widths of the objects of interest in the XYZ directions\n"
                       "       (Due to the missing-wedge artifact they might differ.\n");
      }
      num_arguments_deleted = 4;
    } //if (vArgs[i] == "-blob1-aniso")



    else if ((vArgs[i] == "-spheres-nonmax-overlap") ||
             (vArgs[i] == "-max-volume-overlap") || 
             (vArgs[i] == "-max-overlap"))
    {
      try {
        if (i+1 >= vArgs.size())
          throw invalid_argument("");
        nonmax_max_volume_overlap_large = stof(vArgs[i+1]);
        nonmax_min_radial_separation_ratio = 0.0;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number:\n"
                       "       the maximum overlap allowed between a pair of detected blobs\n"
                       "       (as a fraction of the sum of their radii, estimated using r≈σ√3).");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-spheres-nonmax-overlap")



    else if ((vArgs[i] == "-spheres-nonmax-overlap-small") ||
             (vArgs[i] == "-max-volume-overlap-small") ||
             (vArgs[i] == "-max-overlap-small"))
    {
      try {
        if (i+1 >= vArgs.size())
          throw invalid_argument("");
        nonmax_max_volume_overlap_small = stof(vArgs[i+1]);
        nonmax_min_radial_separation_ratio = 0.0;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number:\n"
                       "       the maximum overlap allowed between a pair of detected blobs\n"
                       "       (as a fraction of the sum of their radii, estimated using r≈σ√3).");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-spheres-nonmax-overlap-small")



    else if ((vArgs[i] == "-spheres-nonmax-overlap-radial") ||
             (vArgs[i] == "-max-overlap-radial"))
    {
      try {
        if (i+1 >= vArgs.size())
          throw invalid_argument("");
        float max_overlap = stof(vArgs[i+1]);
        nonmax_min_radial_separation_ratio = (1.0 - max_overlap);
        //REMOVE THIS CRUFT:
        //nonmax_min_radial_separation_ratio = (1.0 - max_overlap) * sqrt(3.0);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number:\n"
                       "       the maximum overlap allowed between a pair of detected blobs\n"
                       "       (as a fraction of the sum of their radii, estimated using r≈σ√3).");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-spheres-nonmax-overlap-radial")

    

    // REMOVE THIS CRUFT
    //else if ((vArgs[i] == "-spheres-nonmax-separation-sigma") ||
    //         (vArgs[i] == "-blob-s-separation"))
    //{
    //  try {
    //    if (i+1 >= vArgs.size())
    //      throw invalid_argument("");
    //    nonmax_min_radial_separation_ratio = stof(vArgs[i+1]);
    //    //REMOVE THIS CRUFT
    //    //find_extrema_occlusion_ratio = stof(vArgs[i+1]);
    //  }
    //  catch (invalid_argument& exc) {
    //    throw InputErr("Error: The " + vArgs[i] + 
    //                   " argument must be followed by a number:\n"
    //                   "       the minimum allowed separation between a pair of detected blobs\n"
    //                   "       (as a fraction of the sum of their Gaussian widths (σ1+σ2)).");
    //  }
    //  num_arguments_deleted = 2;
    //} //if (vArgs[i] == "-blob-s-separation")




    else if ((vArgs[i] == "-blob-separation") ||
             (vArgs[i] == "-blob-r-separation") ||
             (vArgs[i] == "-blobr-separation") ||
             (vArgs[i] == "-radial-separation") ||
	     (vArgs[i] == "-spheres-nonmax-separation-radius"))
    {
      try {
        if (i+1 >= vArgs.size())
          throw invalid_argument("");
        nonmax_min_radial_separation_ratio = stof(vArgs[i+1]);
        //REMOVE THIS CRUFT
        //nonmax_min_radial_separation_ratio = stof(vArgs[i+1]) * sqrt(3.0);
        //find_extrema_occlusion_ratio = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number:\n"
                       "       the minimum allowed separation between a pair of detected blobs\n"
                       "       (as a fraction of the sum of their radii, estimated using r≈σ√3).");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-blob-separation")




    // REMOVE THIS CRUFT
    //else if ((vArgs[i] == "-spheres-nonmax-separation-diameter") ||
    //         (vArgs[i] == "-blob-d-separation") ||
    //         (vArgs[i] == "-blob-separation"))
    //{
    //  try {
    //    if (i+1 >= vArgs.size())
    //      throw invalid_argument("");
    //    nonmax_min_radial_separation_ratio = stof(vArgs[i+1]) * 2.0;
    //    //REMOVE THIS CRUFT
    //    //nonmax_min_radial_separation_ratio = stof(vArgs[i+1]) * 2.0 * sqrt(3.0);
    //    //find_extrema_occlusion_ratio = stof(vArgs[i+1]);
    //  }
    //  catch (invalid_argument& exc) {
    //    throw InputErr("Error: The " + vArgs[i] + 
    //                   " argument must be followed by a number:\n"
    //                   "       the minimum allowed separation between a pair of detected blobs\n"
    //                   "       (as a fraction of the sum of their diameters, estimated using d=2r≈2σ√3).");
    //  }
    //  num_arguments_deleted = 2;
    //} //if (vArgs[i] == "-blobd-separation")



    else if ((vArgs[i] == "-blob") ||
             (vArgs[i] == "-blob-sigma") ||
             (vArgs[i] == "-blob-s") ||
             (vArgs[i] == "-blob-radii") ||
             (vArgs[i] == "-blob-r") ||
             (vArgs[i] == "-blobr") ||
             (vArgs[i] == "-blob-diameters") ||
             (vArgs[i] == "-blob-d"))
    {
      try {
        if ((i+4 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-') ||
            (vArgs[i+4] == "") || (vArgs[i+4][0] == '-'))
          throw invalid_argument("");

        string blob_file_name_base = vArgs[i+1];
        //if (EndsWith(blob_file_name_base, ".txt"))
        //  blob_file_name_base =
        //    blob_file_name_base.substr(0,
        //                               blob_file_name_base.length()-4);

        blob_minima_file_name = blob_file_name_base + string(".minima.txt");
        blob_maxima_file_name = blob_file_name_base + string(".maxima.txt");

        float blob_width_min = stof(vArgs[i+2]);
        float blob_width_max = stof(vArgs[i+3]);
        if ((blob_width_min <= 0.0) ||
            (blob_width_max <= 0.0) ||
            (blob_width_min >= blob_width_max))
          throw invalid_argument("");

        // The old version of the code allows user to set "N" directly
        //int N = stof(vArgs[i+4]);
        //if (N < 3)
        //  throw InputErr("The 4th parameter to the \"" + vArgs[i] + "\" argument must be >= 3.");
        //float growth_ratio = pow(blob_width_max / blob_width_min, 1.0/(N-1));

        // In the new version of the code, the user specifies "growth_ratio"
        float growth_ratio = stof(vArgs[i+4]);
        if (growth_ratio <= 1.0)
          throw invalid_argument("");
        int N = 1+ceil(log(blob_width_max/blob_width_min) / log(growth_ratio));
        growth_ratio = pow(blob_width_max/blob_width_min, 1.0/N);
        blob_diameters.resize(N);

        // REMOVE THIS CRUFT:
        // Optional:
        // Make sure the difference between the width of the Gaussians used
        // to approximate the Laplacian-of-Gaussian (delta_sigma_over_sigma), 
        // is no larger than the difference between the successive widths
        // in the list of Gaussian blurs we will apply to our source image.
        // (Actually, I make sure it is no larger than 1/3 this difference.)
        //if (((growth_ratio-1.0)/3 < delta_sigma_over_sigma)
        //    &&
        //    (! delta_set_by_user))
        //  delta_sigma_over_sigma = (growth_ratio-1.0)/3;

        blob_width_multiplier = 1.0;
        if ((vArgs[i] == "-blob-sigma") || (vArgs[i] == "-blob-s"))
          // (for a solid uniform 3-D sphere)
          blob_width_multiplier = 2.0*sqrt(3.0);
          //REMOVE THIS CRUFT
          //blob_width_multiplier = 1.0/sqrt(3.0);
        if ((vArgs[i] == "-blob-radii") || (vArgs[i] == "-blob-r") || (vArgs[i] == "-blobr"))
          // (for a solid uniform 3-D sphere)
          blob_width_multiplier = 2.0;
          //REMOVE THIS CRUFT
          //blob_width_multiplier = 1.0/sqrt(3.0);
        else if ((vArgs[i] == "-blob-diameters") || (vArgs[i] == "-blob-d"))
          // (for a solid uniform 3-D sphere)
          blob_width_multiplier = 1.0;
          //REMOVE THIS CRUFT
          //blob_width_multiplier = 1.0/(2.0*sqrt(3.0));
        blob_diameters[0] = blob_width_min * blob_width_multiplier;
        for (int n = 1; n < N; n++) {
          blob_diameters[n] = blob_diameters[n-1] * growth_ratio;
        }
        m_exp = 2.0; // (<--not necessary, we ignore these parameters anyway)
        n_exp = 2.0; // (<--not necessary, we ignore these parameters anyway)
        filter_type = BLOB;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a file name, and 3 positive numbers:\n"
                       "      file_name.txt  blob_width_min  blob_width_max  N\n");
      }
      num_arguments_deleted = 5;
    } //if (vArgs[i] == "-blob")


    else if ((vArgs[i] == "-score-upper-bound") ||
             (vArgs[i] == "-blob-minima-threshold")) {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == ""))
          throw invalid_argument("");
        score_upper_bound = stof(vArgs[i+1]);
        score_lower_bound = std::numeric_limits<float>::infinity();
        score_bounds_are_ratios = false;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number number.\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-blob-minima-threshold")


    else if ((vArgs[i] == "-score-lower-bound") ||
             (vArgs[i] == "-blob-maxima-threshold")) {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == ""))
          throw invalid_argument("");
        score_lower_bound = stof(vArgs[i+1]);
        score_upper_bound = -std::numeric_limits<float>::infinity();
        score_bounds_are_ratios = false;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number number.\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-blob-maxima-threshold")
             

    else if ((vArgs[i] == "-score-lower-bound-ratio") ||
             (vArgs[i] == "-blob-minima-ratio")) {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == ""))
          throw invalid_argument("");
        score_upper_bound = stof(vArgs[i+1]);
        score_lower_bound = std::numeric_limits<float>::infinity();
        score_bounds_are_ratios = true;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number number.\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-blob-minima_ratio")



    else if ((vArgs[i] == "-score-upper-bound-ratio") ||
             (vArgs[i] == "-blob-maxima-ratio")) {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == ""))
          throw invalid_argument("");
        score_lower_bound = stof(vArgs[i+1]);
        score_upper_bound = -std::numeric_limits<float>::infinity();
        score_bounds_are_ratios = true;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number number.\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-blob-maxima-ratio")


    // REMOVE THIS CRUFT:
    //else if (vArgs[i] == "-score-upper-bound") {
    //  try {
    //    if ((i+1 >= vArgs.size()) ||
    //        (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
    //      throw invalid_argument("");
    //    score_upper_bound = stof(vArgs[i+1]);
    //  }
    //  catch (invalid_argument& exc) {
    //    throw InputErr("Error: The " + vArgs[i] + 
    //                   " argument must be followed by a number:\n"
    //                   "       scores above this (usually negative) number will be ignored\n");
    //  }
    //  num_arguments_deleted = 2;
    //}
    //
    //else if (vArgs[i] == "-score-lower-bound") {
    //  try {
    //    if ((i+1 >= vArgs.size()) ||
    //        (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
    //      throw invalid_argument("");
    //    score_lower_bound = stof(vArgs[i+1]);
    //  }
    //  catch (invalid_argument& exc) {
    //    throw InputErr("Error: The " + vArgs[i] + 
    //                   " argument must be followed by a number:\n"
    //                   "       Scores below this number will be ignored.\n");
    //  }
    //  num_arguments_deleted = 2;
    //}



    else if ((vArgs[i] == "-discard-blobs") or
	     (vArgs[i] == "-blob-nonmax") or
             (vArgs[i] == "-blobs-nonmax")) {
      try {
        if ((i+2 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+1] == vArgs[i+2]))
          throw invalid_argument("");
        in_coords_file_name = vArgs[i+1];
        out_coords_file_name = vArgs[i+2];
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by two different file names\n");
      }
      filter_type = SPHERE_NONMAX_SUPPRESSION;
      num_arguments_deleted = 3;
    } // if (vArgs[i] == "-spheres-nonmax")



    else if ((vArgs[i] == "-spheres-nonmax-radii-range") ||
             (vArgs[i] == "-sphere-nonmax-radii-range")) {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "-") || (vArgs[i+1][0] == '-') ||
            (i+2 >= vArgs.size()) || (vArgs[i+2] == "-") || (vArgs[i+2][0] == '-'))
          throw invalid_argument("");
        sphere_diameters_lower_bound = stof(vArgs[i+1]);
        sphere_diameters_upper_bound = stof(vArgs[i+2]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number number.\n");
      }
      num_arguments_deleted = 3;
    }  //if (vArgs[i] == "-spheres-nonmax-radii-range")



    else if ((vArgs[i] == "-spheres-nonmax-score-range") ||
             (vArgs[i] == "-sphere-nonmax-score-range")) {
      try {
        if (i+2 >= vArgs.size())
          throw invalid_argument("");
        score_lower_bound = stof(vArgs[i+1]);
        score_upper_bound = stof(vArgs[i+2]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by two numbers.\n");
      }
      num_arguments_deleted = 3;
    }  //if (vArgs[i] == "-spheres-nonmax-score-range")




    else if (vArgs[i] == "-dog-delta") {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        delta_sigma_over_sigma = stof(vArgs[i+1]);
        delta_set_by_user = true;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 1 positive number.\n");
      }
      num_arguments_deleted = 2;
    }



    else if ((vArgs[i] == "-exponents") ||
             (vArgs[i] == "-gdog-exponents"))
    {
      try {
        if ((i+2 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-'))
          throw invalid_argument("");
        m_exp = stof(vArgs[i+1]);
        n_exp = stof(vArgs[i+2]);
        #ifndef DISABLE_TEMPLATE_MATCHING
        template_background_exponent = n_exp;
        #endif
        exponents_set_by_user = true;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by two positive numbers.\n");
      }
      num_arguments_deleted = 3;
    } //if (vArgs[i] == "-gdog-exponents")


    else if ((vArgs[i] == "-gauss-exponent") ||
             (vArgs[i] == "-exponent"))
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        m_exp = stof(vArgs[i+1]);
        n_exp = m_exp;
        exponents_set_by_user = true;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a positive number.\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-gauss-exponent")


    //Commenting out:  The user cannot set filter_truncate_halfwidth directly anymore:
    //else if (vArgs[i] == "-window") {
    //  if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
    //    throw InputErr("Error: The " + vArgs[i] + 
    //                   " argument must be followed by either 1 or 3 positive integers.\n");
    //  filter_truncate_halfwidth[0] = stoi(vArgs[i+1]);
    //  num_arguments_deleted = 2;
    //  if ((i+4 >= vArgs.size()) && 
    //      (vArgs[i+2][0] != '-') && (vArgs[i+2] != "") &&
    //      (vArgs[i+3][0] != '-') && (vArgs[i+3] != "")) {
    //    filter_truncate_halfwidth[1] = stoi(vArgs[i+2]);
    //    filter_truncate_halfwidth[2] = stoi(vArgs[i+3]);
    //    num_arguments_deleted = 4;
    //  }
    //  else {
    //    filter_truncate_halfwidth[1] = filter_truncate_halfwidth[0];
    //    filter_truncate_halfwidth[2] = filter_truncate_halfwidth[0];
    //  }
    //} //if (vArgs[i] == "-window")


    else if (vArgs[i] == "-window-ratio")
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        filter_truncate_ratio = stof(vArgs[i+1]);
        filter_truncate_threshold=-1.0; //(disables)override any filter_truncate_threshold settings
        //filter_truncate_ratio_exp = exp(-pow(filter_truncate_ratio,
        //                                   n_exp));
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number.\n");
      }
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-window-ratio")


    else if ((vArgs[i] == "-cutoff") || (vArgs[i] == "-window-cutoff"))
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        filter_truncate_threshold = stof(vArgs[i+1]);
        filter_truncate_ratio = -1.0; //(disables) override any filter_truncate_ratio settings
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number.\n");
      }      
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-cutoff")


    else if (vArgs[i] == "-thresh") {
      try {
        if (i+1 >= vArgs.size())
          throw invalid_argument("");
        use_thresholds = true;
        use_dual_thresholds = false;
        out_threshold_01_a = stof(vArgs[i+1]);
        out_threshold_01_b = out_threshold_01_a;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 1 number.\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-thresh2") {
      try {
        if (i+2 >= vArgs.size())
          throw invalid_argument("");
        use_thresholds = true;
        use_dual_thresholds = false;
        out_threshold_01_a = stof(vArgs[i+1]);
        out_threshold_01_b = stof(vArgs[i+2]);
        out_thresh2_use_clipping = false;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 2 numbers.\n");
      }
      num_arguments_deleted = 3;
    }


    else if ((vArgs[i] == "-clip") ||
             (vArgs[i] == "-cl")) {
      try {
        if (i+2 >= vArgs.size())
          throw invalid_argument("");
        use_thresholds = true;
        use_dual_thresholds = false;
        out_threshold_01_a = stof(vArgs[i+1]);
        out_threshold_01_b = stof(vArgs[i+2]);
        out_thresh2_use_clipping = true;
        if (vArgs[i] == "-cl")
          out_thresh2_use_clipping_sigma = true;
        else
          out_thresh2_use_clipping_sigma = false;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 2 numbers.\n");
      }
      num_arguments_deleted = 3;
    }


    else if (vArgs[i] == "-thresh4") {
      try {
        if (i+4 >= vArgs.size())
          throw invalid_argument("");
        use_thresholds = true;
        use_dual_thresholds = true;
        out_threshold_01_a = stof(vArgs[i+1]);
        out_threshold_01_b = stof(vArgs[i+2]);
        out_threshold_10_a = stof(vArgs[i+3]);
        out_threshold_10_b = stof(vArgs[i+4]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 4 numbers.\n");
      }
      num_arguments_deleted = 5;
    }


    else if (vArgs[i] == "-bs") {
      #ifndef DISABLE_BOOTSTRAPPING
      try {
        if ((i+7 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-') ||
            (vArgs[i+4] == "") || (vArgs[i+4][0] == '-') ||
            (vArgs[i+5] == "") || (vArgs[i+5][0] == '-') ||
            (vArgs[i+6] == "") || (vArgs[i+6][0] == '-') ||
            (vArgs[i+7] == "") || (vArgs[i+7][0] == '-'))
          throw invalid_argument("");
        filter_type = BOOTSTRAP_DOGG;
        width_a[0] = stof(vArgs[i+1]);
        width_a[1] = width_a[0];
        width_a[2] = width_a[0];
        width_b[0] = stof(vArgs[i+2]);
        width_b[1] = width_b[0];
        width_b[2] = width_b[0];

        bs_ntests = stoi(vArgs[i+3]);
        float r = stof(vArgs[i+4]);
        f_bs_scramble_radius[0] = r;
        f_bs_scramble_radius[1] = r;
        f_bs_scramble_radius[2] = r;
        bs_threshold = stof(vArgs[i+5]);
        bs_threshold_sign = stof(vArgs[i+6]);
        bs_random_seed = stof(vArgs[i+7]);
        //mask_out = 1.0;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 7 numbers:\n"
                       "       ntests radius threshold sign seed\n");
      }
      num_arguments_deleted = 8;
      #else
      throw InputErr("Error: The " + vArgs[i] + 
                     " feature has been disabled in this version.\n"
                     "       To enable it, edit the \"settings.h\" file and comment out this line:\n"
                     "       \"#define DISABLE_BOOTSTRAPPING\"\n"
                     "       Then recompile.\n");
      #endif //#ifndef DISABLE_BOOTSTRAPPING
    }



    else if (vArgs[i] == "-template-gauss") {
      #ifndef DISABLE_TEMPLATE_MATCHING
      try {
        if ((i+2 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-'))
          throw invalid_argument("");
        filter_type = TEMPLATE_GAUSS;
        mask_out = 0.0;
        use_mask_out = true;
        width_a[0] = stof(vArgs[i+1]);
        width_a[1] = width_a[0];
        width_a[2] = width_a[0];
        float r = stof(vArgs[i+2]);
        template_background_radius[0] = r;
        template_background_radius[1] = r;
        template_background_radius[2] = r;
        //rescale01_out = false;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 2 numbers:\n"
                       "       template_radius background_radius\n");
      }
      num_arguments_deleted = 3;
      #else
      throw InputErr("Error: The " + vArgs[i] + 
                     " feature has been disabled in this version.\n"
                     "       To enable it, edit the \"settings.h\" file and comment out this line:\n"
                     "       \"#define DISABLE_TEMPLATE_MATCHING\"\n"
                     "       Then recompile.\n");
      #endif //#ifndef DISABLE_TEMPLATE_MATCHING
    }



    else if (vArgs[i] == "-template-gauss-aniso") {
      #ifndef DISABLE_TEMPLATE_MATCHING
      try {
        if ((i+6 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-') ||
            (vArgs[i+4] == "") || (vArgs[i+4][0] == '-') ||
            (vArgs[i+5] == "") || (vArgs[i+5][0] == '-') ||
            (vArgs[i+6] == "") || (vArgs[i+6][0] == '-'))
          throw invalid_argument("");
        filter_type = TEMPLATE_GAUSS;
        mask_out = 0.0;
        use_mask_out = true;
        width_a[0] = stof(vArgs[i+1]);
        width_a[1] = stof(vArgs[i+2]);
        width_a[2] = stof(vArgs[i+3]);
        template_background_radius[0] = stof(vArgs[i+4]);
        template_background_radius[1] = stof(vArgs[i+5]);
        template_background_radius[2] = stof(vArgs[i+6]);
        //rescale01_out = false;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 6 positive numbers.\n");
      }
      num_arguments_deleted = 7;
      #else
      throw InputErr("Error: The " + vArgs[i] + 
                     " feature has been disabled in this version.\n"
                     "       To enable it, edit the \"settings.h\" file and comment out this line:\n"
                     "       \"#define DISABLE_TEMPLATE_MATCHING\"\n"
                     "       Then recompile.\n");
      #endif //#ifndef DISABLE_TEMPLATE_MATCHING
    }


    else if ((vArgs[i] == "-fluct-aniso") ||
             (vArgs[i] == "-fluctuation-aniso") ||
             (vArgs[i] == "-fluctuations-aniso")) {
      #ifndef DISABLE_TEMPLATE_MATCHING
      try {
        if ((i+3 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-'))
          throw invalid_argument("");
        filter_type = LOCAL_FLUCTUATIONS;
        mask_out = 0.0;
        use_mask_out = true;
        template_background_radius[0] = stof(vArgs[i+1]);
        template_background_radius[1] = stof(vArgs[i+2]);
        template_background_radius[2] = stof(vArgs[i+3]);
        //rescale01_out = false;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 3 positive numbers.\n");
      }
      num_arguments_deleted = 4;
      #else
      throw InputErr("Error: The " + vArgs[i] + 
                     " feature has been disabled in this version.\n"
                     "       To enable it, edit the \"settings.h\" file and comment out this line:\n"
                     "       \"#define DISABLE_TEMPLATE_MATCHING\"\n"
                     "       Then recompile.\n");
      #endif //#ifndef DISABLE_TEMPLATE_MATCHING
    }


    else if ((vArgs[i] == "-fluct") ||
             (vArgs[i] == "-fluctuation") ||
             (vArgs[i] == "-fluctuations")) {
      #ifndef DISABLE_TEMPLATE_MATCHING
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        filter_type = LOCAL_FLUCTUATIONS;
        mask_out = 0.0;
        use_mask_out = true;
        template_background_radius[0] = stof(vArgs[i+1]);
        template_background_radius[1] = template_background_radius[0];
        template_background_radius[2] = template_background_radius[0];
        //rescale01_out = false;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a positive number:\n"
                       "       the radius of the region of nearby voxels over which we will perform\n"
                       "       the average before calculating fluctuations around that average.\n");
      }
      num_arguments_deleted = 2;
      #else
      throw InputErr("Error: The " + vArgs[i] + 
                     " feature has been disabled in this version.\n"
                     "       To enable it, edit the \"settings.h\" file and comment out this line:\n"
                     "       \"#define DISABLE_TEMPLATE_MATCHING\"\n"
                     "       Then recompile.\n");
      #endif //#ifndef DISABLE_TEMPLATE_MATCHING
    }


    else if (vArgs[i] == "-find-minima") {
      try {
        if (i+1 >= vArgs.size())
          throw invalid_argument("");
        find_minima = true;
        find_minima_file_name = vArgs[i+1];
      }
      catch (invalid_argument& exc) {
          throw InputErr("Error: The " + vArgs[i] + 
                         " argument must be followed by a number.\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-find-maxima") {
      try {
        if (i+1 >= vArgs.size())
          throw invalid_argument("");
        find_maxima = true;
        find_maxima_file_name = vArgs[i+1];
      }
      catch (invalid_argument& exc) {
          throw InputErr("Error: The " + vArgs[i] + 
                         " argument must be followed by a number.\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-distance") {
      if (i+1 >= vArgs.size())
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a file name.\n");
      filter_type = MIN_DISTANCE;
      in_coords_file_name = vArgs[i+1];
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-spheres") {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        filter_type = SPHERE_DECALS;
        sphere_decals_diameter = -1.0;
        in_coords_file_name = vArgs[i+1];
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a file name\n");
      }
      num_arguments_deleted = 2;
    }


    else if ((vArgs[i] == "-sphere-diameter") ||
             (vArgs[i] == "-sphere-diameters") ||
             (vArgs[i] == "-spheres-diameter") ||
             (vArgs[i] == "-spheres-diameters")) {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        sphere_decals_diameter = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number\n");
      }
      num_arguments_deleted = 2;
    }


    else if ((vArgs[i] == "-sphere-radius") ||
             (vArgs[i] == "-sphere-radii") ||
             (vArgs[i] == "-spheres-radius") ||
             (vArgs[i] == "-spheres-radii")) {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        sphere_decals_diameter = stof(vArgs[i+1]) * 2.0;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number\n");
      }
      num_arguments_deleted = 2;
    }


    else if ((vArgs[i] == "-spheres-scale") ||
             (vArgs[i] == "-sphere-scale")) {
      sphere_decals_foreground_use_score = false;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        sphere_decals_scale = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number:\n"
                       "       the ratio of the displyed sphere size to the diameter detected (usually 1).\n");
      }
      num_arguments_deleted = 2;
    }


    else if ((vArgs[i] == "-spheres-shell-ratio") ||
             (vArgs[i] == "-sphere-shell-ratio")) {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        sphere_decals_shell_thickness_is_ratio = true;
        sphere_decals_shell_thickness     = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a numbers:\n"
                       "       -the ratio of the shell thickness to the sphere diameter\n");
      }
      num_arguments_deleted = 2;
    }


    else if ((vArgs[i] == "-sphere-shell-thickness-min") ||
             (vArgs[i] == "-sphere-shell-thicknesses-min") ||
             (vArgs[i] == "-spheres-shell-thickness-min") ||
             (vArgs[i] == "-spheres-shell-thicknesses-min")) {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        sphere_decals_shell_thickness_min = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number\n");
      }
      num_arguments_deleted = 2;
    }


    else if ((vArgs[i] == "-sphere-shell-thickness") ||
             (vArgs[i] == "-sphere-shell-thicknesses") ||
             (vArgs[i] == "-spheres-shell-thickness") ||
             (vArgs[i] == "-spheres-shell-thicknesses")) {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        sphere_decals_shell_thickness = stof(vArgs[i+1]);
        sphere_decals_shell_thickness_is_ratio = false;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number\n");
      }
      num_arguments_deleted = 2;
    }



    else if ((vArgs[i] == "-spheres-score") ||
             (vArgs[i] == "-sphere-score")) {
      sphere_decals_foreground_use_score = true;
      num_arguments_deleted = 1;
    }

    else if ((vArgs[i] == "-spheres-background") ||
             (vArgs[i] == "-sphere-background")) {
      //sphere_decals_background_use_orig = false;
      sphere_decals_background_scale = 0.0;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        sphere_decals_background = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number:\n"
                       "       the voxel intensity value outside the sphere (normally 0).\n");
      }
      num_arguments_deleted = 2;
    }


    else if ((vArgs[i] == "-spheres-background-scale") ||
             (vArgs[i] == "-sphere-background-scale")) {
      //sphere_decals_background_use_orig = true;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        sphere_decals_background_scale = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number, usually between 0 and 1:\n"
                       "       how much to supress fluctuations in the original background image.\n");
      }
      num_arguments_deleted = 2;
    }


    else if ((vArgs[i] == "-spheres-foreground") ||
             (vArgs[i] == "-sphere-foreground")) {
      sphere_decals_foreground_use_score = false;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        sphere_decals_foreground = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number:\n"
                       "       the voxel intensity value outside the sphere (normally 0).\n");
      }
      num_arguments_deleted = 2;
    }


    else if ((vArgs[i] == "-spheres-normalize") ||
             (vArgs[i] == "-sphere-normalize")) {
      sphere_decals_foreground_norm = true;
      num_arguments_deleted = 1;
    }


    else if ((vArgs[i] == "-spheres01") || (vArgs[i] == "-spheres-01") ||
             (vArgs[i] == "-sphere01") || (vArgs[i] == "-sphere-01")) {
      sphere_decals_foreground_norm = false;
      num_arguments_deleted = 1;
    }



    else if (vArgs[i] == "-watershed")
    {
      filter_type = WATERSHED;
      try {
        if (i+1 < vArgs.size()) {
          if ((vArgs[i+1] == "min") || (vArgs[i+1] == "minima")) {
            watershed_use_minima = true;
            watershed_threshold = std::numeric_limits<float>::infinity();
          }
          else if ((vArgs[i+1] == "max") || (vArgs[i+1] == "maxima")) {
            watershed_use_minima = false;
            watershed_threshold = -std::numeric_limits<float>::infinity();
          }
          throw invalid_argument("");
        }
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be either \"minima\" or \"maxima\".\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-watershed-threshold")
    {
      filter_type = WATERSHED;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        watershed_threshold = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number\n");
      }
      
      num_arguments_deleted = 2;
    }



    else if ((vArgs[i] == "-planar") ||
             (vArgs[i] == "-plane") ||
             (vArgs[i] == "-membrane"))
    {
      filter_type = RIDGE_PLANAR;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        width_a[0] = stof(vArgs[i+1]) / sqrt(3.0);
        width_a[1] = width_a[0];
        width_a[2] = width_a[0];
        width_b[0] = 0.0;
        width_b[1] = 0.0;
        width_b[2] = 0.0;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number: the width\n"
                       "       (σ) of the Gaussian used beforehand for smoothing (in physical units).\n"
                       "       (This selects the scale or resolution of the filter.)\n");
      }
      num_arguments_deleted = 2;
    }


    else if ((vArgs[i] == "-planar-background") ||
             (vArgs[i] == "-plane-background") ||
             (vArgs[i] == "-membrane-background"))
    {
      filter_type = RIDGE_PLANAR;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        width_b[0] = stof(vArgs[i+1]);
        width_b[1] = width_b[0];
        width_b[2] = width_b[0];
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number: the width\n"
                       "       (σ) of the Gaussian used beforehand for smoothing (in physical units).\n"
                       "       (This selects the scale or resolution of the filter.)\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-planar-threshold") {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == ""))
          throw invalid_argument("");
        planar_hessian_score_threshold = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number, the threshold\n"
                       "       for deciding whether to compute the surface normal there.\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-planar-tv") {
      if (filter_type != RIDGE_PLANAR)
        throw InputErr("Error: the -planar-tv argument must appear in the argument list\n"
                       "       after the -planar argument.  It does not make sense to use\n"
                       "       -planar-tv refinement if you are not using the -planar filter.\n");
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == ""))
          throw invalid_argument("");
        planar_tv_sigma = stof(vArgs[i+1]);
        if (planar_tv_sigma < 1.0)
          throw InputErr("Error: the -planar-tv argument must be >= 1.0\n");
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number, the ratio\n"
                       "       of widths used in directional blurring (typically 3)\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-planar-tv-angle-exponent") {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == ""))
          throw invalid_argument("");
        planar_tv_exponent = stoi(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number, the threshold\n"
                       "       for deciding whether to compute the surface normal there.\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-planar-tv-threshold") {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == ""))
          throw invalid_argument("");
        planar_tv_score_threshold = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number, the threshold\n"
                       "       for deciding whether to compute the surface normal there.\n");
      }
      num_arguments_deleted = 2;
    }

    else if (vArgs[i] == "-planar-tv-iter") {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == ""))
          throw invalid_argument("");
        planar_tv_num_iters = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number, the threshold\n"
                       "       for deciding whether to compute the surface normal there.\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-planar-bnpts") {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        out_normals_fname = vArgs[i+1];
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a file name.\n");
      }
      num_arguments_deleted = 2;
    }

    else if (vArgs[i] == "-np") {
      #ifdef DISABLE_OPENMP
      throw InputErr("Error: The " + vArgs[i] + 
                     " argument is only available if program was compiled\n"
                     " with support for OpenMP (multiprocessor support).\n");
      #else
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        int num_threads = stoi(vArgs[i+1]);
        omp_set_num_threads(num_threads);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a positive integer, the\n"
                       "       number of threads (processors) requested.\n");
      }
      #endif //#ifdef DISABLE_OPENMP
      num_arguments_deleted = 2;
    }


    //else if (vArgs[i] == "-template-ggauss") {
    //  #ifndef DISABLE_TEMPLATE_MATCHING
    //  if ((i+2 >= vArgs.size()) ||
    //      (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
    //      (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
    //      (vArgs[i+3] == "") || (vArgs[i+3][0] == '-') ||
    //      (vArgs[i+4] == "") || (vArgs[i+4][0] == '-'))
    //    throw InputErr("Error: The " + vArgs[i] + 
    //                   " argument must be followed by 2 numbers:\n"
    //                   "       template_radius background_radius\n");
    //  filter_type = TEMPLATE_GGAUSS;
    //  mask_out = 0.0;
    //  use_mask_out = true;
    //  width_a[0] = stof(vArgs[i+1]);
    //  width_a[1] = width_a[0];
    //  width_a[2] = width_a[0];
    //  m_exp = stof(vArgs[i+2]);
    //  float r = stof(vArgs[i+3]);
    //  template_background_radius[0] = r;
    //  template_background_radius[1] = r;
    //  template_background_radius[2] = r;
    //  template_background_exponent = stof(vArgs[i+4]);
    //  rescale01_out = false;
    //  num_arguments_deleted = 5;
    //  #else
    //  throw InputErr("Error: The " + vArgs[i] + 
    //                 " feature has been disabled in this version.\n"
    //                 "       To enable it, edit the \"settings.h\" file and comment out this line:\n"
    //                 "       \"#define DISABLE_TEMPLATE_MATCHING\"\n"
    //                 "       Then recompile.\n");
    //  #endif //#ifndef DISABLE_TEMPLATE_MATCHING
    //}


    //else if (vArgs[i] == "-template-background-exponent") {
    //  if (i+2 >= vArgs.size())
    //    throw InputErr("Error: The " + vArgs[i] + 
    //                   " argument must be followed by 1 number.\n");
    //  template_background_exponent = stof(vArgs[i+1]);
    //}
    //
    //else if (vArgs[i] == "-template-compare-exponent") {
    //  if (i+2 >= vArgs.size())
    //    throw InputErr("Error: The " + vArgs[i] + 
    //                   " argument must be followed by 1 number.\n");
    //  template_compare_exponent = stof(vArgs[i+1]);
    //}

    // NOT USED YET:
    //else if (vArgs[i] == "-xwedge") {
    //  if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
    //    throw InputErr("Error: The " + vArgs[i] + 
    //                   " argument must be followed by 1 or 3 positive numbers.\n");
    //  missing_wedge_min[0] = stof(vArgs[i+1]);
    //  num_arguments_deleted = 2;
    //  if ((i+4 >= vArgs.size()) && 
    //      (vArgs[i+2][0] != '-') && (vArgs[i+2] != "")) {
    //    missing_wedge_max[0] = stof(vArgs[i+2]);
    //    num_arguments_deleted = 3;
    //  }
    //  else {
    //    missing_wedge_max[0] = missing_wedge_min[0];
    //    missing_wedge_min[0] = -missing_wedge_min[0];
    //  }
    //} //if (vArgs[i] == "-xwedge")
    //
    //else if ((vArgs[i] == "-ywedge") || (vArgs[i] == "-wedge")) {
    //  if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
    //    throw InputErr("Error: The " + vArgs[i] + 
    //                   " argument must be followed by 1 or 3 positive numbers.\n");
    //  missing_wedge_min[1] = stof(vArgs[i+1]);
    //  num_arguments_deleted = 2;
    //  if ((i+4 >= vArgs.size()) && 
    //      (vArgs[i+2][0] != '-') && (vArgs[i+2] != "")) {
    //    missing_wedge_max[1] = stof(vArgs[i+2]);
    //    num_arguments_deleted = 3;
    //  }
    //  else {
    //    missing_wedge_max[1] = missing_wedge_min[1];
    //    missing_wedge_min[1] = -missing_wedge_min[1];
    //  }
    //} //if (vArgs[i] == "-ywedge")



    //Delete all the arguments we have processed in this iteration (if any)
    if (num_arguments_deleted > 0)
    {
      vArgs.erase(vArgs.begin()+i, vArgs.begin()+i+num_arguments_deleted);
      --i; //rewind by one so that when we invoke ++i
      //at the end of this loop, i's value will not change.
      //This will point us to the next unread argument.
    }

  } // loop over arguments "for (int i=1; i < vArgs.size(); ++i)"


  if (vArgs.size() > 1) {
    stringstream err_msg;
    err_msg << "Error: Unrecognized argument: \""
            << vArgs[1] << "\"\n";
    throw InputErr(err_msg.str());
  }


  // ----------
  //if (filter_type == NONE)
  //  throw InputErr("Error: You must specify the kind of filter you want to use\n"
  //                 "       (Ie., the function you wish to convolve with your image, for example\n"
  //                 "        by using the \"-gauss\",\"-gdog\",\"-dog\",\"-file\",... arguments)\n");

  // ----------
  if ((in_file_name.size() == 0) &&
      (filter_type != SPHERE_NONMAX_SUPPRESSION))
  {
    throw InputErr("Error: You must specify the name of the tomogram you want to read\n"
                   "       using the \"-in SOURCE_FILE\" argument\n"
                   "       (These files usually end in \".mrc\" or \".rec\" .)\n");
  }
  if ((out_file_name.size() == 0) &&
      ((!
        ((filter_type == NONE) ||
         (filter_type == BLOB) ||
         (filter_type == SPHERE_NONMAX_SUPPRESSION)))
       ||
       use_thresholds ||
       invert_output))
  {
    throw InputErr("Error: You must specify the name of the tomogram you want to create\n"
                   "       using the \"-out DESTINATION_FILE\" argument\n"
                   "       (These files usually end in \".mrc\" or \".rec\" .)\n");
  }
  if ((out_file_name == in_file_name) && (out_file_name.size() != 0))
  {
    throw InputErr("Error: Input and Output files cannot be the same.\n");
  }


  // ----------
  if ((filter_type == GAUSS) || (filter_type == GGAUSS)) {
    // If the exponent equals 2, then the function we are convolving
    // with is an ordinary Gaussians.  Ordinary Gaussians are are separable
    // and can be convolved successively in 1-D in the x,y,z directions.
    // Hence we can save a great deal of time by doing this compared to
    // using a more general filter function which would force us to perform
    // the full 3-D convolution.
    if (exponents_set_by_user) {
      filter_type = GGAUSS;
      if (m_exp == 2.0) {
        filter_type = GAUSS; // <-- use fast seperable Gaussians
        width_a[0] /= sqrt(2.0);  // compensate for using
        width_a[1] /= sqrt(2.0);  //  exp(-0.5*(r/width_a)^2) instead of
        width_a[2] /= sqrt(2.0);  //  exp(-(r/width_a)^2)
      } //if (m_exp == 2.0)
    } //if (exponents_set_by_user)
  } //if ((filter_type == GAUSS) || (filter_type == GGAUSS))


  // ----------
  if ((filter_type == DOG) || (filter_type == DOGG)) {
    // If the exponents equal 2, then the functions we are convolving
    // with are ordinary Gaussians.  Ordinary Gaussians are are separable
    // and can be convolved successively in 1-D in the x,y,z directions.
    // Hence we can save a great deal of time by doing this compared to
    // using a more general filter function which would force us to perform
    // the full 3-D convolution.
    if (exponents_set_by_user) {
      filter_type = DOGG;
      if ((m_exp == 2.0) && (n_exp == 2.0)) {
        filter_type = DOG; // <-- use fast seperable Gaussians

        width_a[0] /= sqrt(2.0);  // compensate for using
        width_a[1] /= sqrt(2.0);  //  exp(-0.5*(r/width_a)^2) instead of
        width_a[2] /= sqrt(2.0);  //  exp(-(r/width_a)^2)
        width_b[0] /= sqrt(2.0);
        width_b[1] /= sqrt(2.0);
        width_b[2] /= sqrt(2.0);
      } //if ((m_exp == 2.0) && (n_exp == 2.0))
      else
        filter_type == DOGG;
    } //if (exponents_set_by_user) {
  } //if ((filter_type == DOG) || (filter_type == DOGG))


  if ((filter_type == TEMPLATE_GAUSS) || (filter_type == TEMPLATE_GGAUSS)) {
    // If the exponents equal 2, then the functions we are convolving
    // with are ordinary Gaussians.  Ordinary Gaussians are are separable
    // and can be convolved successively in 1-D in the x,y,z directions.
    // Hence we can save a great deal of time by doing this compared to
    // using a more general filter function which would force us to perform
    // the full 3-D convolution.
    if (exponents_set_by_user) {
      filter_type = TEMPLATE_GGAUSS;
      if ((m_exp == 2.0) && (n_exp == 2.0)) {
        filter_type = TEMPLATE_GAUSS; // <-- use fast Gaussians instead of DOGG which is slow

        width_a[0] /= sqrt(2.0);  // compensate for using
        width_a[1] /= sqrt(2.0);  //  exp(-0.5*(r/width_a)^2) instead of
        width_a[2] /= sqrt(2.0);  //  exp(-(r/width_a)^2)
        width_b[0] /= sqrt(2.0);
        width_b[1] /= sqrt(2.0);
        width_b[2] /= sqrt(2.0);
        template_background_radius[0] /= sqrt(2.0);// compensate for using:
        template_background_radius[1] /= sqrt(2.0);//  exp(-0.5*(r/w)^2)
        template_background_radius[2] /= sqrt(2.0);// instead of exp(-(r/w)^2)
      } //if ((m_exp == 2.0) && (n_exp == 2.0))
      else
        filter_type == TEMPLATE_GGAUSS;
    } //if (exponents_set_by_user) {
  } //if ((filter_type == DOG) || (filter_type == DOGG))


  if (filter_type == BLOB) {
    if ((blob_minima_file_name != "") &&
        (score_lower_bound == std::numeric_limits<float>::infinity()))
    {
      // If the user is not interested in local maxima, then we only
      // need to create one file (for the minima, as opposed to a separate
      // file for both minima and maxima).  In that case, remove the
      // ".minima.txt" which was earlier automatically added to the end
      // of that file's name.
      string blob_file_name_base = blob_minima_file_name;
      if (EndsWith(blob_file_name_base, ".minima.txt"))
        blob_file_name_base =
          blob_file_name_base.substr(0,
                                     blob_file_name_base.length()-11);
      blob_minima_file_name = blob_file_name_base;
    }
    if ((blob_maxima_file_name != "") &&
        (score_upper_bound == -std::numeric_limits<float>::infinity()))
    {
      // If the user is not interested in local minima, then we only
      // need to create one file (for the maxima, as opposed to a separate
      // file for both minima and maxima).  In that case, remove the
      // ".maxima.txt" which was earlier automatically added to the end
      // of that file's name.
      string blob_file_name_base = blob_maxima_file_name;
      if (EndsWith(blob_file_name_base, ".maxima.txt"))
        blob_file_name_base =
          blob_file_name_base.substr(0,
                                     blob_file_name_base.length()-11);
      blob_maxima_file_name = blob_file_name_base;
    }
  } //if (filter_type == BLOB)


  if (filter_type == RIDGE_PLANAR) {
    float sigma = width_a[0];
    // The parameter entered by the user is a ratio, not an absolute number.
    // Now that we know what sigma is, multiply planar_tv_sigma by sigma.
    planar_tv_sigma *= sigma;
  }

} // Settings::ParseArgs()



bool EndsWith(string const &fullString, string const &ending) {
  if (fullString.length() >= ending.length()) {
    return (0 ==
            fullString.compare(fullString.length() - ending.length(),
                               ending.length(),
                               ending));
  }
  else return false;
}
