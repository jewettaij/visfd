#include <sstream>
#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>
#include <set>
#include <array>
using namespace std;

#ifndef DISABLE_OPENMP
#include <omp.h>       // (OpenMP-specific)
#endif
#include "err.hpp"  // defines InputErr exception class
#include "file_io.hpp"
#include "settings.hpp"



Settings::Settings() {
  // Default settings
  voxel_width = 0.0;  //How many Angstroms per voxel? (if 0 then read from file)
  voxel_width_divide_by_10 = false; //Use nm instead of Angstroms?
  filter_type = NONE;

  in_file_name = "";
  in_set_image_size[0] = 0;
  in_set_image_size[1] = 0;
  in_set_image_size[2] = 0;
  invert_output = false;
  out_file_name = "";
  out_file_overwrite = false;
  mask_file_name = "";
  normalize_near_boundaries = true;
  mask_select = 1;
  use_mask_select = false;
  specify_masked_brightness = true;
  masked_voxel_brightness = 0.0;
  undefined_voxel_brightness = -1.0;
  undefined_voxels_are_max = true;
  vMaskRegions.resize(0);
  is_mask_crds_in_voxels = true; //are mask crds in units of voxels or physical distance?

  resize_with_binning = 0; //width of bin to use to reduce image size?

  save_intermediate_fname_base = "";
  load_intermediate_fname_base = "";

  median_radius = 0.0;
  morphology_r = 0.0;
  morphology_rmax = 0.0;
  morphology_bmax = 0.0;

  // The code is a little bit confusing because I tried to cram as much
  // functionality into the smallest number of parameters possible:
  // The "difference of gaussians" (dog) filter used here requires
  //  8 parameters: a_x, a_y, a_z, b_x, b_y, b_z, m, n
  // For details, see the comment at the beginning of "settings.h".
  //
  // Default values of a_x, a_y, a_z, b_x, b_y, b_z, m, n are:

  width_a[0] = -1.0;  // a_x = gaussian width in x direction(physical units)
  width_a[1] = -1.0;  // a_y = gaussian width in y direction
  width_a[2] = -1.0;  // a_z = gaussian width in z direction
  width_b[0] = -1.0;  // b_x = 2nd (outer) gaussian width in x direction
  width_b[1] = -1.0;  // b_y = gaussian width in y direction
  width_b[2] = -1.0;  // b_z = gaussian width in z direction(<0 disables)
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

  log_width[0] = 0.0;
  log_width[1] = 0.0;
  log_width[2] = 0.0;
  delta_sigma_over_sigma = 0.02;

  find_minima = false;
  find_maxima = false;
  find_minima_file_name = "";
  find_maxima_file_name = "";
  neighbor_connectivity = 3;
  extrema_on_boundary = true;
  //REMOVE THIS CRUFT
  //find_extrema_occlusion_ratio = 1.0;
  in_coords_file_name = "";
  out_coords_file_name = "";
  sphere_decals_diameter = -1.0;
  sphere_decals_foreground = 1.0;
  sphere_decals_background = 0.0;
  sphere_decals_background_scale = 1.0;
  sphere_decals_foreground_use_score = true;
  //sphere_decals_background_use_orig = true;
  sphere_decals_background_norm = false;
  sphere_decals_foreground_norm = false;
  sphere_decals_scale = 1.0;
  //sphere_decals_shell_thickness = -1.0;
  sphere_decals_shell_thickness = 1.0;
  sphere_decals_shell_thickness_is_ratio = true;
  sphere_decals_shell_thickness_min = 1.0;
  score_lower_bound = -std::numeric_limits<float>::infinity();
  score_upper_bound = std::numeric_limits<float>::infinity();
  score_bounds_are_ratios = false;
  sphere_diameters_lower_bound = -std::numeric_limits<float>::infinity();
  sphere_diameters_upper_bound = std::numeric_limits<float>::infinity();
  training_pos_crds.clear();
  training_neg_crds.clear();
  is_training_pos_in_voxels = false;
  is_training_neg_in_voxels = false;
  blob_width_multiplier = 1.0;
  nonmax_min_radial_separation_ratio = 0.0;
  nonmax_max_volume_overlap_small = std::numeric_limits<float>::infinity();
  nonmax_max_volume_overlap_large = std::numeric_limits<float>::infinity();
  auto_thresh_score = false;

  // ---- parameters used by the curve and surface detectors ----
  out_normals_fname = "";
  ridges_are_maxima = false;
  //max_distance_to_feature = std::numeric_limits<float>::infinity();
  //max_distance_to_feature = std::sqrt(3.0)/2;
  max_distance_to_feature = 1.3; // chosen after some experimentation. seems ok
  surface_normal_curve_ds = 0.2;
  surface_find_ridge = true;
  hessian_score_threshold = 0.05; //discard voxels which are not the best 5%
  hessian_score_threshold_is_a_fraction = true;
  tv_score_threshold = 0.0;
  tv_sigma = 0.0;
  tv_exponent = 4;
  tv_truncate_ratio = sqrt(2.0);

  // ---- parameters for watershed segmentation (and clustering) ----
  clusters_begin_at_maxima = false;
  watershed_boundary_label = 0.0;
  watershed_show_boundaries = true;
  watershed_threshold = std::numeric_limits<float>::infinity();
  watershed_markers_filename = "";
  cluster_connected_voxels = false;
  connect_threshold_saliency = std::numeric_limits<float>::infinity();
  must_link_constraints.clear();
  must_link_constraint_directions.clear();
  is_must_link_constraints_in_voxels = false;
  select_cluster = 0;

  // ---- parameters for directional watershed segmentation ----
  //connect_threshold_vector_saliency = -std::numeric_limits<float>::infinity();
  //connect_threshold_vector_neighbor = -std::numeric_limits<float>::infinity();
  //connect_threshold_tensor_saliency = -std::numeric_limits<float>::infinity();
  //connect_threshold_tensor_neighbor = -std::numeric_limits<float>::infinity();
  connect_threshold_vector_saliency = std::cos(M_PI*15/180.0);//15 degree change
  connect_threshold_vector_neighbor = std::cos(M_PI*15/180.0);//allowed between
  connect_threshold_tensor_saliency = std::cos(M_PI*15/180.0);//neighboring
  connect_threshold_tensor_neighbor = std::cos(M_PI*15/180.0);//voxels

  // --- parameters for intensity maps and thresholding ---
  use_intensity_map = false;
  use_dual_thresholds = false;
  in_threshold_01_a = 0.0;
  in_threshold_01_b = 0.0;
  in_threshold_10_a = 0.0;
  in_threshold_10_b = 0.0;
  out_thresh2_use_clipping = false;
  out_thresh2_use_clipping_sigma = false;
  out_thresh_a_value = 0.0;
  out_thresh_b_value = 1.0;
  rescale_min_max_in = false;
  rescale_min_max_out = false;
  in_rescale_min = 0.0;
  in_rescale_max = 1.0;
  out_rescale_min = 0.0;
  out_rescale_max = 1.0;
  use_rescale_multiply = false;
  out_rescale_multiply = 1.0;
  out_rescale_offset = 0.0;
  use_gauss_thresholds = false;
  out_thresh_gauss_x0 = 0.0;
  out_thresh_gauss_sigma = 1.0;
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
  template_background_radius[0] = -1.0;         //impossible value
  template_background_radius[1] = -1.0;         //impossible value
  template_background_radius[2] = -1.0;         //impossible value
  template_background_exponent = 2.0;         //default value
  //template_compare_exponent = 2.0;            //default value
  #endif //#ifndef DISABLE_TEMPLATE_MATCHING

  #ifndef DISABLE_INTENSITY_PROFILES
  blob_profiles_file_name_base = "";
  blob_profiles_center_criteria = BlobCenterCriteria::CENTER;
  #endif

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
  bool user_set_thickness_manually = false;
  bool user_set_background_scale_manually = false;
  bool user_set_exponents_manually = false;
  bool user_set_delta_manually = false;
  bool user_set_watershed_threshold_manually = false;

  // training data for a single images
  string training_pos_fname = "";
  string training_neg_fname = "";
  // training data in multiple independent files (for multiple images)
  vector<string> multi_training_neg_fnames;
  vector<string> multi_training_pos_fnames;

  // PASS 1: Parse all the arguments that alter the argument list.
  for (int i=1; i < vArgs.size(); ++i)
  {
    if (vArgs[i] == "-save-progress")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        save_intermediate_fname_base = vArgs[i+1];
        int nf = save_intermediate_fname_base.size();
        string fname_suffix = "";
        if (nf > 4)
          fname_suffix = save_intermediate_fname_base.substr(nf-4, nf);
        string fname_base = save_intermediate_fname_base;
        if ((fname_suffix == ".mrc") || (fname_suffix == ".rec") ||
            (fname_suffix == ".MRC") || (fname_suffix == ".REC"))
          fname_base = save_intermediate_fname_base.substr(0, nf-4);
        save_intermediate_fname_base = fname_base;
        string save_intermediate_fname_info = fname_base + "_info.txt";
        // Now save all the current command line arguments into the INFO file.
        fstream f;
        f.open(save_intermediate_fname_info, ios::out);
        for (int j = 1; j < vArgs.size(); j++) {
          if ((j == i) || (j == i+1))
            continue;  // omit the "-save-progress filename" arguments
          else if ((vArgs[j] == "-cl") ||
                   ((j>1) && (vArgs[j-1] == "-cl")) ||
                   ((j>2) && (vArgs[j-2] == "-cl")))
            continue;  // omit the "-cl a b" arguments
          else if ((vArgs[j] == "-clip") ||
                   ((j>1) && (vArgs[j-1] == "-clip")) ||
                   ((j>2) && (vArgs[j-2] == "-clip")))
            continue;  // omit the "-clip a b" arguments
          f << vArgs[j] << "\n";
        }
        f.close();
        vArgs.erase(vArgs.begin()+i,
                    vArgs.begin()+i+2);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a file name\n"
                       "       (The suffix at the end of the file (eg \".mrc\" or \".rec\") is optional.)\n");
      }
    } //if (vArgs[i] == "-save-progress")

    else if (vArgs[i] == "-load-progress")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        load_intermediate_fname_base = vArgs[i+1];
        int nf = load_intermediate_fname_base.size();
        string fname_suffix = "";
        if (nf > 4)
          fname_suffix = load_intermediate_fname_base.substr(nf-4, nf);
        string fname_base = load_intermediate_fname_base;
        if ((fname_suffix == ".mrc") || (fname_suffix == ".rec") ||
            (fname_suffix == ".MRC") || (fname_suffix == ".REC"))
          fname_base = load_intermediate_fname_base.substr(0, nf-4);
        load_intermediate_fname_base = fname_base;
        string load_intermediate_fname_info = fname_base + "_info.txt";
        // Now load the arguments that were in the INFO file.
        fstream f;
        f.open(load_intermediate_fname_info, ios::in);
        vector<string> new_args;
        string line;
        while (getline(f, line, '\n'))
          if (line.size() > 0)
            new_args.push_back(line);
        f.close();
        vArgs.erase(vArgs.begin() + i,
                    vArgs.begin() + i + 2);
        vArgs.insert(vArgs.begin() + 1,
                     new_args.begin(),
                     new_args.end());
      } // try
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a file name\n");
      }
    } // else if (vArgs[i] == "-load-progress")
  } // for (int i=1; i < vArgs.size(); ++i)  (PASS 1)




  // PASS 2: Parse all of the ordinary arguments.
  //         (that don't alter with the argument list)
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



    else if (vArgs[i] == "-image-size")
    {
      try {
        if ((i+3 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-'))
          throw invalid_argument("");
        in_set_image_size[0] = stoi(vArgs[i+1]);
        in_set_image_size[1] = stoi(vArgs[i+2]);
        in_set_image_size[2] = stoi(vArgs[i+3]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by 6 integers.\n");
      }
      num_arguments_deleted = 4;
    } // if (vArgs[i] == "-image-size")



    else if ((vArgs[i] == "-out") || (vArgs[i] == "-o"))
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a file name.\n");
      out_file_name = vArgs[i+1];
      num_arguments_deleted = 2;
    } // if ((vArgs[i] == "-out") || (vArgs[i] == "-o"))
   


    else if ((vArgs[i] == "-outf") || (vArgs[i] == "-out-force"))
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a file name.\n");
      out_file_name = vArgs[i+1];
      out_file_overwrite = true;
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-outf")



    else if (vArgs[i] == "-mask")
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a file name.\n");
      mask_file_name = vArgs[i+1];
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-mask")



    else if (vArgs[i] == "-normalize-filters")
    {
      InputErr arg_err("Error: The " + vArgs[i] +
                      " argument must be followed by \"yes\" or \"no\".\n");
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw arg_err;
      if (vArgs[i+1] == "yes")
        normalize_near_boundaries = true;
      if (vArgs[i+1] == "no")
        normalize_near_boundaries = false;
      else
        throw arg_err;
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-mask")

    

    else if (vArgs[i] == "-mask-select")
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == ""))
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



    else if ((vArgs[i] == "-mask-rect") ||
             (vArgs[i] == "-mask-rectangle"))
    {
      try {
        if ((i+6 >= vArgs.size()) ||
            (vArgs[i+1] == "") ||
            (vArgs[i+2] == "") ||
            (vArgs[i+3] == "") ||
            (vArgs[i+4] == "") ||
            (vArgs[i+5] == "") ||
            (vArgs[i+6] == ""))
          throw invalid_argument("");
        SimpleRegion<float> region;
        region.type = SimpleRegion<float>::RECT;
        region.value = 1;
        region.data.rect.xmin = stof(vArgs[i+1]);
        region.data.rect.xmax = stof(vArgs[i+2]);
        region.data.rect.ymin = stof(vArgs[i+3]);
        region.data.rect.ymax = stof(vArgs[i+4]);
        region.data.rect.zmin = stof(vArgs[i+5]);
        region.data.rect.zmax = stof(vArgs[i+6]);
        vMaskRegions.push_back(region);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by 6 numbers.\n");
      }
      num_arguments_deleted = 7;
    } // if (vArgs[i] == "-mask-rect")



    else if ((vArgs[i] == "-mask-rect-subtract") ||
             (vArgs[i] == "-mask-rectangle-subtract"))
    {
      try {
        if ((i+6 >= vArgs.size()) ||
            (vArgs[i+1] == "") ||
            (vArgs[i+2] == "") ||
            (vArgs[i+3] == "") ||
            (vArgs[i+4] == "") ||
            (vArgs[i+5] == "") ||
            (vArgs[i+6] == ""))
          throw invalid_argument("");
        SimpleRegion<float> region;
        region.type = SimpleRegion<float>::RECT;
        region.value = -1;
        region.data.rect.xmin = stof(vArgs[i+1]);
        region.data.rect.xmax = stof(vArgs[i+2]);
        region.data.rect.ymin = stof(vArgs[i+3]);
        region.data.rect.ymax = stof(vArgs[i+4]);
        region.data.rect.zmin = stof(vArgs[i+5]);
        region.data.rect.zmax = stof(vArgs[i+6]);
        vMaskRegions.push_back(region);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by 6 numbers.\n");
      }
      num_arguments_deleted = 7;
    } // if (vArgs[i] == "-mask-rect-subtract")



    else if (vArgs[i] == "-mask-sphere")
    {
      try {
        if ((i+4 >= vArgs.size()) ||
            (vArgs[i+1] == "") ||
            (vArgs[i+2] == "") ||
            (vArgs[i+3] == "") ||
            (vArgs[i+4] == ""))
          throw invalid_argument("");
        SimpleRegion<float> region;
        region.type = SimpleRegion<float>::SPHERE;
        region.value = 1;
        region.data.sphere.x0 = stof(vArgs[i+1]);
        region.data.sphere.y0 = stof(vArgs[i+2]);
        region.data.sphere.z0 = stof(vArgs[i+3]);
        region.data.sphere.r = stof(vArgs[i+4]);
        vMaskRegions.push_back(region);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by 4 numbers.\n");
      }
      num_arguments_deleted = 5;
    } // if (vArgs[i] == "-mask-sphere")



    else if (vArgs[i] == "-mask-sphere-subtract")
    {
      try {
        if ((i+4 >= vArgs.size()) ||
            (vArgs[i+1] == "") ||
            (vArgs[i+2] == "") ||
            (vArgs[i+3] == "") ||
            (vArgs[i+4] == ""))
          throw invalid_argument("");
        SimpleRegion<float> region;
        region.type = SimpleRegion<float>::SPHERE;
        region.value = -1;
        region.data.sphere.x0 = stof(vArgs[i+1]);
        region.data.sphere.y0 = stof(vArgs[i+2]);
        region.data.sphere.z0 = stof(vArgs[i+3]);
        region.data.sphere.r = stof(vArgs[i+4]);
        vMaskRegions.push_back(region);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by 4 numbers.\n");
      }
      num_arguments_deleted = 5;
    } // if (vArgs[i] == "-mask-sphere-subtract")



    else if ((vArgs[i] == "-mask-crds-units") ||
             (vArgs[i] == "-mask-coords-units") ||
             (vArgs[i] == "-mask-coordinates-units"))
    {
      try {
        if ((i+1 >= vArgs.size()) && (vArgs[i+1] == "voxels"))
          is_mask_crds_in_voxels = false;
        if ((i+1 >= vArgs.size()) &&
            ((vArgs[i+1] == "distance") ||
             (vArgs[i+1] == "physical") ||
             (vArgs[i+1] == "Ansgroms") ||
             (vArgs[i+1] == "angstroms") ||
             (vArgs[i+1] == "nm") ||
             (vArgs[i+1] == "nanometers")))
          is_mask_crds_in_voxels = true;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by either \"voxels\" or \"distance\".\n");
      }
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-mask-rect-units")



    else if (vArgs[i] == "-mask-out")
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == ""))
          throw invalid_argument("");
        specify_masked_brightness = true;
        masked_voxel_brightness = stof(vArgs[i+1]);
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
    } // if (vArgs[i] == "-w")



    else if ((vArgs[i] == "-a2nm") || (vArgs[i] == "-ang-to-nm"))
    {
      voxel_width_divide_by_10 = true;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-a2nm")



    else if (vArgs[i] == "-bin")
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        resize_with_binning = stoi(vArgs[i+1]);
        if (resize_with_binning < 1)
          throw invalid_argument("");
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a positive integer.\n");
      }
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-bin")



    else if ((vArgs[i] == "-dilation") ||
             (vArgs[i] == "-dilate"))
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        morphology_r = stof(vArgs[i+1]);
        filter_type = DILATION;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a nonnegative number\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-dilation")



    else if ((vArgs[i] == "-erosion") ||
             (vArgs[i] == "-erode"))
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        morphology_r = stof(vArgs[i+1]);
        filter_type = EROSION;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a nonnegative number\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-erosion")



    else if ((vArgs[i] == "-dilation-binary-soft") ||
             (vArgs[i] == "-dilate-binary-soft"))
    {
      try {
        if ((i+2 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-'))
          throw invalid_argument("");
        morphology_r = stof(vArgs[i+1]);
        morphology_rmax = stof(vArgs[i+2]);
        morphology_bmax = stof(vArgs[i+3]);
        filter_type = DILATION;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by nonnegative numbers\n");
      }
      num_arguments_deleted = 4;
    } //(vArgs[i] == "-dilation-binary-soft")



    else if ((vArgs[i] == "-erosion-binary-soft") ||
             (vArgs[i] == "-erode-binary-soft"))
    {
      try {
        if ((i+3 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-'))
          throw invalid_argument("");
        morphology_r = stof(vArgs[i+1]);
        morphology_rmax = stof(vArgs[i+2]);
        morphology_bmax = stof(vArgs[i+3]);
        filter_type = EROSION;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by nonnegative numbers\n");
      }
      num_arguments_deleted = 4;
    } //(vArgs[i] == "-erosion-binary-soft")



    else if ((vArgs[i] == "-dilation-gauss") ||
             (vArgs[i] == "-dilate-gauss") ||
             (vArgs[i] == "-erosion-gauss") ||
             (vArgs[i] == "-erode-gauss"))
    {
      float blur_distance;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        blur_distance = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a nonnegative number\n");
      }
      filter_type = GAUSS;
      // First, specify the "thickness" of the opening or closing operation.
      // This is the same as the "blur_distance" used in Gaussian blurring.
      width_a[0] = blur_distance;
      width_a[1] = width_a[0];
      width_a[2] = width_a[0];
      // Then append a post-processing threshold filter:
      use_intensity_map = true;
      if ((vArgs[i] == "-dilation-gauss") ||
          (vArgs[i] == "-dilate-gauss"))
        in_threshold_01_a = 0.1572992070502851; // ≈ 1-erf(1)
      else if ((vArgs[i] == "-erosion-gauss") ||
               (vArgs[i] == "-erode-gauss"))
        in_threshold_01_a = 0.8427007929497149;  // ≈ erf(1)
      in_threshold_01_b = in_threshold_01_a;

      num_arguments_deleted = 2;
    } // if ((vArgs[i] == "-erode-gauss") || (vArgs[i] == "-dilate-gauss"))



    else if ((vArgs[i] == "-opening") ||
             (vArgs[i] == "-open"))
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        morphology_r = stof(vArgs[i+1]);
        filter_type = OPENING;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a nonnegative number\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-opening")



    else if ((vArgs[i] == "-closing") ||
             (vArgs[i] == "-close"))
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        morphology_r = stof(vArgs[i+1]);
        filter_type = CLOSING;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a nonnegative number\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-closing")



    else if (vArgs[i] == "-top-hat-white")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        morphology_r = stof(vArgs[i+1]);
        filter_type = TOP_HAT_WHITE;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a nonnegative number\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-top-hat-white")



    else if (vArgs[i] == "-top-hat-black")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        morphology_r = stof(vArgs[i+1]);
        filter_type = TOP_HAT_BLACK;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a nonnegative number\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-top-hat-black")



    else if (vArgs[i] == "-truncate")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
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
    } // if (vArgs[i] == "-truncate")



    else if (vArgs[i] == "-truncate-threshold")
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
    } // if (vArgs[i] == "-truncate-thresold")



    else if (vArgs[i] == "-rescale")
    {
      try {
        if ((i+2 >= vArgs.size()) ||
            (vArgs[i+1] == "") ||
            (vArgs[i+2] == ""))
          throw invalid_argument("");
        use_intensity_map = true;
        use_rescale_multiply = true;
        out_rescale_multiply = stof(vArgs[i+1]);
        out_rescale_offset = stof(vArgs[i+2]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by 2 numbers:\n"
                       " outA  outB\n"
                       "  (the desired minimum and maximum voxel intensity values for the final image)\n");
      }
      num_arguments_deleted = 3;
    } // if (vArgs[i] == "-rescale")



    else if (vArgs[i] == "-fill")
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == ""))
          throw invalid_argument("");
        use_intensity_map = true;
        use_rescale_multiply = true;
        out_rescale_multiply = 0.0;
        out_rescale_offset = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number.");
      }
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-rescale")



    else if ((vArgs[i] == "-thresh-range") ||
             (vArgs[i] == "-thresh-range-out"))
    {
      try {
        if ((i+2 >= vArgs.size()) ||
            (vArgs[i+1] == "") ||
            (vArgs[i+2] == ""))
          throw invalid_argument("");
        out_thresh_a_value = stof(vArgs[i+1]);
        out_thresh_b_value = stof(vArgs[i+2]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by 2 numbers:\n"
                       " outA  outB\n"
                       "  (the desired minimum and maximum voxel intensity values for the final image)\n");
      }
      num_arguments_deleted = 3;
    } // if (vArgs[i] == "-outab")



    else if (vArgs[i] == "-rescale-min-max")
    {
      try {
        if ((i+2 >= vArgs.size()) ||
            (vArgs[i+1] == "") ||
            (vArgs[i+2] == ""))
          throw invalid_argument("");
        rescale_min_max_out = true;
        out_rescale_max = stof(vArgs[i+1]);
        out_rescale_min = stof(vArgs[i+2]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by 2 numbers:\n"
                       " outA  outB\n"
                       "  (the desired minimum and maximum voxel intensity values for the final image)\n");
      }
      num_arguments_deleted = 3;
    } // if (vArgs[i] == "-rescale-min-max")



    else if ((vArgs[i] == "-no-rescale") || (vArgs[i] == "-norescale"))
    {
      rescale_min_max_out = false;
      in_threshold_01_a = 1.0;
      in_threshold_01_b = 1.0;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-norescale")



    else if ((vArgs[i] == "-invert") ||
             (vArgs[i] == "-inv"))
    {
      invert_output = true;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-invert")



    else if ((vArgs[i] == "-thresh") ||
             (vArgs[i] == "-thresh-out")) {
      try {
        if (i+1 >= vArgs.size())
          throw invalid_argument("");
        use_intensity_map = true;
        use_dual_thresholds = false;
        in_threshold_01_a = stof(vArgs[i+1]);
        in_threshold_01_b = in_threshold_01_a;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by 1 number.\n");
      }
      num_arguments_deleted = 2;
    }



    else if ((vArgs[i] == "-thresh2") ||
             (vArgs[i] == "-thresh2-out")) {
      try {
        if (i+2 >= vArgs.size())
          throw invalid_argument("");
        use_intensity_map = true;
        use_dual_thresholds = false;
        in_threshold_01_a = stof(vArgs[i+1]);
        in_threshold_01_b = stof(vArgs[i+2]);
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
        use_intensity_map = true;
        use_dual_thresholds = false;
        in_threshold_01_a = stof(vArgs[i+1]);
        in_threshold_01_b = stof(vArgs[i+2]);
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



    else if ((vArgs[i] == "-thresh4") ||
             (vArgs[i] == "-thresh4-out")) {
      try {
        if (i+4 >= vArgs.size())
          throw invalid_argument("");
        use_intensity_map = true;
        use_dual_thresholds = true;
        in_threshold_01_a = stof(vArgs[i+1]);
        in_threshold_01_b = stof(vArgs[i+2]);
        in_threshold_10_a = stof(vArgs[i+3]);
        in_threshold_10_b = stof(vArgs[i+4]);
        bool all_increasing = ((in_threshold_01_a <= in_threshold_01_b) &&
                               (in_threshold_01_b <= in_threshold_10_a) &&
                               (in_threshold_10_a <= in_threshold_10_b));
        bool all_decreasing = ((in_threshold_01_a >= in_threshold_01_b) &&
                               (in_threshold_01_b >= in_threshold_10_a) &&
                               (in_threshold_10_a >= in_threshold_10_b));
        if (! (all_increasing || all_decreasing))
          throw invalid_argument("");
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by 4 numbers\n"
                       "       (These numbers must be either in increasing or decreasing order.)\n");
      }
      num_arguments_deleted = 5;
    }



    else if ((vArgs[i] == "-thresh-interval") ||
             (vArgs[i] == "-thresh-interval-out")) {
      try {
        if (i+2 >= vArgs.size())
          throw invalid_argument("");
        use_intensity_map = true;
        use_dual_thresholds = true;
        in_threshold_01_a = stof(vArgs[i+1]);
        in_threshold_01_b = stof(vArgs[i+1]);
        in_threshold_10_a = stof(vArgs[i+2]);
        in_threshold_10_b = stof(vArgs[i+2]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by 4 numbers.\n");
      }
      num_arguments_deleted = 3;
    }


    else if ((vArgs[i] == "-thresh-gauss") ||
             (vArgs[i] == "-thresh-gauss-out")) {
      try {
        if (i+2 >= vArgs.size())
          throw invalid_argument("");
        use_intensity_map = true;
        use_gauss_thresholds = true;
        out_thresh_gauss_x0 = stof(vArgs[i+1]);
        out_thresh_gauss_sigma = stof(vArgs[i+2]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by 4 numbers.\n");
      }
      num_arguments_deleted = 3;
    }





    else if (vArgs[i] == "-median")
    {
      #ifndef CXX17_UNSUPPORTED
      try {
        if ((i+1 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        median_radius = stof(vArgs[i+1]);
        filter_type = MEDIAN;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a nonnegative number\n");
      }
      num_arguments_deleted = 2;
      #else //#ifndef CXX17_UNSUPPORTED
      throw InputErr("Error: The " + vArgs[i] +
                     " argument is not supported\n"
                     "       by this version of " + g_program_name + ".\n"
                     "       You must obtain a fully C++17 compliant C++ compiler, and\n"
                     "       edit your compiler settings (eg, in the \"setup_...\" files), delete the\n"
                     "       \"-DCXX17_UNSUPPORTED\" compiler flag, and recompile the source code.\n"
                     "       This feature has never been tested and it might not work.\n");
      #endif
    } //if (vArgs[i] == "-median")



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



    #ifndef DISABLE_DOGGXY
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
    #endif //#ifndef DISABLE_DOGGXY



    else if ((vArgs[i] == "-log") ||
             (vArgs[i] == "-log-r") ||
             (vArgs[i] == "-log-d"))
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        float blob_sigma_multiplier = 1.0;
        blob_sigma_multiplier = 1.0;
        if (vArgs[i] == "-log")
          blob_sigma_multiplier = 1.0;
        if (vArgs[i] == "-log-r")
          // (for a solid uniform 3-D sphere)
          blob_sigma_multiplier = 1.0/sqrt(3.0);
        if (vArgs[i] == "-log-d")
          // (for a solid uniform 3-D sphere)
          blob_sigma_multiplier = 1.0/(2.0*sqrt(3.0));
        log_width[0] = stof(vArgs[i+1]) * blob_sigma_multiplier;
        log_width[1] = log_width[0];
        log_width[2] = log_width[0];
        m_exp = 2.0;
        n_exp = 2.0;
        filter_type = LOG_DOG;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a positive number:\n"
                       "       the (approximate) width of the objects of interest.\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-log")



    else if (vArgs[i] == "-log-aniso")
    {
      try {
        if ((i+3 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-'))
          throw invalid_argument("");
        log_width[0] = stof(vArgs[i+1]);
        log_width[1] = stof(vArgs[i+2]);
        log_width[2] = stof(vArgs[i+3]);
        m_exp = 2.0;
        n_exp = 2.0;
        filter_type = LOG_DOG;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by 3 positive numbers:\n"
                       "       the approximate widths of the objects of interest in the XYZ directions\n"
                       "       (Due to the missing-wedge artifact they might differ.\n");
      }
      num_arguments_deleted = 4;
    } //if (vArgs[i] == "-log-aniso")







    else if (vArgs[i] == "-dog-delta") {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        delta_sigma_over_sigma = stof(vArgs[i+1]);
        user_set_delta_manually = true;
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
        user_set_exponents_manually = true;
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
        #ifndef DISABLE_TEMPLATE_MATCHING
        template_background_exponent = n_exp;
        #endif
        user_set_exponents_manually = true;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a positive number.\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-gauss-exponent")



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



    else if ((vArgs[i] == "-radial-separation") ||
             (vArgs[i] == "-blob-separation") ||
             (vArgs[i] == "-blob-r-separation") ||
             (vArgs[i] == "-blobr-separation") ||
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
    } //if (vArgs[i] == "-radial-separation")



    else if ((vArgs[i] == "-blob") ||
             (vArgs[i] == "-blob-sigma") ||
             (vArgs[i] == "-blob-s") ||
             (vArgs[i] == "-blobs") ||
             (vArgs[i] == "-blob-radii") ||
             (vArgs[i] == "-blob-r") ||
             (vArgs[i] == "-blobr") ||
             (vArgs[i] == "-blob-diameters") ||
             (vArgs[i] == "-blob-d") ||
             (vArgs[i] == "-blob"))
    {
      try {
        if ((i+5 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-') ||
            (vArgs[i+4] == "") || (vArgs[i+4][0] == '-') ||
            (vArgs[i+5] == "") || (vArgs[i+5][0] == '-'))
          throw invalid_argument("");

        string blob_file_name_base = vArgs[i+2];
        //if (EndsWith(blob_file_name_base, ".txt"))
        //  blob_file_name_base =
        //    blob_file_name_base.substr(0,
        //                               blob_file_name_base.length()-4);

        if ((vArgs[i+1] == "minima") || (vArgs[i+1] == "min")) {
          blob_minima_file_name = blob_file_name_base;
          blob_maxima_file_name = "";  // disable searching for maxima
          score_upper_bound = 0.0;     // disable searching for maxima
        }
        else if ((vArgs[i+1] == "maxima") || (vArgs[i+1] == "max")) {
          blob_maxima_file_name = blob_file_name_base;
          blob_minima_file_name = "";   // disable searching for minima
          score_lower_bound = 0.0;      // disable searching for maxima
        }
        else if (vArgs[i+1] == "all") {
          blob_minima_file_name = blob_file_name_base + string(".minima.txt");
          blob_maxima_file_name = blob_file_name_base + string(".maxima.txt");
          if (score_lower_bound == 0.0)
            // if searching for minima was disabled, then enable it
            score_lower_bound = -std::numeric_limits<float>::infinity();
          if (score_upper_bound == 0.0)
            // if searching for maxima was disabled, then enable it
            score_upper_bound = std::numeric_limits<float>::infinity();
        }
        else {
          throw InputErr("Error: The 1st parameter to the \"" + vArgs[i] + "\" argument must be one of:\n"
                         "            \"minima\", \"maxima\", (or) \"all\"\n"
                         "\n"
                         "       (It indicates whether you are try to detect\n"
                         "        dark blobs, bright blobs, or both.)\n");
        }

        float blob_width_min = stof(vArgs[i+3]);
        float blob_width_max = stof(vArgs[i+4]);
        if ((blob_width_min <= 0.0) ||
            (blob_width_max <= 0.0) ||
            (blob_width_min >= blob_width_max))
          throw invalid_argument("");

        // The old version of the code allows user to set "N" directly
        //int N = stof(vArgs[i+5]);
        //if (N < 3)
        //  throw InputErr("The 4th parameter to the \"" + vArgs[i] + "\" argument must be >= 3.");
        //float growth_ratio = pow(blob_width_max / blob_width_min, 1.0/(N-1));

        // In the new version of the code, the user specifies "growth_ratio"
        float growth_ratio = stof(vArgs[i+5]);
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
        //    (! user_set_delta_manually))
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
        else if ((vArgs[i] == "-blob-diameters") || (vArgs[i] == "-blob-d"))
          // (for a solid uniform 3-D sphere)
          blob_width_multiplier = 1.0;
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
                       " argument must be followed by a type (either\n"
                       "       \"minima\", \"maxima\", or \"all\")\n"
                       "       as well as a file name, followed by 3 positive numbers.\n"
                       "For example:\n"
                       "           minima file_name.txt 200.0 280.0 1.02\n");
      }
      num_arguments_deleted = 6;
    } //if (vArgs[i] == "-blob")




    else if ((vArgs[i] == "-discard-blobs") ||
             (vArgs[i] == "-blob-nonmax") ||
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
      filter_type = BLOB_NONMAX_SUPPRESSION;
      num_arguments_deleted = 3;
    } //if (vArgs[i] == "-discard-blobs")




    else if (vArgs[i] == "-auto-thresh")
    {
      auto_thresh_score = true;
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == ""))
          throw invalid_argument("");
        string list_of_attributes = vArgs[i+1];
        if (list_of_attributes != "score")
          throw invalid_argument("");
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by \"score\"\n"
                       "       (In the future, it may be followed by a list of attributes surrounded\n"
                       "        in quotes.  But currently only the \"score\" attribute is supported.\n"
                       "        Thresholds for these attributes will be determined automatically.)\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-auto-thresh")



    else if (vArgs[i] == "-supervised")
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") ||
            (i+2 >= vArgs.size()) || (vArgs[i+2] == ""))
          throw invalid_argument("");
        training_pos_fname = vArgs[i+1];
        training_neg_fname = vArgs[i+2];
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by two file names:\n"
                       "       FILE_POS.txt  FILE_NEG.txt\n"
                       "       containing positive and negative training data, respectively.\n");
      }
      num_arguments_deleted = 3;
    } //if (vArgs[i] == "-supervised")



    else if (vArgs[i] == "-supervised-multi")
    {
      filter_type = BLOB_NONMAX_SUPERVISED_MULTI;
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == ""))
          throw invalid_argument("");
        string fname = vArgs[i+1];

        { //Read multi_training_meta_fname
          fstream f;
          f.open(fname, ios::in);
          const char comment_char = '#';
          string strLine;
          size_t i_line = 1;
          while (getline(f, strLine))
          {
            // ignore text after comments
            size_t ic = strLine.find(comment_char);
            if (ic != string::npos)
              strLine = strLine.substr(0, ic);
            stringstream ssLine(strLine);
            vector<string> tokens;
            string x;
            try {
              while (ssLine >> x) {
                tokens.push_back(x);
              }
            } // try {
            catch ( ... ) {
              stringstream err_msg;
              err_msg << "Error: Read error (invalid entry?) in file:\n"
                      << "       \"" << fname << "\", line: " << i_line << "\n";
              throw VisfdErr(err_msg.str());
            }
            if (tokens.size() == 0)
              continue;
            else if (tokens.size() != 3) {
              stringstream err_msg;
              err_msg << "Error: Format error in file:\n"
                      << "       \"" << fname << "\", line: " << i_line << "\n"
                      << "\n"
                      << "       Expected 3 file names in every line in this file:\n"
                      << "\n"
                      << "       training_pos_file  training_neg_file  blob_info_list.txt\n";
              throw VisfdErr(err_msg.str());
            }
            multi_training_pos_fnames.push_back(tokens[0]);
            multi_training_neg_fnames.push_back(tokens[1]);
            multi_in_coords_file_names.push_back(tokens[2]);
          } //while (getline(f, strLine))
          f.close();
        } //Read multi_training_meta_fname
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed a file name. Example: FILE_TRAINING_SETS.txt\n"
                       "\n"
                       "File format details:\n"
                       " This should be a 3-column text file containing file names.  Example:\n"
                       "   training_pos_1.txt  training_crds_neg_1.txt  blob_info_1.txt\n"
                       "   training_pos_2.txt  training_crds_neg_2.txt  blob_info_2.txt\n"
                       "   training_pos_3.txt  training_crds_neg_3.txt  blob_info_3.txt\n"
                       "                :                        :                :\n"
                       "   training_pos_N.txt  training_crds_neg_N.txt  blob_info_N.txt\n"
                       "\n"
                       " containing positive and negative training data, and a list of detected blobs\n"
                       " for each of the N images you want to simultaneously analyze.\n");
      }

      num_arguments_deleted = 2;

    } //if (vArgs[i] == "-supervised-multi")



    else if ((vArgs[i] == "-minima-threshold") ||
             (vArgs[i] == "-score-upper-bound"))
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == ""))
          throw invalid_argument("");
        score_upper_bound = stof(vArgs[i+1]);
        //score_lower_bound = -std::numeric_limits<float>::infinity();
        score_bounds_are_ratios = false;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number number.\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-minima-threshold")



    else if ((vArgs[i] == "-maxima-threshold") ||
             (vArgs[i] == "-score-lower-bound"))
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == ""))
          throw invalid_argument("");
        score_lower_bound = stof(vArgs[i+1]);
        //score_upper_bound = std::numeric_limits<float>::infinity();
        score_bounds_are_ratios = false;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number number.\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-maxima-threshold")



    else if ((vArgs[i] == "-minima-ratio") ||
             (vArgs[i] == "-score-lower-bound-ratio"))
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == ""))
          throw invalid_argument("");
        score_upper_bound = stof(vArgs[i+1]);
        //score_lower_bound = -std::numeric_limits<float>::infinity();
        score_bounds_are_ratios = true;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number number.\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-minima_ratio")



    else if ((vArgs[i] == "-maxima-ratio") ||
             (vArgs[i] == "-score-upper-bound-ratio"))
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == ""))
          throw invalid_argument("");
        score_lower_bound = stof(vArgs[i+1]);
        //score_upper_bound = std::numeric_limits<float>::infinity();
        score_bounds_are_ratios = true;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number number.\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-maxima-ratio")



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
        //masked_voxel_brightness = 1.0;
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
        masked_voxel_brightness = 0.0;
        specify_masked_brightness = true;
        width_a[0] = stof(vArgs[i+1]);
        width_a[1] = width_a[0];
        width_a[2] = width_a[0];
        float r = stof(vArgs[i+2]);
        template_background_radius[0] = r;
        template_background_radius[1] = r;
        template_background_radius[2] = r;
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
        masked_voxel_brightness = 0.0;
        specify_masked_brightness = true;
        width_a[0] = stof(vArgs[i+1]);
        width_a[1] = stof(vArgs[i+2]);
        width_a[2] = stof(vArgs[i+3]);
        template_background_radius[0] = stof(vArgs[i+4]);
        template_background_radius[1] = stof(vArgs[i+5]);
        template_background_radius[2] = stof(vArgs[i+6]);
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
        masked_voxel_brightness = 0.0;
        specify_masked_brightness = true;
        template_background_radius[0] = stof(vArgs[i+1]);
        template_background_radius[1] = stof(vArgs[i+2]);
        template_background_radius[2] = stof(vArgs[i+3]);
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
        masked_voxel_brightness = 0.0;
        specify_masked_brightness = true;
        template_background_radius[0] = stof(vArgs[i+1]);
        template_background_radius[1] = template_background_radius[0];
        template_background_radius[2] = template_background_radius[0];
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
      filter_type = FIND_EXTREMA;
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
      filter_type = FIND_EXTREMA;
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


    else if (vArgs[i] == "-neighbor-connectivity") {
      try {
        if (i+1 >= vArgs.size())
          throw invalid_argument("");
        neighbor_connectivity = stoi(vArgs[i+1]);
        if (neighbor_connectivity <= 0)
          throw invalid_argument("");
      }
      catch (invalid_argument& exc) {
          throw InputErr("Error: The " + vArgs[i] +
                         " argument must be followed by a positive integer.\n");
      }
      num_arguments_deleted = 2;
    }

    
    else if (vArgs[i] == "-boundary-extrema") {
      extrema_on_boundary = true;
      num_arguments_deleted = 1;
    }


    else if (vArgs[i] == "-ignore-boundary-extrema") {
      extrema_on_boundary = false;
      num_arguments_deleted = 1;
    }


    else if (vArgs[i] == "-distance-points") {
      if (i+1 >= vArgs.size())
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a file name.\n");
      filter_type = DISTANCE_TO_POINTS;
      in_coords_file_name = vArgs[i+1];
      num_arguments_deleted = 2;
    }


    else if ((vArgs[i] == "-draw-spheres") ||
             (vArgs[i] == "-spheres")) {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        filter_type = DRAW_SPHERES;
        in_coords_file_name = vArgs[i+1];
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a file name\n");
      }
      num_arguments_deleted = 2;
    }

    else if (vArgs[i] == "-draw-hollow-spheres") {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        filter_type = DRAW_SPHERES;
        in_coords_file_name = vArgs[i+1];
        if (! user_set_thickness_manually) {
          sphere_decals_shell_thickness = 0.05;
          sphere_decals_shell_thickness_is_ratio = true;
          sphere_decals_shell_thickness_min = 1.0;
        }
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a file name\n");
      }
      num_arguments_deleted = 2;
    }

    else if ((vArgs[i] == "-diameters") ||
             (vArgs[i] == "-diameter") ||
             (vArgs[i] == "-sphere-diameter") ||
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


    else if ((vArgs[i] == "-radii") ||
             (vArgs[i] == "-radius") ||
             (vArgs[i] == "-sphere-radius") ||
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


    else if ((vArgs[i] == "-sphere-shell-ratio") ||
             (vArgs[i] == "-spheres-shell-ratio")) {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        sphere_decals_shell_thickness_is_ratio = true;
        sphere_decals_shell_thickness     = stof(vArgs[i+1]);
        user_set_thickness_manually = true;
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
        user_set_thickness_manually = true;
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
        user_set_thickness_manually = true;
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

    else if ((vArgs[i] == "-background") ||
             (vArgs[i] == "-spheres-background") ||
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


    else if ((vArgs[i] == "-background-scale") ||
             (vArgs[i] == "-spheres-background-scale") ||
             (vArgs[i] == "-sphere-background-scale")) {
      //sphere_decals_background_use_orig = true;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        sphere_decals_background_scale = stof(vArgs[i+1]);
        user_set_background_scale_manually = true;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number, usually between 0 and 1:\n"
                       "       how much to supress fluctuations in the original background image.\n");
      }
      num_arguments_deleted = 2;
    }


    else if ((vArgs[i] == "-foreground") ||
             (vArgs[i] == "-spheres-foreground") ||
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

    else if (vArgs[i] == "-background-auto") {
      sphere_decals_background_norm = true;
      if (! user_set_background_scale_manually) {
        sphere_decals_background_scale = 0.35;
      }
      num_arguments_deleted = 1;
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
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        if ((vArgs[i+1] == "min") || (vArgs[i+1] == "minima")) {
          clusters_begin_at_maxima = false;
          if (! user_set_watershed_threshold_manually)
            watershed_threshold = std::numeric_limits<float>::infinity();
        }
        else if ((vArgs[i+1] == "max") || (vArgs[i+1] == "maxima")) {
          clusters_begin_at_maxima = true;
          if (! user_set_watershed_threshold_manually)
            watershed_threshold = -std::numeric_limits<float>::infinity();
        }
        else
          throw invalid_argument("");
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by an argument:  \"type\"  \"width\"\n"
                       "       The \"type\" argument must be either \"minima\" or \"maxima\".\n"
                       "       (It depends on whether you want to detect dark or bright objects.)\n");
      }
      num_arguments_deleted = 2;
    }



    else if (vArgs[i] == "-watershed-threshold")
    {
      filter_type = WATERSHED;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == ""))
          throw invalid_argument("");
        user_set_watershed_threshold_manually = true;
        watershed_threshold = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number\n");
      }      
      num_arguments_deleted = 2;
    }



    else if (vArgs[i] == "-watershed-show-boundaries")
    {
      filter_type = WATERSHED;
      watershed_show_boundaries = true;
      num_arguments_deleted = 1;
    }



    else if (vArgs[i] == "-watershed-hide-boundaries")
    {
      filter_type = WATERSHED;
      watershed_show_boundaries = false;
      num_arguments_deleted = 1;
    }



    else if (vArgs[i] == "-watershed-boundary")
    {
      filter_type = WATERSHED;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == ""))
          throw invalid_argument("");
        watershed_boundary_label = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number\n");
      }
      num_arguments_deleted = 2;
    }

    

    else if (vArgs[i] == "-markers")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == ""))
          throw invalid_argument("");
        watershed_markers_filename = vArgs[i+1];
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by an image file name\n");
      }      
      num_arguments_deleted = 2;
    }



    else if (vArgs[i] == "-undefined-out")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == ""))
          throw invalid_argument("");
        if (vArgs[i+1] == "max")
          undefined_voxels_are_max = true;
        else {
          undefined_voxels_are_max = false;
          undefined_voxel_brightness = stof(vArgs[i+1]);
        }
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number or \"max\"\n");
      }
      num_arguments_deleted = 2;
    }



    else if ((vArgs[i] == "-surface") ||
             (vArgs[i] == "-planar"))
    {
      throw InputErr("Error: As of 2021-7-15, the " + vArgs[i] +          
                     " argument has been renamed.\n"
                     "       It is now called \"-membrane\".\n"
                     "       See documentation for details.\n");
    }



    else if (vArgs[i] == "--membrane-normals-file") {
      throw InputErr("Error: As of 2019-4-11, the " + vArgs[i] +
                     " argument has been renamed.\n"
                     "       It is now called \"-normals-file\".\n"
                     "       See documentation for details.\n");
    }



    else if (vArgs[i] == "-planar-tv") {
      throw InputErr("Error: As of 2019-4-11, the " + vArgs[i] +
                     " argument has been renamed.\n"
                     "       It is now called \"-tv\"\n");
    }



    else if ((vArgs[i] == "-membrane") || (vArgs[i] == "-surface-ridge") ||
             (vArgs[i] == "-edge") || (vArgs[i] == "-surface-edge") ||
             (vArgs[i] == "-curve"))
    {
      if ((vArgs[i] == "-membrane") ||
          (vArgs[i] == "-surface-ridge"))
        filter_type = SURFACE_RIDGE;
      else if ((vArgs[i] == "-edge") || (vArgs[i] == "-surface-edge"))
        filter_type = SURFACE_EDGE;
      else if (vArgs[i] == "-curve")
        filter_type = CURVE;

      try {
        if ((i+2 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-'))
          throw invalid_argument("");
        if ((vArgs[i+1] == "min") || (vArgs[i+1] == "minima"))
          ridges_are_maxima = false;
        else if ((vArgs[i+1] == "max") || (vArgs[i+1] == "maxima"))
          ridges_are_maxima = true;
        else
          throw invalid_argument("");
        float thickness = stof(vArgs[i+2]);
        //max_distance_to_feature = 1.5*thickness;

        // Now choose the "sigma" parameter of the Gaussians we will use:
        float sigma;
        switch (filter_type) {
        case SURFACE_EDGE:
          sigma = thickness;
          break;
        case SURFACE_RIDGE:
          sigma = thickness / sqrt(3.0);
          break;
        case CURVE:
          sigma = thickness / sqrt(3.0);  //CONTINUEHERE: <-- IS THIS CORRECT?
          break;
        default:
          assert(false);
          break;
        }
        width_a[0] = sigma;
        width_a[1] = width_a[0];
        width_a[2] = width_a[0];
        width_b[0] = 0.0;
        width_b[1] = 0.0;
        width_b[2] = 0.0;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by 2 arguments:  \"type\"  \"width\"\n"
                       "       The \"type\" argument must be either \"minima\" or \"maxima\"\n"
                       "       (depending on whether you want to detect dark or bright objects)\n"
                       "       The \"width\" argument should be target thickness of the thin\n"
                       "       surface-like object you are interest in detecting.\n");
      }
      num_arguments_deleted = 3;
    }


    else if ((vArgs[i] == "-detection-background") ||
             (vArgs[i] == "-membrane-background") ||
             (vArgs[i] == "-curve-background"))
    {
      filter_type = SURFACE_RIDGE;
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
                       "       (This selects the low frequency variations in the background that\n"
                       "        you wish to ignore.)\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-detection-threshold") {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        hessian_score_threshold = stof(vArgs[i+1]);
        hessian_score_threshold_is_a_fraction = false;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number, the threshold\n"
                       "       for deciding whether to compute the surface normal there.\n");
      }
      num_arguments_deleted = 2;
    }


    else if ((vArgs[i] == "-tv-best") ||
             (vArgs[i] == "-best"))
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        hessian_score_threshold = stof(vArgs[i+1]);
        hessian_score_threshold_is_a_fraction = true;
        if (! ((0.0 <= hessian_score_threshold) &&
               (hessian_score_threshold <= 1.0)))
          throw invalid_argument("");
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number between 0 and 1:  The fraction of the\n"
                       "       number of voxels you wish to keep (not discard) after ridge detection\n"
                       "       (giving preference to the most \"ridge-like\" voxels).\n"
                       "       The \"ridgyness\" (saliency) of other voxels will be set to 0\n"
                       "\n"
                       "       If there are additional steps in the computation (such as tensor voting\n"
                       "       planar templates or ant-colony-optimization), then\n"
                       "       the computational cost of these subsequent steps will be proportional\n"
                       "       to this number.  It must lie in the range from 0 to 1.  (Default: 1)\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-tv")
    {
      if ((filter_type != SURFACE_RIDGE) &&
          (filter_type != SURFACE_EDGE) &&
          (filter_type != CURVE))
        throw InputErr("Error: The -tv argument must appear in the argument list after the\n"
                       "       -membrane, -edge, or -curve arguments.  (It does not make sense to use\n"
                       "       tensor voting unless you are using one of these directional detectors.)\n");
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        tv_sigma = stof(vArgs[i+1]);
        //if (tv_sigma < 1.0)
        //  throw InputErr("Error: the -tv argument must be >= 1.0\n");
        //  (This is probably a good idea, but I no longer require it.)
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number, the ratio\n"
                       "       of widths used in directional blurring (typically 3)\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-tv-angle-exponent")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        tv_exponent = stoi(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number, the threshold\n"
                       "       for deciding whether to compute the surface normal there.\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-tv-threshold")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == ""))
          throw invalid_argument("");
        tv_score_threshold = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number, the threshold\n"
                       "       for deciding whether to compute the surface normal there.\n");
      }
      num_arguments_deleted = 2;
    }

    else if (vArgs[i] == "-tv-truncate-ratio")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        tv_truncate_ratio = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a positive number: the\n"
                       "       distance (in units σ_tv) over which tensor votes are cast. (See docs.)\n");
      }
      num_arguments_deleted = 2;
    }

    // CRUFT
    //   I can't remember what this code was supposed to be used for.
    // COMMENTING OUT:
    //else if (vArgs[i] == "-tv-iter") {
    //  try {
    //    if ((i+1 >= vArgs.size()) ||
    //        (vArgs[i+1] == ""))
    //      throw invalid_argument("");
    //    tv_num_iters = stof(vArgs[i+1]);
    //  }
    //  catch (invalid_argument& exc) {
    //    throw InputErr("Error: The " + vArgs[i] +
    //                   " argument must be followed by a number.");
    //  }
    //  num_arguments_deleted = 2;
    //}


    else if ((vArgs[i] == "-normals-file") ||
             (vArgs[i] == "-surface-normals-file"))
    {
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



    else if ((vArgs[i] == "-max-voxels-to-feature") ||
             (vArgs[i] == "-max-voxels-to-surface") ||
             (vArgs[i] == "-max-voxels-to-membrane") ||
             (vArgs[i] == "-max-voxels-to-edge") ||
             (vArgs[i] == "-max-voxels-to-curve"))
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        if ((vArgs[i+1] == "inf") ||
            (vArgs[i+1] == "infinity") ||
            (vArgs[i+1] == "disable"))
          // max_distance_to_feature = std::numeric_limits<float>::infinity()
          max_distance_to_feature = 0.0; // either 0 or infinity() disables this
        else
          max_distance_to_feature = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a positive number\n");
      }
      num_arguments_deleted = 2;
    }



    else if ((vArgs[i] == "-max-distance-to-feature") ||
             (vArgs[i] == "-max-distance-to-surface") ||
             (vArgs[i] == "-max-distance-to-membrane") ||
             (vArgs[i] == "-max-distance-to-edge") ||
             (vArgs[i] == "-max-distance-to-curve"))
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        if ((vArgs[i+1] == "inf") ||
            (vArgs[i+1] == "infinity") ||
            (vArgs[i+1] == "disable"))
          // max_distance_to_feature = std::numeric_limits<float>::infinity()
          max_distance_to_feature = 0.0; // either 0 or infinity() disables this
        else
          max_distance_to_feature = -stof(vArgs[i+1]); //(will flip sign later)
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a positive number\n");
      }
      num_arguments_deleted = 2;
    }



    else if ((vArgs[i] == "-connect") ||
             (vArgs[i] == "-connect-bright") ||
             (vArgs[i] == "-connect-saliency"))
    {
      cluster_connected_voxels = true;
      clusters_begin_at_maxima = true;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == ""))
          throw invalid_argument("");
        connect_threshold_saliency = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-connect-dark")
    {
      cluster_connected_voxels = true;
      clusters_begin_at_maxima = false;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == ""))
          throw invalid_argument("");
        connect_threshold_saliency = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-connect-angle")
    {
      cluster_connected_voxels = true;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        double theta = stof(vArgs[i+1]);
        connect_threshold_vector_saliency = cos(theta*M_PI/180.0);
        connect_threshold_vector_neighbor = cos(theta*M_PI/180.0);
        connect_threshold_tensor_saliency = cos(theta*M_PI/180.0);
        connect_threshold_tensor_neighbor = cos(theta*M_PI/180.0);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a nonnegative number\n");
      }
      num_arguments_deleted = 2;
    }


    else if ((vArgs[i] == "-connect-vector-saliency") ||
             (vArgs[i] == "-cvs"))
    {
      cluster_connected_voxels = true;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        connect_threshold_vector_saliency = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a nonnegative number\n");
      }
      num_arguments_deleted = 2;
    }

    else if ((vArgs[i] == "-connect-vector-neighbor") ||
             (vArgs[i] == "-cvn"))
    {
      cluster_connected_voxels = true;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        connect_threshold_vector_neighbor = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a nonnegative number\n");
      }
      num_arguments_deleted = 2;
    }

    else if ((vArgs[i] == "-connect-tensor-saliency") ||
             (vArgs[i] == "-cts"))
    {
      cluster_connected_voxels = true;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        connect_threshold_tensor_saliency = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a nonnegative number\n");
      }
      num_arguments_deleted = 2;
    }

    else if ((vArgs[i] == "-connect-tensor-neighbor") ||
             (vArgs[i] == "-ctn"))
    {
      cluster_connected_voxels = true;
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        connect_threshold_tensor_neighbor = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a nonnegative number\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-select-cluster")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        select_cluster = stoi(vArgs[i+1]);
        if (select_cluster < 0)
          throw invalid_argument("");
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a nonnegative integer\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-must-link")
    {
      cluster_connected_voxels = true;
      try
      {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == ""))
          throw invalid_argument("");
        must_link_filename = vArgs[i+1];
      } // try {
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed by a number\n");
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


    else if ((vArgs[i] == "-blob-intensity-vs-radius") ||
             (vArgs[i] == "-blob-radial-intensity"))
    {
      #ifndef DISABLE_INTENSITY_PROFILES
      try {
        if ((i+3 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-') ||
            (vArgs[i+2] == vArgs[i+3]))
          throw invalid_argument("");
        if ((vArgs[i+1] == "min") || (vArgs[i+1] == "minima"))
          blob_profiles_center_criteria = BlobCenterCriteria::MINIMA;
        else if ((vArgs[i+1] == "max") || (vArgs[i+1] == "maxima"))
          blob_profiles_center_criteria = BlobCenterCriteria::MAXIMA;
        else if ((vArgs[i+1] == "center") || (vArgs[i+1] == "cen"))
          blob_profiles_center_criteria = BlobCenterCriteria::CENTER;
        else
          throw invalid_argument("");
        in_coords_file_name = vArgs[i+2];
        blob_profiles_file_name_base = vArgs[i+3];
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] +
                       " argument must be followed\n"
                       "       by 3 additional arguments:\n"
                       "       CENTER_TYPE  input_coords_file  output_file_base_name \n"
                       "       where \"CENTER_TYPE\" is one of these 3 choices: \"min\",\"max\",\"center\",\n"
                       "       and the other two arguments are file names.\n");
      }
      filter_type = BLOB_RADIAL_INTENSITY;
      num_arguments_deleted = 4;
      #else
      throw InputErr("Error: The " + vArgs[i] + " argument is not supported by this\n"
                     "       version of " + g_program_name + ".  Recompile the code to enable this feature.\n");
      #endif //#ifndef DISABLE_INTENSITY_PROFILES
    } // else if (vArgs[i] == "-blob-intensity-vs-radius")


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
    //  masked_voxel_brightness = 0.0;
    //  specify_masked_brightness = true;
    //  width_a[0] = stof(vArgs[i+1]);
    //  width_a[1] = width_a[0];
    //  width_a[2] = width_a[0];
    //  m_exp = stof(vArgs[i+2]);
    //  float r = stof(vArgs[i+3]);
    //  template_background_radius[0] = r;
    //  template_background_radius[1] = r;
    //  template_background_radius[2] = r;
    //  template_background_exponent = stof(vArgs[i+4]);
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
      ((in_set_image_size[0] == 0) ||
       (in_set_image_size[1] == 0) ||
       (in_set_image_size[2] == 0))
      &&
      (filter_type != BLOB_NONMAX_SUPPRESSION))
  {
    throw InputErr("Error: Unknown image size.\n"
                   "       You can either specify the name of the image file you want to read\n"
                   "       using the \"-in SOURCE_FILE\" argument, ...\n"
                   "         (these files usually end in \".mrc\" or \".rec\")\n"
                   "       Alternatively you can start with a blank image using the\n"
                   "       \"-image-size nx ny nz\" argument.  Either way, you must\n"
                   "       indicate the size of the image that you want to create or analyze.\n");
  }
  if ((out_file_name.size() == 0) &&
      ((!
        ((out_normals_fname != "") ||
         (filter_type == NONE) ||
         (filter_type == BLOB) ||
         (filter_type == FIND_EXTREMA) ||
         (filter_type == BLOB_NONMAX_SUPPRESSION) ||
         (filter_type == BLOB_NONMAX_SUPERVISED_MULTI) ||
         (filter_type == BLOB_RADIAL_INTENSITY)))
       ||
       use_intensity_map ||
       invert_output))
  {
    throw InputErr("Error: You must specify the name of the tomogram you want to create\n"
                   "       using the \"-out DESTINATION_FILE\" argument\n"
                   "       (These files usually end in \".mrc\" or \".rec\" .)\n");
  }
  if ((out_file_name == in_file_name) && (out_file_name.size() != 0) &&
      (! out_file_overwrite))
  {
    throw InputErr("Error: Input and Output image files cannot be the same.  (The purpose of this\n"
                   "       error is to prevent erasing your original source image by accident.)\n"
                   "\n"
                   "       To override this, use \"-outf\" instead of \"-out\".\n");
  }


  // REMOVE THIS CRUFT
  // ----------
  //if ((filter_type == GAUSS) || (filter_type == GGAUSS)) {
  //  // If the exponent equals 2, then the function we are convolving
  //  // with is an ordinary Gaussians.  Ordinary Gaussians are are separable
  //  // and can be convolved successively in 1-D in the x,y,z directions.
  //  // Hence we can save a great deal of time by doing this compared to
  //  // using a more general filter function which would force us to perform
  //  // the full 3-D convolution.
  //  if (user_set_exponents_manually) {
  //    filter_type = GGAUSS;
  //    if (m_exp == 2.0) {
  //      filter_type = GAUSS; // <-- use fast seperable Gaussians
  //      width_a[0] /= sqrt(2.0);  // compensate for using
  //      width_a[1] /= sqrt(2.0);  //  exp(-0.5*(r/width_a)^2) instead of
  //      width_a[2] /= sqrt(2.0);  //  exp(-(r/width_a)^2)
  //    } //if (m_exp == 2.0)
  //  } //if (user_set_exponents_manually)
  //} //if ((filter_type == GAUSS) || (filter_type == GGAUSS))


  // REMOVE THIS CRUFT
  // ----------
  //if ((filter_type == DOG) || (filter_type == DOGG)) {
  //  // If the exponents equal 2, then the functions we are convolving
  //  // with are ordinary Gaussians.  Ordinary Gaussians are are separable
  //  // and can be convolved successively in 1-D in the x,y,z directions.
  //  // Hence we can save a great deal of time by doing this compared to
  //  // using a more general filter function which would force us to perform
  //  // the full 3-D convolution.
  //  if (user_set_exponents_manually) {
  //    filter_type = DOGG;
  //    if ((m_exp == 2.0) && (n_exp == 2.0)) {
  //      filter_type = DOG; // <-- use fast seperable Gaussians
  //
  //      width_a[0] /= sqrt(2.0);  // compensate for using
  //      width_a[1] /= sqrt(2.0);  //  exp(-0.5*(r/width_a)^2) instead of
  //      width_a[2] /= sqrt(2.0);  //  exp(-(r/width_a)^2)
  //      width_b[0] /= sqrt(2.0);
  //      width_b[1] /= sqrt(2.0);
  //      width_b[2] /= sqrt(2.0);
  //    } //if ((m_exp == 2.0) && (n_exp == 2.0))
  //    else
  //      filter_type = DOGG;
  //  } //if (user_set_exponents_manually) {
  //} //if ((filter_type == DOG) || (filter_type == DOGG))


  // REMOVE THIS CRUFT
  //if ((filter_type == TEMPLATE_GAUSS) || (filter_type == TEMPLATE_GGAUSS)) {
  //  // If the exponents equal 2, then the functions we are convolving
  //  // with are ordinary Gaussians.  Ordinary Gaussians are are separable
  //  // and can be convolved successively in 1-D in the x,y,z directions.
  //  // Hence we can save a great deal of time by doing this compared to
  //  // using a more general filter function which would force us to perform
  //  // the full 3-D convolution.
  //  if (user_set_exponents_manually) {
  //    filter_type = TEMPLATE_GGAUSS;
  //    if ((m_exp == 2.0) && (n_exp == 2.0)) {
  //      filter_type = TEMPLATE_GAUSS; // <-- use fast Gaussians instead of DOGG which is slow

  //      width_a[0] /= sqrt(2.0);  // compensate for using
  //      width_a[1] /= sqrt(2.0);  //  exp(-0.5*(r/width_a)^2) instead of
  //      width_a[2] /= sqrt(2.0);  //  exp(-(r/width_a)^2)
  //      width_b[0] /= sqrt(2.0);
  //      width_b[1] /= sqrt(2.0);
  //      width_b[2] /= sqrt(2.0);
  //      template_background_radius[0] /= sqrt(2.0);// compensate for using:
  //      template_background_radius[1] /= sqrt(2.0);//  exp(-0.5*(r/w)^2)
  //      template_background_radius[2] /= sqrt(2.0);// instead of exp(-(r/w)^2)
  //    } //if ((m_exp == 2.0) && (n_exp == 2.0))
  //    else
  //      filter_type = TEMPLATE_GGAUSS;
  //  } //if (user_set_exponents_manually) {
  //} //if ((filter_type == DOG) || (filter_type == DOGG))



  // REMOVE THIS CRUFT.
  //   (THIS CODE IS NO LONGER NECESSARY AFTER 2019-3-03)
  //if (filter_type == BLOB) {
  //  if ((blob_minima_file_name != "") &&
  //      (score_lower_bound == -std::numeric_limits<float>::infinity()))
  //  {
  //    // If the user is not interested in local maxima, then we only
  //    // need to create one file (for the minima, as opposed to a separate
  //    // file for both minima and maxima).  In that case, remove the
  //    // ".minima.txt" which was earlier automatically added to the end
  //    // of that file's name.
  //    string blob_file_name_base = blob_minima_file_name;
  //    if (EndsWith(blob_file_name_base, ".minima.txt"))
  //      blob_file_name_base =
  //        blob_file_name_base.substr(0,
  //                                   blob_file_name_base.length()-11);
  //    blob_minima_file_name = blob_file_name_base;
  //  }
  //  if ((blob_maxima_file_name != "") &&
  //      (score_upper_bound == std::numeric_limits<float>::infinity()))
  //  {
  //    // If the user is not interested in local minima, then we only
  //    // need to create one file (for the maxima, as opposed to a separate
  //    // file for both minima and maxima).  In that case, remove the
  //    // ".maxima.txt" which was earlier automatically added to the end
  //    // of that file's name.
  //    string blob_file_name_base = blob_maxima_file_name;
  //    if (EndsWith(blob_file_name_base, ".maxima.txt"))
  //      blob_file_name_base =
  //        blob_file_name_base.substr(0,
  //                                   blob_file_name_base.length()-11);
  //    blob_maxima_file_name = blob_file_name_base;
  //  }
  //} //if (filter_type == BLOB)



  if (filter_type == SURFACE_RIDGE) {
    float sigma = width_a[0];
    // The parameter entered by the user is a ratio, not an absolute number.
    // Now that we know what sigma is, multiply tv_sigma by sigma.
    tv_sigma *= sigma;
  }


  if (cluster_connected_voxels &&
      (filter_type != SURFACE_RIDGE) &&
      (filter_type != SURFACE_EDGE) &&
      (filter_type != CURVE))
  {
    assert(connect_threshold_saliency!=std::numeric_limits<float>::infinity());
    filter_type = LABEL_CONNECTED;
  }


  // Now read various files which the user asked us to read
  // and save the resulting information somewhere.

  if (training_pos_fname != "")
    // read a list of locations in the image from a file
    // (positive training data for supervised learning)
    is_training_pos_in_voxels =
      ReadCoordinates(training_pos_fname,
                      training_pos_crds,
                      '#');
  if (training_neg_fname != "")
    // read a list of locations in the image from a file
    // (negative training data for supervised learning)
    is_training_neg_in_voxels =
      ReadCoordinates(training_neg_fname,
                      training_neg_crds,
                      '#');

  // If there are multiple independent training set data files, read them too
  if ((multi_training_pos_fnames.size() > 0) ||
      (multi_training_neg_fnames.size() > 0))
  {
    if ((in_coords_file_name != "") &&
        (out_coords_file_name != ""))
      throw InputErr("Error: This program is too stupid to automatically discard blobs from multiple\n"
                     "       images in a single invocation.  For this reason,\n"
                     "       you may not use the \"-supervised-multi\" argument simultaneously with\n"
                     "       any other arguments which also perform nonmax-suppression of blobs\n"
                     "       (such as \"-discard-blobs\").  In particular, you must make sure that\n"
                     "       any overlapping blobs in any of the reference lists of blobs must be\n"
                     "       removed in advance.  (These blobs are contained in the files referred\n"
                     "       to in the 3rd column of the file supplied to \"-supervised-multi\".\n"

                     "       They can be removed by using this program with\"-discard-blobs\"\n"
                     "       together with \"-radial-separation\", or \"-max-volume-overlap\", and/or\n"
                     "       \"-max-volume-overlap-small\".  This must be done separately\n"
                     "       for each of those files in the 3rd column.)\n");
    if (mask_file_name != "")
      throw InputErr("Error: You may not use the \"-supervised-multi\" argument simultaneously with\n"
                     "       the \"-mask\" argument.  If you wish to ignore blobs which lie outside\n"
                     "       the mask, then you must delete the corresponding lines from the\n"
                     "       blob lists beforehand.  (These are stored in the files indicated\n"
                     "       by the 3rd column of the file supplied to \"-supervised-multi\"\n"
                     "       You can do this by running this program using the \"-discard-blobs\"\n"
                     "       together with the \"-mask\" argument.  This must be done separately\n"
                     
                     "       on each of these lists of blobs beforehand.)\n");


    int Nsets = multi_training_pos_fnames.size();
    multi_training_pos_crds.resize(Nsets);
    multi_is_training_pos_in_voxels.resize(Nsets);
    for (int I = 0; I < multi_training_pos_fnames.size(); ++I) {
      // read a list of locations in the image from a file
      // (positive training data for supervised learning)
      multi_is_training_pos_in_voxels[I] =
        ReadCoordinates(multi_training_pos_fnames[I],
                        multi_training_pos_crds[I],
                        '#');
    }
    multi_training_neg_crds.resize(Nsets);
    multi_is_training_neg_in_voxels.resize(Nsets);
    for (int I = 0; I < multi_training_neg_fnames.size(); ++I) {
      // read a list of locations in the image from a file
      // (negative training data for supervised learning)
      multi_is_training_neg_in_voxels[I] =
        ReadCoordinates(multi_training_neg_fnames[I],
                        multi_training_neg_crds[I],
                        '#');
    }
  } //if ((multi_training_pos_fnames.size() > 0) || ...


  if (must_link_filename != "")
    // read a list of locations that the user thinks should
    // belong to the same cluster. (useful for clustering)
    is_must_link_constraints_in_voxels =
      ProcessLinkConstraints(must_link_filename,
                             must_link_constraints,
                             must_link_constraint_directions);


  if (filter_type == SURFACE_EDGE) {
    // When clustering voxels together into larger blobs, we sometimes
    // consider the orientation of the hessian.  But this is not always useful.
    // If we are using an edge detector, then the hessian (2nd derivative)
    // of the brightness should be close to zero near the edge of a light-dark
    // boundary.  At these boundary locations, the gradient is large, but
    // the second derivative is approximately zero.  Consequently, we should
    // not be comparing the directions of the principal axes of the hessian with
    // the tensor axes, (again, ince the hessian there is approximately zero).
    // To disable that comparision, set these thresholds to zero:
    connect_threshold_vector_saliency = -std::numeric_limits<float>::infinity();
    connect_threshold_tensor_saliency = -std::numeric_limits<float>::infinity();
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
