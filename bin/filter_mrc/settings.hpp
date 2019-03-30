#ifndef _SETTINGS_HPP
#define _SETTINGS_HPP

#define DISABLE_BOOTSTRAPPING
//#define DISABLE_TEMPLATE_MATCHING


#include<vector>
#include<string>
using namespace std;
#include <mrc_simple.hpp>
#include <filter3d.hpp>
#include "filter3d_unsupported.hpp"


extern string g_program_name;
extern string g_version_string;
extern string g_date_string;


/// @brief  The "Settings" class contains a list of parameters specified by the
///         user which control the behavior of the "filter_mrc" program.

class Settings {
 public:

  Settings();

  // The user can select one of the following algorithms when running
  // this program.  (Some of the algorithms share the same parameters.)

  typedef enum eFilterType {
    NONE,    //(the user has not selected a filter, or is using thresholding)
    GAUSS,   //3D Gaussian filter
    GGAUSS,  //3D Generalized Gaussian filter
    DOG,     //3D Difference of Gaussians (arbitrary widths,
             //                            fast, seperable Gaussians)
    DOGG,    //3D Generaliws Difference-of-Generalized-Gaussians with arbitrary
             //   widths and arbitrary exponents, slow)
    DOGGXY,  //2D Generalized Difference-of-Generalized-Gaussians with arbitrary
             //   widths andarbitrary exponents, slow)
             //   and a 1D Gaussian filter (in the Z direction)
    DOG_SCALE_FREE,  // scale-free version of the DOG filter (similar widths)
    TEMPLATE_GAUSS,  // Perform a least-squares fit between a Gaussian
                     // and the surrounding voxel brightnesses
    TEMPLATE_GGAUSS, // Perform a least-squares fit between a generalized
                     // Gaussian and the surrounding voxel brightnesses
    BOOTSTRAP_DOGG,  // DOGG filter with bootstrapping (significance testing)
    BLOB,            // blob detection
    RIDGE_PLANAR,    // a ridge detector for surfaces (ie a membrane detector)
    WATERSHED,       // Watershed segmentation
    LOCAL_FLUCTUATIONS, // Report the fluctuation of nearby voxel intensities
    SPHERE_NONMAX_SUPPRESSION,//throw away overlapping objects(detected earlier)
    DISTANCE_TO_POINTS, //voxel intensity = distance to nearest object
    SPHERE_DECALS,   //voxel intensity = 1 if within R from the nearest object
    BLOB_RADIAL_INTENSITY // Report intensity-vs-radius for each blob
  } FilterType; 
  
  FilterType filter_type; // The user must select one of these filter types


  // ---- parameters describing the physical size of each voxel ----
  //      (shared by all filter types)

  float voxel_width;  //physical width of each voxel (for example in Angstroms)
  bool voxel_width_divide_by_10; //Use nm instead of Angstroms?


  // -- parameters which determine the image files we will be reading/writing --
  //      (shared by all filter types)
  
  string in_file_name; // name of the image file we want to read
  string out_file_name; // name of the image file we want to create
  // Mask parameters are used to select (ignore) voxels from the original image.
  string mask_file_name; // name of an image file used for masking
  bool use_mask_select; // do we select voxels with a specific value?
  int mask_select; // select only voxels from the input with this value
  bool use_mask_out; // should we specify the brightness of ignored voxels?
  int mask_out;//what brightness do we set these voxels?(to go in out_file_name)
  bool mask_rectangle_in_voxels; // Is the mask-region a rectangular box?
  float mask_rectangle_xmin; // if so, what is the shape of that box?
  float mask_rectangle_xmax; // :
  float mask_rectangle_ymin;
  float mask_rectangle_ymax;
  float mask_rectangle_zmin;
  float mask_rectangle_zmax;


  // ---- parameters for "difference of generalized gaussians" filters ----
  //
  // A significant number of the filters used by this program
  // are all of the type "difference of generalized gaussians".
  // (Specifically, these are: GAUSS, GGAUSS, DOG, DOGG, DOGGXY, DOG_SCALE_FREE,
  //                           TEMPLATE_GAUSS, TEMPLATE_GGAUSS, BOOTSTRAP_DOGG)
  // 
  // The settings for these filters all use the same parameters:
  // ("width_a", and "width_b", "m_exp", "n_exp")
  // They are described below:
  //
  //   h(x,y,z) = h_a(x,y,z) - h_b(x,y,z)
  //
  // where both h_a() and h_b() are generalized normal distributions (Gaussians)
  // https://en.wikipedia.org/wiki/Generalized_normal_distribution
  //
  // If b_x>0, b_y>0, b_z>0, then:
  //
  //   h_a(x,y,z) = A*exp(-r_a^m)
  //          r_a = sqrt((x/a_x)*2 + (y/a_y)*2 + (z/a_z)*2)
  //   h_b(x,y,z) = B*exp(-r_b^n)
  //          r_b = sqrt((x/b_x)*2 + (y/b_y)*2 + (z/b_z)*2)
  //
  // Note: The a_x, a_y, a_z, b_x, b_y, b_z parameters in this formula are
  //       represented by the "width_a" and "width_b" data members, and the
  //       "m" and "n" parameters are represented by "m_exp" and "n_exp"
  //       (See below.)
  //
  // The "A" and "B" parameters are chosen so that sum_{x,y,z} h(x,y,z) = 0
  // (In the special case below where h_b=0, then sum_{x,y,z} h(x,y,z) = 1)
  // We might wish to use a "difference of gaussians" approach in some
  // directions (eg Z) , and ordinary gaussians" in other directions (eg X,Y):
  // For each dimension where the corresponding "b" parameter is < 0,
  // the program uses a "separable" filter which is Gaussian in that dimension.
  //    (REALLY? I should check this. -Andrew 2019-2-23)
  // For example, If b_z < 0, then both h_a and h_b become:
  //             functions of x,y, multiplied by a Gaussian in the Z direction:
  //   h_a(x,y,z) = A * exp(-r_a^m)   *  C * exp(-(1/2)z^2/a_z)
  //          r_a = sqrt((x/a_x)*2 + (y/a_y)*2
  //   h_b(x,y,z) = B * exp(-r_b^n)   *  C * exp(-(1/2)z^2/a_z)
  //          r_b = sqrt((x/b_x)*2 + (y/b_y)*2
  //
  // Thus, a simple 3-D Gaussian (whose principle axes in the x,y,z directions)
  // can be implemented setting b_x = b_y = b_z = -1.0.
  //   h(x,y,z) =
  //     A * exp(-(1/2)(x/a_x)^2) * exp(-(1/2)(y/a_y)^2) * exp(-(1/2)(z/a_z)^2)
  // "Seperable" filters like this are faster to compute,
  // requiring time  O( Nx * Ny * Nz * max(a_x, a_y, a_z) )
  //     instead of  O( Nx * Ny * Nz * a_x * a_y * a_z) )
  //                 (where Nx, Ny, Nz denote the source image size in voxels)
  //
  // Special case:
  //       In the specific case where the exponents m = n = 2, this filter
  //       can be implemented as the sum of two terms, each 
  //       of which are seperable in x,y,z and can be computed quickly:
  //   h(x,y,z) = h_a(x,y,z) - h_b(x,y,z)
  //   h_a(x,y,z) =
  //     A * exp(-(1/2)(x/a_x)^2) * exp(-(1/2)(y/a_y)^2) * exp(-(1/2)(z/a_z)^2)
  //   h_b(x,y,z) =
  //     B * exp(-(1/2)(x/b_x)^2) * exp(-(1/2)(y/b_y)^2) * exp(-(1/2)(z/b_z)^2)

  float width_a[3];    //"a" parameter in formula above
  float width_b[3];    //"b" parameter in formula above
  float m_exp;         //"m" parameter in formula above
  float n_exp;         //"n" parameter in formula above

  // How wide to make the filter?
  //Commenting out:The user can no longer set filter_truncate_halfwidth directly
  //float filter_truncate_halfwidth[3];  Instead they can set these parameters:
  float filter_truncate_ratio;     // ignore voxels further away than this*width
  float filter_truncate_threshold; // ignore voxels if filter falls below this


  // --- parameters for intensity maps and thresholding ---
  //
  // The following parameters are used whenever the user wants to
  // use an extremely simple filter that computes the brightness of each
  // voxel by only considering its own brightness
  // (without considering the brightness of any nearby voxels).
  // Examples of these filters include thresholding filters, clipping filters,
  // and inversions (reversing the light and dark voxels).
  // (As of 2019-2-23, these filters can happen regardless
  //  of the eFilterType setting selected by the user.)

  bool use_intensity_map;
  bool rescale01_in;
  bool rescale01_out;
  bool invert_output;

  bool  use_dual_thresholds;
  float out_threshold_01_a;
  float out_threshold_01_b;
  float out_threshold_10_a;
  float out_threshold_10_b;
  bool  out_thresh2_use_clipping;
  bool  out_thresh2_use_clipping_sigma;
  float out_thresh_a_value;
  float out_thresh_b_value;
  bool  use_rescale_multiply;
  float out_rescale_multiply;
  float out_rescale_offset;
  bool  use_gauss_thresholds;
  float out_thresh_gauss_x0;
  float out_thresh_gauss_sigma;

  // ---- parameters for watershed segmentation ----

  bool watershed_use_minima;
  bool watershed_show_boundaries;
  float watershed_boundary_label;
  float watershed_threshold;

  // ---- parameters for the clustering of connected voxels into islands ----

  bool cluster_connected_voxels;
  float connect_threshold_saliency;
  float connect_threshold_vector_saliency;
  float connect_threshold_vector_neighbor;
  float connect_threshold_tensor_saliency;
  float connect_threshold_tensor_neighbor;
  size_t select_cluster;
  string must_link_filename;
  vector<vector<array<int, 3> > > must_link_constraints;

  // ---- parameters used by the ridge detector used for detecting surfaces ----

  string out_normals_fname;
  float planar_hessian_score_threshold;
  bool  planar_hessian_score_threshold_is_a_fraction;
  float planar_tv_score_threshold;
  float planar_tv_sigma;
  int   planar_tv_exponent = 4;
  int   planar_tv_num_iters = 1;
  float planar_tv_truncate_ratio;
 

  // ---- parameters for scale free blob detection ----

  float dogsf_width[3];       // width of the "DOG_SCALE_FREE" filter in x,y,z directions
  float delta_sigma_over_sigma;  // Approximation parameter used in "DOG_SCALE_FREE" filter
                                 // This is the square root of "delta_t/t". See:
                                 // https://en.wikipedia.org/wiki/Blob_detection#The_Laplacian_of_Gaussian
                                 // for details.


  vector<float> blob_diameters;    // blob widths considered for scale free blob detection
  float blob_width_multiplier;
  string blob_minima_file_name;
  string blob_maxima_file_name;
  #ifndef DISABLE_INTENSITY_PROFILES
  string blob_profiles_file_name_base;
  BlobCenterCriteria blob_profiles_center_criteria;
  #endif

  // ---- parameters for performing non-max suppression on detected blobs ----

  float nonmax_min_radial_separation_ratio;
  float nonmax_max_volume_overlap_small;
  float nonmax_max_volume_overlap_large;

  float score_upper_bound;   //Ignore blobs with scores above this (usually <0)
  float score_lower_bound;   //Ignore blobs with scores less than this (>0)
  float score_bounds_are_ratios; //Are these ratios relative to min & max score?
  float sphere_diameters_lower_bound;//Ignore blobs outside a desired range of sizes
  float sphere_diameters_upper_bound;

  // ---- parameters for detecting local minima and local maxima ----

  bool find_minima = false;
  bool find_maxima = false;
  string find_minima_file_name;
  string find_maxima_file_name;
  int neighbor_connectivity = 3;
  bool extrema_on_boundary = true;

  // ---- parameters for creating images showing detected blobs ----

  string in_coords_file_name;
  string out_coords_file_name;
  float sphere_decals_diameter;
  float sphere_decals_scale;
  float sphere_decals_foreground;      //intensity of voxels in each sphere
  float sphere_decals_background;      //intensity of voxels not in any spheres
  float sphere_decals_background_scale; //multiply background voxel intensities
  bool sphere_decals_foreground_use_score; //color each sphere by its "score"?
                                       //This overrides sphere_decals_foreground
  bool sphere_decals_foreground_norm;  //Divide foreground intensity by the
                                       //number of voxels in each sphere
  float sphere_decals_shell_thickness;
  float sphere_decals_shell_thickness_min;
  bool sphere_decals_shell_thickness_is_ratio;

  // ---- parameters for depreciated features (which might go away soon) ----

  #ifndef DISABLE_TEMPLATE_MATCHING
  float template_background_radius[3];
  float template_background_exponent;
  //float template_compare_exponent;
  #endif //#ifndef DISABLE_TEMPLATE_MATCHING

  #ifndef DISABLE_BOOTSTRAPPING
  int bs_ntests;            //number of bootstrap tests to perform
  float bs_threshold;       //only consider voxels whose filtered intensity
  float bs_threshold_sign;  //(multiplied by bs_threshold_sign) exceeds this
  int   bs_scramble_radius[3];  //scramble the image range (in voxels)
  float f_bs_scramble_radius[3];//scramble the image range (physical dist)
  int bs_random_seed;  //seed for random number generator for bootstrapping
  #endif //#ifndef DISABLE_BOOTSTRAPPING



  /// @brief
  /// Determine all of these settings from the list of command-line
  /// arguments supplied by the user:

  void ParseArgs(int argc, char **argv);

 private:
  void ParseArgs(vector<string>& vArgs);
  void ConvertArgvToVectorOfStrings(int argc, 
                                    char **argv, 
                                    vector<string>& dest);

}; // class Settings



//The "EndsWith()" function checks if the last few characters of a string
//match another string.  (This is useful for manipulating file names.)
bool EndsWith(string const &fullString, string const &ending);



#endif //#ifndef _SETTINGS_HPP
