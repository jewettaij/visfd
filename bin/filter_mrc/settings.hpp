#ifndef _SETTINGS_HPP
#define _SETTINGS_HPP

#define DISABLE_BOOTSTRAPPING
//#define DISABLE_TEMPLATE_MATCHING


#include<vector>
#include<string>
using namespace std;
#include <mrc_simple.hpp>
#include <visfd.hpp>
using namespace visfd;
#include "feature_unsupported.hpp"


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
    NONE,     //(The user has not selected a filter, or is using thresholding)
    FIND_EXTREMA, // find local minima or maxima in the image
    MEDIAN,   // a median filter (https://en.wikipedia.org/wiki/Median_filter)
    DILATION, // grayscale image dilation
    EROSION,  // grayscale image erosion
    OPENING,  // grayscale image opening
    CLOSING,  // grayscale image closing
    TOP_HAT_WHITE, // top hat transforms
    TOP_HAT_BLACK, // (see https://en.wikipedia.org/wiki/Top-hat_transform)
    GAUSS,   //3D Gaussian filter
    GGAUSS,  //3D Generalized Gaussian filter (slow)
    DOG,     //3D Difference of Gaussians (arbitrary widths)
    DOGG,    //3D Generalized Difference-of-Generalized-Gaussians with arbitrary
             //   widths and arbitrary exponents, slow)
    LOG_DOG,         // Approximation to Laplacian-of-Guassian using DOG (fast)
    LOCAL_FLUCTUATIONS, // Report the fluctuation of nearby voxel intensities
    WATERSHED,       // Watershed segmentation
    LABEL_CONNECTED, // Connected component labeling (count islands in a sea)
    CURVE,           // Detect 1D curve-like ridges (a filament detector)
    SURFACE_EDGE,    // Detect the edge of a light-dark boundary (gradient)
    SURFACE_RIDGE,   // Detect 2D surface-like ridges (a thin membrane detector)
    BLOB,            // Scale-free blob detection
    //   Analyze a list of objects supplied by the user
    BLOB_NONMAX_SUPPRESSION,//throw away overlapping objects(detected earlier)
    BLOB_NONMAX_SUPERVISED_MULTI,//use training data to throw away bad scoring blobs in multiple different images
    //   Annotate the image by drawing simple objects in the image
    DRAW_SPHERES,    // draw spheres onto the image

    //   DEPRECIATED:
    DOGGXY,  //2D Generalized Difference-of-Generalized-Gaussians with arbitrary
             //   widths and arbitrary exponents in XY, and Gaussian blur in Z.
    DISTANCE_TO_POINTS, //voxel intensity=distance to nearest point in a set
    DISTANCE_TO_VOXELS, //report distance of each point to selected voxels
    RANDOM_SPHERES,  // generate a list of randomly positioned spheres
    TEMPLATE_GAUSS,  // Perform a least-squares fit between a Gaussian
                     // and the surrounding voxel brightnesses
    TEMPLATE_GGAUSS, // Perform a least-squares fit between a generalized
                     // Gaussian and the surrounding voxel brightnesses
    BOOTSTRAP_DOGG,  // DOGG filter with bootstrapping (significance testing)
    BLOB_RADIAL_INTENSITY // brightness-vs-distance-to-center for each blob
  } FilterType; 
  
  FilterType filter_type; // The user must select one of these filter types


  // ---- parameters describing the physical size of each voxel ----
  //      (shared by all filter types)

  float voxel_width;  //physical width of each voxel (for example in Angstroms)
  bool voxel_width_divide_by_10; //Use nm instead of Angstroms?
  int image_size_orig[3]; //number of voxels in the source image before binning
  float cellA_orig[3]; //backup of the voxel width in the source image.


  // -- parameters which determine the image files we will be reading/writing --
  //      (shared by all filter types)
  
  string in_file_name; // name of the image file we want to read
  int in_set_image_size[3];//image size (if the user did not supply a file name)
  string out_file_name; // name of the image file we want to create
  bool out_file_overwrite; // allow out_file_name to equal in_file_name?
  // Mask parameters are used to select (ignore) voxels from the original image.
  bool normalize_near_boundaries; // compensate when mask boundary cuts filter?
  string mask_file_name; // name of an image file used for masking
  bool use_mask_select; // do we select voxels with a specific value?
  int mask_select; // select only voxels from the input with this value
  bool specify_masked_brightness; //did the user specify the brightness of masked voxels?
  float masked_voxel_brightness;//what brightness do we assign to masked voxels?
  float undefined_voxel_brightness; //what brightness do we assign to other kinds of undefined voxels?
  bool undefined_voxels_are_max; // set undefined voxels to the maximum image brightness + 1
  bool is_mask_crds_in_voxels; // Are mask crds in units of voxels or physical distance?
  vector<SimpleRegion<float> > vMaskRegions; // set mask=one or more regions in space
  int resize_with_binning; //width of bin to use to reduce image size?
  bool resize_with_binning_explicit; //did the user explictly ask for binning?

  string save_intermediate_fname_base = ""; //save progress of slow calculations
  string load_intermediate_fname_base = ""; //load progress from a previous run?

  float median_radius;       // radius used for the median filter
  float morphology_r;        // distance parameter for morphological operations
  float morphology_rmax;     // distance parameter for morphological operations
  float morphology_bmax;     // maximum b value used in soft structure factors

  // ---- parameters for "difference of generalized gaussians" filters ----
  //
  // A significant number of the filters used by this program
  // are all of the type "difference of generalized gaussians".
  // (Specifically, these are: GAUSS, GGAUSS, DOG, DOGG, DOGGXY, LOG_DOG,
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
  bool rescale_min_max_in;
  bool rescale_min_max_out;
  float in_rescale_min;
  float in_rescale_max;
  float out_rescale_min;
  float out_rescale_max;
  bool rescale01_in;
  bool rescale01_out;
  bool invert_output;

  bool  use_dual_thresholds;
  float in_threshold_01_a;
  float in_threshold_01_b;
  float in_threshold_10_a;
  float in_threshold_10_b;
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

  bool watershed_show_boundaries;
  float watershed_boundary_label;
  float watershed_threshold;
  string watershed_markers_filename;

  // ---- parameters for connected component analysis ----

  bool clusters_begin_at_maxima;
  bool cluster_connected_voxels;//NOTE TO SELF:I should eliminate or rename this
  float connect_threshold_saliency; // I should eliminate or rename this
  size_t select_cluster;            // I should rename this variable
  string must_link_filename;
  vector<vector<array<float, 3> > > must_link_constraints;
  vector<vector<DirectionPairType> > must_link_constraint_directions;
  bool is_must_link_constraints_in_voxels = false;

  // ---- parameters for directional connected component analysis ----

  float connect_threshold_vector_saliency;
  float connect_threshold_vector_neighbor;
  float connect_threshold_tensor_saliency;
  float connect_threshold_tensor_neighbor;

  // ---- parameters used by the curve and surface detectors ----

  string out_normals_fname;
  bool  ridges_are_maxima;
  float max_distance_to_feature;
  float surface_normal_curve_ds;
  bool surface_find_ridge;
  float hessian_score_threshold;
  bool  hessian_score_threshold_is_a_fraction;
  float tv_score_threshold;
  float tv_sigma;
  int   tv_exponent;
  //int   tv_num_iters = 1;  <-- CRUFT.  REMOVE THIS LINE LATER
  float tv_truncate_ratio;

  // ---- parameters for scale free blob detection ----

  float log_width[3];            // width of the "LOG_DOG" filter in x,y,z directions
  float delta_sigma_over_sigma;  // Approximation parameter used in "LOG_DOG" filter
                                 // This is the square root of "delta_t/t". See:
                                 // https://en.wikipedia.org/wiki/Blob_detection#The_Laplacian_of_Gaussian
                                 // for details.


  vector<float> blob_diameters;    // blob widths considered for scale free blob detection
  float blob_width_multiplier;
  string blob_minima_file_name;
  string blob_maxima_file_name;
  float blob_aspect_ratio[3];
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

  bool auto_thresh_score = false;
  vector<array<float, 3> > training_pos_crds;
  vector<array<float, 3> > training_neg_crds;

  // training data for a single images
  bool is_training_pos_in_voxels = false;
  bool is_training_neg_in_voxels = false;

  // generate a random point cloud?
  int    rand_crds_n;
  double rand_crds_diameter;
  int    rand_crds_seed;

  // measure distances in image
  string out_distances_file_name;

  // training data in multiple independent files (for multiple images)
  vector<string> multi_in_crds_file_names;
  vector<bool> multi_is_training_pos_in_voxels;
  vector<bool> multi_is_training_neg_in_voxels;
  vector<vector<array<float, 3> > > multi_training_pos_crds;
  vector<vector<array<float, 3> > > multi_training_neg_crds;
  // ---- parameters for detecting local minima and local maxima ----

  bool find_minima = false;
  bool find_maxima = false;
  string find_minima_file_name;
  string find_maxima_file_name;
  int neighbor_connectivity = 3;
  bool extrema_on_boundary = true;

  // ---- parameters for creating images showing detected blobs ----

  vector<string> in_crds_file_names;
  string out_crds_file_name;
  float sphere_decals_diameter;
  bool sphere_decals_diameter_in_voxels;
  float sphere_decals_scale;
  float sphere_decals_foreground;      //intensity of voxels in each sphere
  float sphere_decals_background;      //intensity of voxels not in any spheres
  float sphere_decals_background_scale; //multiply background voxel intensities
  bool sphere_decals_foreground_use_score; //color each sphere by its "score"?
                                       //This overrides sphere_decals_foreground
  bool sphere_decals_background_norm;  //Automatically change shift and rescale
                                       //brightness of the background voxels?
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
