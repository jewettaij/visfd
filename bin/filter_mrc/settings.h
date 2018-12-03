#ifndef _SETTINGS_H
#define _SETTINGS_H

#define DISABLE_BOOTSTRAPPING
//#define DISABLE_TEMPLATE_MATCHING


#include<vector>
#include<string>
using namespace std;

// Generally speaking, the filter used is of the form of a
// "difference of gaussians":
//
//   h(x,y,z) = h_a(x,y,z) - h_b(x,y,z)
//
// Both h_a() and h_b() are generalized normal distributions (Gaussians)
// https://en.wikipedia.org/wiki/Generalized_normal_distribution
//
// If b_x>0, b_y>0, b_z>0, then:
//
//   h_a(x,y,z) = A*exp(-r_a^m)
//          r_a = sqrt((x/a_x)*2 + (y/a_y)*2 + (z/a_z)*2)
//   h_b(x,y,z) = B*exp(-r_b^n)
//          r_b = sqrt((x/b_x)*2 + (y/b_y)*2 + (z/b_z)*2)
//
// The "A" and "B" parameters are chosen so that sum_{x,y,z} h(x,y,z) = 0
// (In the special case below where h_b=0, then sum_{x,y,z} h(x,y,z) = 1)
// You may wish to use a "difference of gaussians" approach in some
// directions, and ordinary gaussians" in other directions.
//
// If b_z < 0, then both h_a and h_b become "seperable" filters, IE
//             functions of x,y, multiplied by a Gaussian in the Z direction:
//   h_a(x,y,z) = A * exp(-r_a^m)   *  C * exp(-(1/2)z^2/a_z)
//          r_a = sqrt((x/a_x)*2 + (y/a_y)*2
//   h_b(x,y,z) = B * exp(-r_b^n)   *  C * exp(-(1/2)z^2/a_z)
//          r_b = sqrt((x/b_x)*2 + (y/b_y)*2
// Seperable filters are faster to compute.
// (In this case, "C" is chosen so that C * sum_z exp(-(1/2)z^2/a_z) = 1)
//
// Likewise for each dimension where the corresponding "b" parameter is < 0,
// we use a "separable" filter which is Gaussian in that dimension.
// The second "B" term is only calculated in dimensions for which t>0.
// In this way, the function above is general.
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


class Settings {
 public:

  float voxel_width;  //physical width of each voxel (for example in Angstroms)
  bool voxel_width_divide_by_10; //Use nm instead of Angstroms?

  string in_file_name;
  string out_file_name;
  string mask_file_name;
  int mask_select;
  bool use_mask_select;
  int mask_out;
  bool use_mask_out;

  typedef enum eFilterType {
    NONE,    //(the user has not selected a filter)
    GAUSS,   //3D Gaussian filter
    GGAUSS,  //3D Generalized Gaussian filter
    BLOB,    //3D Blob detection of arbitrary size blobs
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
    LOCAL_FLUCTUATIONS, // Report the fluctuation of nearby voxel intensities
    BOOTSTRAP_DOGG,  // DOGG filter with bootstrapping (significance testing)
    MIN_DISTANCE,    //voxel intensity = distance to nearest object
    SPHERE_DECALS,   //voxel intensity = 1 if within R from the nearest object
    SPHERE_NONMAX_SUPPRESSION,
    RIDGE_PLANAR
  } FilterType; 
  
  FilterType filter_type;

  float width_a[3];    //"a" parameter in formula above
  float width_b[3];    //"b" parameter in formula above
  float m_exp;         //"m" parameter in formula above
  float n_exp;         //"n" parameter in formula above

  float dogsf_width[3];       // width of the "DOG_SCALE_FREE" filter in x,y,z directions
  float delta_sigma_over_sigma;  // Approximation parameter used in "DOG_SCALE_FREE" filter
                                 // This is the square root of "delta_t/t". See:
                                 // https://en.wikipedia.org/wiki/Blob_detection#The_Laplacian_of_Gaussian
                                 // for details.


  vector<float> blob_diameters;    // blob widths considered for scale free blob detection
  float blob_width_multiplier;
  string blob_minima_file_name;
  string blob_maxima_file_name;
  float nonmax_min_radial_separation_ratio;
  float nonmax_max_volume_overlap_small;
  float nonmax_max_volume_overlap_large;

  float score_upper_bound;   //Ignore blobs with scores above this (usually <0)
  float score_lower_bound;   //Ignore blobs with scores less than this (>0)
  float score_bounds_are_ratios; //Are these ratios relative to min & max score?
  float sphere_diameters_lower_bound;//Ignore blobs outside a desired range of sizes
  float sphere_diameters_upper_bound;

  float filter_truncate_threshold;
  float filter_truncate_ratio;
  //Commenting out:The user can no longer set filter_truncate_halfwidth directly
  //float filter_truncate_halfwidth[3];

  bool rescale01_in;
  bool rescale01_out;
  bool invert_output;

  bool use_thresholds;
  bool use_dual_thresholds;
  float out_threshold_01_a;
  float out_threshold_01_b;
  float out_threshold_10_a;
  float out_threshold_10_b;
  bool  out_thresh2_use_clipping;
  bool  out_thresh2_use_clipping_sigma;
  float out_thresh_a_value;
  float out_thresh_b_value;

  #ifndef DISABLE_BOOTSTRAPPING
  int bs_ntests;            //number of bootstrap tests to perform
  float bs_threshold;       //only consider voxels whose filtered intensity
  float bs_threshold_sign;  //(multiplied by bs_threshold_sign) exceeds this
  int   bs_scramble_radius[3];  //scramble the image range (in voxels)
  float f_bs_scramble_radius[3];//scramble the image range (physical dist)
  int bs_random_seed;  //seed for random number generator for bootstrapping
  #endif //#ifndef DISABLE_BOOTSTRAPPING

  #ifndef DISABLE_TEMPLATE_MATCHING
  float template_background_radius[3];
  float template_background_exponent;
  //float template_compare_exponent;
  #endif //#ifndef DISABLE_TEMPLATE_MATCHING

  bool find_minima = false;
  bool find_maxima = false;
  string find_minima_file_name;
  string find_maxima_file_name;
  //REMOVE THIS CRUFT:
  //float find_extrema_occlusion_ratio; // nonmax_min_radial_separation_ratio

  string in_coords_file_name;
  string out_coords_file_name;
  float sphere_decals_diameter;
  float sphere_decals_scale;
  float sphere_decals_foreground;      //intensity of voxels in each sphere
  float sphere_decals_background;      //intensity of voxels not in any spheres
  float sphere_decals_background_scale; //multiply background voxel intensities
  bool sphere_decals_foreground_use_score; //color each sphere by its "score"?
                                       //This overrides sphere_decals_foreground
  //bool sphere_decals_background_use_orig;  //superimpose spheres on background image
                                       //This overrides sphere_decals_background
  bool sphere_decals_foreground_norm;  //Divide foreground intensity by the
                                       //number of voxels in each sphere
  float sphere_decals_shell_thickness;
  float sphere_decals_shell_thickness_min;
  bool sphere_decals_shell_thickness_is_ratio;

  string out_normals_fname;
  float planar_hessian_score_threshold;
  float planar_tv_score_threshold;
  float planar_tv_sigma;
  int   planar_tv_exponent = 4;
  int   planar_tv_num_iters = 1;
 
  // not used (yet):
  //float missing_wedge_min[2];  //range of angles sampled (around x and y axis)
  //float missing_wedge_max[2];  //range of angles sampled (around x and y axis)

  Settings();
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



#endif //#ifndef _SETTINGS_H
