#ifndef _FILTER_3D_VARIANTS_HPP
#define _FILTER_3D_VARIANTS_HPP

/// THE CODE IN THIS FILE WAS NOT INTENDED FOR PUBLIC USE
/// This file contains several variants of functions found in <filter3d.hpp>.
/// These versions differ only in that they accept slightly different format
/// of arguments, usually providing multiple (often redundant) ways to
/// specify parameters (such as the filter-window-width).
/// The goal was to give the end user multiple different convenient options for
/// specifying these parameters, along with some reasonable default values.
/// (I don't expect end users to ever want to mess with these defaults.)
/// Unfortunately, these new functions have uglier, clumsier argument lists
/// (which is why I don't include them in <filter3d.hpp>).


#include <cstring>
#include <ostream>
#include <vector>
#include <tuple>
#include <set>
#include <cassert>
#include <queue>
#include <cmath>
#include <limits>
using namespace std;
#include <err_report.hpp>
#include <alloc3d.hpp>
#include <filter1d.hpp>  // defines "Filter1D" (used in "ApplyGauss()")
#include <filter3d.hpp> // defines "AverageArr()" and similar functions
#include <eigen3_simple.hpp>


/// @brief This is a version of the function of the same name defined in 
/// "filter2d.h".  This version allows the user to control the width of the 
/// interval considered for the filter in 2 ways: "filter_truncate_ratio", 
/// which specifies the width of the interval in units of the "width[]"
/// (ie. σ_x, σ_y) parameters, OR the "filter_truncate_threshold" which 
/// specifies how much the filter must decay before the filter is cutoff.
/// This function was not intended for public use.

template<class Scalar>
Filter2D<Scalar, int>
GenFilterGenGauss2D(Scalar width[2],    //!< "s_x", "s_y" parameters
                    Scalar m_exp,       //!< "m" parameter in formula
                    Scalar filter_truncate_ratio, //!< controls window width
                    Scalar filter_truncate_threshold, //!< controls window width
                    Scalar *pA=NULL,    //!< optional:report A coeff to user
                    ostream *pReportEquation = NULL //!< optional: print equation to the terminal
                    )
                    
{
  // choose the filter window width based on the filter_truncate_threshold
  int halfwidth[2];
  int ix = 0;

  if (filter_truncate_ratio < 0.0)
    // Choose the filter domain window based on the "truncate_threshold"
    //    filter_truncate_threshold = exp(-(filter_truncate_ratio)^m_exp);
    //    -> filter__truncate_ratio^m_exp = -log(filter_truncate_threshold)
    filter_truncate_ratio = pow(-log(filter_truncate_threshold), 1.0/m_exp);

  return GenFilterGenGauss2D(width,
                             m_exp,
                             filter_truncate_ratio,
                             pA,
                             pReportEquation);

} //GenFilterGenGauss2D(...,filter_truncate_ratio, filter_truncate_thresh,...)




/// @brief This is a version of the function of the same name defined in 
/// "filter3d.hpp".  This version allows the user to control the width of the 
/// interval considered for the filter in 2 ways: "filter_truncate_ratio", 
/// which specifies the width of the interval in units of the "width[]"
/// (ie. σ_x, σ_y, σ_z) parameters, OR the "filter_truncate_threshold" which 
/// specifies how much the filter must decay before the filter is cutoff.
/// This function was not intended for public use.

template<class Scalar>
Filter3D<Scalar, int>
GenFilterGenGauss3D(Scalar width[3],    //!<"σ_x", "σ_y", "σ_z", parameters
                    Scalar m_exp,       //!<"m" parameter in formula
                    Scalar filter_truncate_ratio, //!<controls window width
                    Scalar filter_truncate_threshold, //!<controls window width
                    Scalar *pA=NULL,    //!<optional:report A coeff to user
                    ostream *pReportEquation = NULL  //!< optional: print equation to the terminal
                    )
{
  // choose the filter window width based on the filter_truncate_threshold
  int halfwidth[3];
  int ix = 0;

  if (filter_truncate_ratio < 0.0)
    // Choose the filter domain window based on the "truncate_threshold"
    //    filter_truncate_threshold = exp(-(filter_truncate_ratio)^m_exp);
    //    -> filter__truncate_ratio^m_exp = -log(filter_truncate_threshold)
    filter_truncate_ratio = pow(-log(filter_truncate_threshold), 1.0/m_exp);

  return GenFilterGenGauss3D(width,
                             m_exp,
                             filter_truncate_ratio,
                             pA,
                             pReportEquation);
} //GenFilterGenGauss3D(...,filter_truncate_ratio, filter_truncate_thresh,...)




/// @brief This is a version of the function of the same name defined in 
/// "filter2d.h".  This version allows the user to control the width of the 
/// interval considered for the filter in 2 ways: "filter_truncate_ratio", 
/// which specifies the width of the interval in units of the "width_a/b[]"
/// (ie. a_x, b_x, ...) parameters, OR the "filter_truncate_threshold" which 
/// specifies how much the filter must decay before the filter is cutoff.
/// This function was not intended for public use.

template<class Scalar>
Filter2D<Scalar, int> 
GenFilterDogg2D(Scalar width_a[2],  //!< "a" parameter in formula
                Scalar width_b[2],  //!< "b" parameter in formula
                Scalar m_exp,       //!< "m" parameter in formula
                Scalar n_exp,       //!< "n" parameter in formula
                Scalar filter_truncate_ratio,//how many sigma before cutoff?
                Scalar filter_truncate_threshold,//cutoff when decay below threshold
                Scalar *pA=NULL, //!< optional:report A,B coeffs to user
                Scalar *pB=NULL, //!< optional:report A,B coeffs to user
                ostream *pReportProgress = NULL  //!< optional: print equation to the terminal
                )
{
  Scalar filter_truncate_ratio_A = filter_truncate_ratio;
  Scalar filter_truncate_ratio_B = filter_truncate_ratio;
  if (filter_truncate_ratio < 0.0) {
    // Choose the filter domain window based on the "filter_truncate_threshold"
    //    filter_truncate_threshold = exp(-filter_truncate_ratio^m_exp);
    //    -> filter_truncate_ratio^m_exp = -log(filter_truncate_threshold)
    filter_truncate_ratio_A = pow(-log(filter_truncate_threshold), 1.0/m_exp);
    // Do the same thing for the other Gaussian:
    filter_truncate_ratio_B = pow(-log(filter_truncate_threshold), 1.0/n_exp);
  }
  Filter2D<Scalar, int> filterXY_A =
    GenFilterGenGauss2D(width_a, //"a_x", "a_y" gaussian width parameters
                        m_exp,   //"m" exponent parameter
                        filter_truncate_ratio_A);

  Filter2D<Scalar, int> filterXY_B =
    GenFilterGenGauss2D(width_b, //"b_x", "b_y" gaussian width parameters
                        n_exp,   //"n" exponent parameter
                        filter_truncate_ratio_B);

  return _GenFilterDogg2D(width_a,
                          width_b,
                          m_exp,
                          n_exp,
                          filterXY_A, filterXY_B,
                          pA,
                          pB,
                          pReportProgress);
} //GenFilterDogg2D(...,filter_truncate_ratio, filter_truncate_threshold,...)




/// @brief This is a version of the function of the same name defined in 
/// "filter3d.hpp".  This version allows the user to control the width of the 
/// interval considered for the filter in 2 ways: "filter_truncate_ratio", 
/// which specifies the width of the interval in units of the "width_a/b[]"
/// (ie. a_x, b_x, ...) parameters, OR the "filter_truncate_threshold" which 
/// specifies how much the filter must decay before the filter is cutoff.
/// This function was not intended for public use.

template<class Scalar>
// Create a 3-D filter and fill it with a difference of (generalized) Gaussians:
Filter3D<Scalar, int> 
GenFilterDogg3D(Scalar width_a[3],  //"a" parameter in formula
                Scalar width_b[3],  //"b" parameter in formula
                Scalar m_exp,       //"m" parameter in formula
                Scalar n_exp,       //"n" parameter in formula
                Scalar filter_truncate_ratio,//how many sigma before cutoff?
                Scalar filter_truncate_threshold,//cutoff when decay below threshold
                Scalar *pA = NULL, //optional:report A,B coeffs to user
                Scalar *pB = NULL, //optional:report A,B coeffs to user
                ostream *pReportProgress = NULL  //!< optional: print equation to the terminal
                )
{
  Scalar filter_truncate_ratio_A = filter_truncate_ratio;
  Scalar filter_truncate_ratio_B = filter_truncate_ratio;
  if (filter_truncate_ratio < 0.0) {
    // Choose the filter domain window based on the "filter_truncate_threshold"
    //    filter_truncate_threshold = exp(-filter_truncate_ratio^m_exp);
    //    -> filter_truncate_ratio^m_exp = -log(filter_truncate_threshold)
    filter_truncate_ratio_A = pow(-log(filter_truncate_threshold), 1.0/m_exp);
    // Do the same thing for the other Gaussian:
    filter_truncate_ratio_B = pow(-log(filter_truncate_threshold), 1.0/n_exp);
  }
  Filter3D<Scalar, int> filter_A =
    GenFilterGenGauss3D(width_a, //"a_x", "a_y", "a_z" gaussian width parameters
                        m_exp,   //"m" exponent parameter
                        filter_truncate_ratio_A);

  Filter3D<Scalar, int> filter_B =
    GenFilterGenGauss3D(width_b, //"b_x", "b_y", "b_z" gaussian width parameters
                        n_exp,   //"n" exponent parameter
                        filter_truncate_ratio_B);

  return _GenFilterDogg3D(width_a,
                          width_b,
                          m_exp,
                          n_exp,
                          filter_A, filter_B,
                          pA,
                          pB,
                          pReportProgress);

} //GenFilterDogg3D(...,filter_truncate_ratio, filter_truncate_threshold,...)









/// @brief This is a version of the function of the same name defined in 
/// "filter3d.hpp".  This version allows the user to control the width of the 
/// interval considered for the filter in 2 ways: "filter_truncate_ratio", 
/// which specifies the width of the interval in units of the "width[]"
/// (ie. σ_x, σ_y, σ_z) parameters, OR the "filter_truncate_threshold" which 
/// specifies how much the filter must decay before the filter is cutoff.
/// This function was not intended for public use.

template<class Scalar>
Scalar
ApplyGauss3D(int const image_size[3], 
             Scalar const *const* const *aaafSource,
             Scalar ***aaafDest,
             Scalar const *const *const *aaafMask,
             Scalar const sigma[3],
             Scalar filter_truncate_ratio,
             Scalar filter_truncate_threshold,
             bool normalize = true,
             ostream *pReportProgress = NULL)
{

  if (filter_truncate_ratio <= 0) {
    assert(filter_truncate_threshold > 0.0);
    //    filter_truncate_threshold = exp(-(1/2)*filter_truncate_ratio^2);
    //    -> filter_truncate_ratio^2 = -2*log(filter_truncate_threshold)
    filter_truncate_ratio = sqrt(-2*log(filter_truncate_threshold));
  }
  return ApplyGauss3D(image_size, 
                      aaafSource,
                      aaafDest,
                      aaafMask,
                      sigma,
                      filter_truncate_ratio,
                      normalize,
                      pReportProgress);

} // ApplyGauss3D(..., filter_truncate_ratio, filter_truncate_threshold, ...)





/// @brief This is a version of the function of the same name defined in 
/// "filter3d.hpp".  This version allows the user to control the width of the 
/// interval considered for the filter in 2 ways: "filter_truncate_ratio", 
/// which specifies the width of the interval in units of the "sigma_a/b[]"
/// (ie. a_x, b_x, ...) parameters, OR the "filter_truncate_threshold" which 
/// specifies how much the filter must decay before the filter is cutoff.
/// This function was not intended for public use.

template<class Scalar>
void
ApplyDog3D(int const image_size[3], //source image size
           Scalar const *const *const *aaafSource,  //source image
           Scalar ***aaafDest,     //save results here
           Scalar const *const *const *aaafMask,  //ignore voxels where mask==0
           Scalar sigma_a[3],
           Scalar sigma_b[3],
           Scalar filter_truncate_ratio,
           Scalar filter_truncate_threshold,
           Scalar *pA = NULL,
           Scalar *pB = NULL,
           ostream *pReportProgress = NULL)
{
  Scalar ***aaafTemp; //temporary array to store partially processed tomogram
  Scalar *afTemp;     //temporary array to store partially processed tomogram

  Alloc3D(image_size,
          &afTemp,
          &aaafTemp);

  Scalar A, B;        // let the user know what A B coefficients were used

  A = ApplyGauss3D(image_size,
                   aaafSource,
                   aaafDest,         // <-- save result here
                   aaafMask,
                   sigma_a,
                   filter_truncate_ratio,
                   filter_truncate_threshold,
                   true,
                   pReportProgress);

  B = ApplyGauss3D(image_size,
                   aaafSource,
                   aaafTemp,         // <-- save result here
                   aaafMask,
                   sigma_b,
                   filter_truncate_ratio,
                   filter_truncate_threshold,
                   true,
                   pReportProgress);

  // Subtract the second convolved signal from the first
  for (int iz = 0; iz < image_size[2]; iz++)
    for (int iy = 0; iy < image_size[1]; iy++)
      for (int ix = 0; ix < image_size[0]; ix++)
        aaafDest[iz][iy][ix] -= aaafTemp[iz][iy][ix];

  // Deallocate the temporary array
  Dealloc3D(image_size,
            &afTemp,
            &aaafTemp);

  // Report the A and B normalization coefficients to the caller?
  if (pA)
    *pA = A;
  if (pB)
    *pB = B;

} // ApplyDog3D(...,filter_truncate_ratio,filter_truncate_threshold,...)





template<class Scalar>
void
ApplyLog3D(int const image_size[3], //source image size
           Scalar const *const *const *aaafSource,   //source image
           Scalar ***aaafDest,     //save results here
           Scalar const *const *const *aaafMask,     //ignore voxels where mask==0
           Scalar const sigma[3],  //Gaussian width in x,y,z drections
           Scalar delta_sigma_over_sigma, //difference in Gauss widths
           Scalar filter_truncate_ratio,
           Scalar filter_truncate_threshold,
           Scalar *pA = NULL,
           Scalar *pB = NULL,
           ostream *pReportProgress = NULL)
{

  if (filter_truncate_ratio < 0)
    //    filter_truncate_threshold = exp(-(1/2)*filter_truncate_ratio^2);
    //    -> filter_truncate_ratio^2 = -2*log(filter_truncate_threshold)
    filter_truncate_ratio = sqrt(-2*log(filter_truncate_threshold));

  ApplyLog3D(image_size,
             aaafSource,
             aaafDest,
             aaafMask,
             sigma,
             delta_sigma_over_sigma,
             filter_truncate_ratio,
             pA,
             pB,
             pReportProgress);

} // ApplyLog3D(..,filter_truncate_ratio,filter_truncate_threshold,...)




/// @brief Calculate the fluctuations in intensity around every voxel lying in
///        a Gaussian-weighted ellipsoidal volume with half-width along x,y,z
///        given by radius[].  (You can approximate a uniform sphere by setting
///        all entries in radius[] equal to eachother, and by setting
///        template_background_exponent argument to a large number.)
///        This version of the function gives the user the option to specify
///        the width of the truncation filter window using a decay threshold
///        (OR alternately by specifying a distance ratio).
///        For details regarding what this function does, see the comments
///        for the other version of this function (located in filter3d.hpp).

template<class Scalar, class Integer>
static void
LocalFluctuationsRadial(Integer const image_size[3], //!< number of voxels in x,y,z directions
                        Scalar const *const *const *aaafSource, //!< original image
                        Scalar ***aaafDest, //!< store filtered image here (fluctuation magnitude)
                        Scalar const *const *const *aaafMask, //!< optional: if not NULL then ignore voxel ix,iy,iz if aaafMask[iz][iy][ix]==0
                        Scalar radius[3],  //!< radius (=sigma/√3) of neighborhooed over which we search in x,y,z directions (ellipsoidal shaped search area)
                        Scalar template_background_exponent=2, //!< exponent controlling sharpness of the (generalized) Gaussian (slow if != 2)
                        Scalar filter_truncate_ratio=2.5, //!< width over which we search is this many times larger than the gaussian width parameter (sigma)
                        Scalar filter_truncate_threshold=0.02, //!< alternatively, specify how much the Gaussian must decay before we truncate it
                        bool normalize = true, //!< normalize the result?
                        ostream *pReportProgress = NULL //!< report progress to the user?
                  )
{
  if (filter_truncate_ratio < 0.0)
    // Choose the filter domain window based on the "truncate_threshold"
    //    filter_truncate_threshold = exp(-(filter_truncate_ratio)^exponent);
    //    -> filter__truncate_ratio^exponent = -log(filter_truncate_threshold)
    filter_truncate_ratio = pow(-log(filter_truncate_threshold),
                                1.0 / template_background_exponent);

  LocalFluctuationsRadial(image_size,
                          aaafSource,
                          aaafDest,
                          aaafMask,
                          radius,
                          template_background_exponent,
                          filter_truncate_ratio,
                          normalize,
                          pReportProgress);

} //LocalFluctuations()



#endif //#ifndef _FILTER_3D_VARIANTS_HPP
