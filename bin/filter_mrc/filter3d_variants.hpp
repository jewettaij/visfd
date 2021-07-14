#ifndef _FILTER_3D_VARIANTS_HPP
#define _FILTER_3D_VARIANTS_HPP




namespace visfd {




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
#include <visfd.hpp>



/// @brief This is a version of the function of the same name defined in 
/// "filter2d.h".  This version allows the user to control the width of the 
/// interval considered for the filter in 2 ways: "filter_truncate_ratio", 
/// which specifies the width of the interval in units of the "width[]"
/// (ie. σ_x, σ_y) parameters, OR the "filter_truncate_threshold" which 
/// specifies how much the filter must decay before the filter is cutoff.
/// This function was not intended for public use.

template<typename Scalar>
Filter2D<Scalar, int>
GenFilterGenGauss2D(const Scalar width[2], //!< "s_x", "s_y" parameters
                    Scalar m_exp,          //!< "m" parameter in formula
                    Scalar filter_truncate_ratio, //!< controls window width
                    Scalar filter_truncate_threshold, //!< controls window width
                    Scalar *pA=nullptr,    //!< optional:report A coeff to user
                    ostream *pReportEquation = nullptr //!< optional: print equation to the terminal
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
                             pA);

} //GenFilterGenGauss2D(...,filter_truncate_ratio, filter_truncate_thresh,...)




/// @brief This is a version of the function of the same name defined in 
/// "filter3d.hpp".  This version allows the user to control the width of the 
/// interval considered for the filter in 2 ways: "filter_truncate_ratio", 
/// which specifies the width of the interval in units of the "width[]"
/// (ie. σ_x, σ_y, σ_z) parameters, OR the "filter_truncate_threshold" which 
/// specifies how much the filter must decay before the filter is cutoff.
/// This function was not intended for public use.

template<typename Scalar>
Filter3D<Scalar, int>
GenFilterGenGauss3D(const Scalar width[3], //!<"σ_x", "σ_y", "σ_z", parameters
                    Scalar m_exp,          //!<"m" parameter in formula
                    Scalar filter_truncate_ratio, //!<controls window width
                    Scalar filter_truncate_threshold, //!<controls window width
                    Scalar *pA=nullptr    //!<optional:report A coeff to user
                    )
{
  // choose the filter window width based on the filter_truncate_threshold
  int halfwidth[3];
  int ix = 0;
  float A;

  if (filter_truncate_ratio < 0.0)
    // Choose the filter domain window based on the "truncate_threshold"
    //    filter_truncate_threshold = exp(-(filter_truncate_ratio)^m_exp);
    //    -> filter__truncate_ratio^m_exp = -log(filter_truncate_threshold)
    filter_truncate_ratio = pow(-log(filter_truncate_threshold), 1.0/m_exp);

  return GenFilterGenGauss3D(width,
                             m_exp,
                             filter_truncate_ratio,
                             pA);

} //GenFilterGenGauss3D(...,filter_truncate_ratio, filter_truncate_thresh,...)




/// @brief
/// Create a 2D filter and fill it with a difference of (generalized) Gaussians:
/// This version requires that the caller has already created individual
/// filters for the two gaussians.
/// All this function does is subtract one filter from the other (and rescale).

template<typename Scalar>

static Filter2D<Scalar, int> 
_GenFilterDogg2D(const Scalar width_a[2],  //!< "a" parameter in formula
                 const Scalar width_b[2],  //!< "b" parameter in formula
                 Filter2D<Scalar, int>& filterXY_A, //!< filters for the two
                 Filter2D<Scalar, int>& filterXY_B, //!< gaussians
                 Scalar *pA=nullptr, //!< optional:report A,B coeffs to user
                 Scalar *pB=nullptr  //!< optional:report A,B coeffs to user
                 )
{

  Scalar A, B;
  //A, B = height of the central peak
  A = filterXY_A.aafH[0][0];
  B = filterXY_B.aafH[0][0];

  // The "difference of gaussians" filter is the difference between
  // these two (generalized) gaussian filters.
  int halfwidth[2];
  halfwidth[0] = std::max(filterXY_A.halfwidth[0], filterXY_B.halfwidth[0]);
  halfwidth[1] = std::max(filterXY_A.halfwidth[1], filterXY_B.halfwidth[1]);
  Filter2D<Scalar, int> filter(halfwidth);

  for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
    for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
      filter.aafH[iy][ix] = 0.0;

      // The two filters may have different widths, so we have to check
      // that ix and iy lie within the domain of these two filters before
      // adding or subtracting their values from the final GDOG filter.
      if (((-filterXY_A.halfwidth[0]<=ix) && (ix<=filterXY_A.halfwidth[0])) &&
          ((-filterXY_A.halfwidth[1]<=iy) && (iy<=filterXY_A.halfwidth[1])))

        filter.aafH[iy][ix] +=filterXY_A.aafH[iy][ix]; // /(A-B); COMMENTING OUT
                         
      // COMMENTING OUT: (The factor of 1/(A-B) insures that the central peak has height 1)

      if (((-filterXY_B.halfwidth[0]<=ix) && (ix<=filterXY_B.halfwidth[0])) &&
          ((-filterXY_B.halfwidth[1]<=iy) && (iy<=filterXY_B.halfwidth[1])))

        filter.aafH[iy][ix] -=filterXY_B.aafH[iy][ix]; // /(A-B); COMMENTING OUT
                         
    } // for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
  } // for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {

  // COMMENTING OUT 1/(A-B)
  //if (pA && pB) {
  //  *pA = A/(A-B); // Rescale A and B numbers returned to the caller
  //  *pB = B/(A-B); // (because we divided the array entries by (A-B) earlier)
  //}

  return filter;

} //_GenFilterDogg2D()



template<typename Scalar>
// Create a 2D filter and fill it with a difference of (generalized) Gaussians:
Filter2D<Scalar, int> 
GenFilterDogg2D(const Scalar width_a[2],  //"a" parameter in formula
                const Scalar width_b[2],  //"b" parameter in formula
                Scalar m_exp,             //"m" parameter in formula
                Scalar n_exp,             //"n" parameter in formula
                const int halfwidth[2],   //size of 2d array to store the filter
                Scalar *pA = nullptr,     //optional:report A,B coeffs to user
                Scalar *pB = nullptr      //optional:report A,B coeffs to user
                )
{
  Filter2D<Scalar, int> filterXY_A =
    GenFilterGenGauss2D(width_a,      //"a_x", "a_y" gaussian width parameters
                        m_exp,        //"n" exponent parameter
                        halfwidth);

  Filter2D<Scalar, int> filterXY_B =
    GenFilterGenGauss2D(width_b,      //"b_x", "b_y" gaussian width parameters
                        n_exp,        //"n" exponent parameter
                        halfwidth);

  return _GenFilterDogg2D(width_a,
                          width_b,
                          filterXY_A, filterXY_B,
                          pA,
                          pB);
} //GenFilterDogg2D(...halfwidth...)



/// @brief This is a version of the function of the same name defined in 
/// "filter2d.h".  This version allows the user to control the width of the 
/// interval considered for the filter in 2 ways: "filter_truncate_ratio", 
/// which specifies the width of the interval in units of the "width_a/b[]"
/// (ie. a_x, b_x, ...) parameters, OR the "filter_truncate_threshold" which 
/// specifies how much the filter must decay before the filter is cutoff.
/// This function was not intended for public use.

template<typename Scalar>
Filter2D<Scalar, int> 
GenFilterDogg2D(const Scalar width_a[2],  //!< "a" parameter in formula
                const Scalar width_b[2],  //!< "b" parameter in formula
                Scalar m_exp,             //!< "m" parameter in formula
                Scalar n_exp,             //!< "n" parameter in formula
                Scalar filter_truncate_ratio,//how many sigma before cutoff?
                Scalar filter_truncate_threshold,//cutoff when decay below threshold
                Scalar *pA=nullptr, //!< optional:report A,B coeffs to user
                Scalar *pB=nullptr  //!< optional:report A,B coeffs to user
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
                          filterXY_A, filterXY_B,
                          pA,
                          pB);
} //GenFilterDogg2D(...,filter_truncate_ratio, filter_truncate_threshold,...)




/// @brief  Create a 3D filter and fill it with a difference of (generalized) Gaussians 
///
/// This version requires that the caller has already created individual
/// filters for the two gaussians.
/// All this function does is subtract one filter from the other (and rescale).
/// @note: DEPRECIATION WARNING: It's not clear if this filter is useful.
///                              I may delete this function in the future.
/// @note: THIS FUNCTION WAS NOT INTENDED FOR PUBLIC USE
template<typename Scalar>

static Filter3D<Scalar, int> 
_GenFilterDogg3D(const Scalar width_a[3],  //!< "a" parameter in formula
                 const Scalar width_b[3],  //!< "b" parameter in formula
                 Scalar m_exp,             //!< "m" parameter in formula
                 Scalar n_exp,             //!< "n" parameter in formula
                 const Filter3D<Scalar, int>& filter_A, //!< filters for the two
                 const Filter3D<Scalar, int>& filter_B, //!< gaussians
                 Scalar *pA=nullptr, //!< optional:report A,B coeffs to user
                 Scalar *pB=nullptr, //!< optional:report A,B coeffs to user
                 ostream *pReportEquation = nullptr //!< optional: report equation to the user
                 )
{
  Scalar A, B;
  //A, B = height of the central peak
  A = filter_A.aaafH[0][0][0];
  B = filter_B.aaafH[0][0][0];


  // The "difference of gaussians" filter is the difference between
  // these two (generalized) gaussian filters.
  int halfwidth[3];
  halfwidth[0] = std::max(filter_A.halfwidth[0], filter_B.halfwidth[0]);
  halfwidth[1] = std::max(filter_A.halfwidth[1], filter_B.halfwidth[1]);
  halfwidth[2] = std::max(filter_A.halfwidth[2], filter_B.halfwidth[2]);
  Filter3D<Scalar, int> filter(halfwidth);

  //FOR DEBUGGING REMOVE EVENTUALLY
  //if (pReportEquation)
  //  *pReportEquation << "Array of 3D filter entries:" << endl;


  for (int iz=-halfwidth[2]; iz<=halfwidth[2]; iz++) {
    for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
      for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {

        filter.aaafH[iz][iy][ix] = 0.0;

        // The two filters may have different widths, so we have to check
        // that ix,iy and iz lie within the domain of these two filters before
        // adding or subtracting their values from the final GDoG filter.
        if (((-filter_A.halfwidth[0]<=ix) && (ix<=filter_A.halfwidth[0])) &&
            ((-filter_A.halfwidth[1]<=iy) && (iy<=filter_A.halfwidth[1])) &&
            ((-filter_A.halfwidth[2]<=iz) && (iz<=filter_A.halfwidth[2])))

          filter.aaafH[iz][iy][ix] +=
            filter_A.aaafH[iz][iy][ix];     //  /(A-B);  COMMENTING OUT
                          
                          

        // COMMENTING OUT: (The factor of 1/(A-B) insures that the central peak has height 1)


        if (((-filter_B.halfwidth[0]<=ix) && (ix<=filter_B.halfwidth[0])) &&
            ((-filter_B.halfwidth[1]<=iy) && (iy<=filter_B.halfwidth[1])) &&
            ((-filter_B.halfwidth[2]<=iz) && (iz<=filter_B.halfwidth[2])))

          filter.aaafH[iz][iy][ix] -=
            filter_B.aaafH[iz][iy][ix];     //  /(A-B);  COMMENTING OUT

        //*pReportEquation << aaafH[iz][iy][ix];
        //                         
        //if (ix == 0) pReportEquation << "\n"; else pReportEquation << " ";

      } // for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
    } // for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
  } // for (int iz=-halfwidth[2]; iz<=halfwidth[2]; iz++) {


  // COMMENTING OUT the factor of 1/(A-B):
  //A = A/(A-B);
  //B = B/(A-B);

  if (pA && pB) {
    *pA = A; // Rescale A and B numbers returned to the caller
    *pB = B; // (because we divided the array entries by (A-B) earlier)
  }

  if (pReportEquation) {
    *pReportEquation << "\n";
    *pReportEquation << " Filter Used:\n"
      " h(x,y,z)   = h_a(x,y,z) - h_b(x,y,z)\n"
      " h_a(x,y,z) = A*exp(-((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2)^(m/2))\n"
      " h_b(x,y,z) = B*exp(-((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2)^(n/2))\n"
      "  ... where      A = " << A << "\n"
      "                 B = " << B << "\n" 
      "                 m = " << m_exp << "\n"
      "                 n = " << n_exp << "\n" 
      "   (a_x, a_y, a_z) = "
                     << "(" << width_a[0]
                     << " " << width_a[1]
                     << " " << width_a[2] << ")\n"
      "   (b_x, b_y, b_z) = "
                     << "(" << width_b[0]
                     << " " << width_b[1]
                     << " " << width_b[2] << ")\n";
    *pReportEquation << " You can plot a slice of this function\n"
                     << "     in the X direction using:\n"
      " draw_filter_1D.py -dogg " << A << " " << B
                     << " " << width_a[0] << " " << width_b[0]
                     << " " << m_exp << " " << n_exp << endl;
    if ((width_a[1] != width_a[0]) || (width_a[2] != width_a[0])) {
      *pReportEquation << " and in the Y direction using:\n"
        " draw_filter_1D.py -dogg " << A << " " << B
                       << " " << width_a[1] << " " << width_b[1]
                       << " " << m_exp << " " << n_exp << endl;
      *pReportEquation << " and in the Z direction using:\n"
        " draw_filter_1D.py -dogg " << A << " " << B
                       << " " << width_a[2] << " " << width_b[2]
                       << " " << m_exp << " " << n_exp << endl;
    }
  } //if (pReportEquation)
  
  return filter;
} //_GenFilterDogg3D()




/// @brief Create a 3D filter and fill it with a difference of (generalized) Gaussians
///   @verbatim h(x,y,z) = A*exp(-(r/a)^m) - B*exp(-(r/b)^n)  @endverbatim
/// where  @verbatim r = sqrt(x^2 + y^2 + z^2) @endverbatim
///   and "A" and "B" are determined by normalization of each term independently
/// DEPRECIATION WARNING:  It's not clear if this type if filter is useful.
///                        I may delete this function in the future.

template<typename Scalar>

Filter3D<Scalar, int> 
GenFilterDogg3D(const Scalar width_a[3],   //!< "a" parameter in formula
                const Scalar width_b[3],   //!< "b" parameter in formula
                Scalar m_exp,              //!< "m" parameter in formula
                Scalar n_exp,              //!< "n" parameter in formula
                const int halfwidth[3],    //!< the width of the filter
                Scalar *pA=nullptr,     //!< optional:report A,B coeffs to user
                Scalar *pB=nullptr,     //!< optional:report A,B coeffs to user
                ostream *pReportEquation = nullptr //!< optional: print params used?
                )
{
  Filter3D<Scalar, int> filter_A =
    GenFilterGenGauss3D(width_a,      //"a_x", "a_y" gaussian width parameters
                        m_exp,        //"m" exponent parameter
                        halfwidth);

  Filter3D<Scalar, int> filter_B =
    GenFilterGenGauss3D(width_b,      //"b_x", "b_y" gaussian width parameters
                        n_exp,        //"n" exponent parameter
                        halfwidth);

  // The next function is defind in "filter3d_implementation.hpp"
  return _GenFilterDogg3D(width_a,
                          width_b,
                          m_exp,
                          n_exp,
                          filter_A, filter_B,
                          pA,
                          pB,
                          pReportEquation);
} //GenFilterDogg3D(...halfwidth...)



/// @brief This is a version of the function of the same name defined in 
/// "filter3d.hpp".  This version allows the user to control the width of the 
/// interval considered for the filter in 2 ways: "filter_truncate_ratio", 
/// which specifies the width of the interval in units of the "width_a/b[]"
/// (ie. a_x, b_x, ...) parameters, OR the "filter_truncate_threshold" which 
/// specifies how much the filter must decay before the filter is cutoff.
/// This function was not intended for public use.

template<typename Scalar>
// Create a 3-D filter and fill it with a difference of (generalized) Gaussians:
Filter3D<Scalar, int> 
GenFilterDogg3D(const Scalar width_a[3],  //"a" parameter in formula
                const Scalar width_b[3],  //"b" parameter in formula
                Scalar m_exp,             //"m" parameter in formula
                Scalar n_exp,             //"n" parameter in formula
                Scalar filter_truncate_ratio,//how many sigma before cutoff?
                Scalar filter_truncate_threshold,//cutoff when decay below threshold
                Scalar *pA = nullptr, //optional:report A,B coeffs to user
                Scalar *pB = nullptr, //optional:report A,B coeffs to user
                ostream *pReportProgress = nullptr  //!< optional: print equation to the terminal
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

template<typename Scalar>
Scalar
ApplyGauss(const int image_size[3], 
           Scalar const *const* const *aaafSource,
           Scalar ***aaafDest,
           Scalar const *const *const *aaafMask,
           const Scalar sigma[3],
           Scalar filter_truncate_ratio,
           Scalar filter_truncate_threshold,
           bool normalize = true,
           ostream *pReportProgress = nullptr)
{

  if (filter_truncate_ratio <= 0) {
    assert(filter_truncate_threshold > 0.0);
    //    filter_truncate_threshold = exp(-(1/2)*filter_truncate_ratio^2);
    //    -> filter_truncate_ratio^2 = -2*log(filter_truncate_threshold)
    filter_truncate_ratio = sqrt(-2*log(filter_truncate_threshold));
  }
  return ApplyGauss(image_size, 
                    aaafSource,
                    aaafDest,
                    aaafMask,
                    sigma,
                    filter_truncate_ratio,
                    normalize,
                    pReportProgress);

} // ApplyGauss(..., filter_truncate_ratio, filter_truncate_threshold, ...)





/// @brief This is a version of the function of the same name defined in 
/// "filter3d.hpp".  This version allows the user to control the width of the 
/// interval considered for the filter in 2 ways: "filter_truncate_ratio", 
/// which specifies the width of the interval in units of the "sigma_a/b[]"
/// (ie. a_x, b_x, ...) parameters, OR the "filter_truncate_threshold" which 
/// specifies how much the filter must decay before the filter is cutoff.
/// This function was not intended for public use.

template<typename Scalar>
void
ApplyDog(const int image_size[3], //source image size
         Scalar const *const *const *aaafSource,  //source image
         Scalar ***aaafDest,     //save results here
         Scalar const *const *const *aaafMask,  //ignore voxels where mask==0
         const Scalar sigma_a[3],
         const Scalar sigma_b[3],
         Scalar filter_truncate_ratio,
         Scalar filter_truncate_threshold,
         Scalar *pA = nullptr,
         Scalar *pB = nullptr,
         ostream *pReportProgress = nullptr)
{
  Scalar ***aaafTemp; //temporary array to store partially processed tomogram

  aaafTemp = Alloc3D(image_size);

  Scalar A, B;        // let the user know what A B coefficients were used

  A = ApplyGauss(image_size,
                 aaafSource,
                 aaafDest,         // <-- save result here
                 aaafMask,
                 sigma_a,
                 filter_truncate_ratio,
                 filter_truncate_threshold,
                 true,
                 pReportProgress);

  B = ApplyGauss(image_size,
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
  Dealloc3D(aaafTemp);

  // Report the A and B normalization coefficients to the caller?
  if (pA)
    *pA = A;
  if (pB)
    *pB = B;

} // ApplyDog(...,filter_truncate_ratio,filter_truncate_threshold,...)





template<typename Scalar>
void
ApplyLog(const int image_size[3], //source image size
         Scalar const *const *const *aaafSource,   //source image
         Scalar ***aaafDest,     //save results here
         Scalar const *const *const *aaafMask,     //ignore voxels where mask==0
         const Scalar sigma[3],  //Gaussian width in x,y,z drections
         Scalar delta_sigma_over_sigma, //difference in Gauss widths
         Scalar filter_truncate_ratio,
         Scalar filter_truncate_threshold,
         Scalar *pA = nullptr,
         Scalar *pB = nullptr,
         ostream *pReportProgress = nullptr)
{

  if (filter_truncate_ratio < 0)
    //    filter_truncate_threshold = exp(-(1/2)*filter_truncate_ratio^2);
    //    -> filter_truncate_ratio^2 = -2*log(filter_truncate_threshold)
    filter_truncate_ratio = sqrt(-2*log(filter_truncate_threshold));

  ApplyLog(image_size,
           aaafSource,
           aaafDest,
           aaafMask,
           sigma,
           delta_sigma_over_sigma,
           filter_truncate_ratio,
           pA,
           pB,
           pReportProgress);

} // ApplyLog(..,filter_truncate_ratio,filter_truncate_threshold,...)




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

template<typename Scalar>
static void
LocalFluctuationsByRadius(const int image_size[3], //!< number of voxels in x,y,z directions
                          Scalar const *const *const *aaafSource, //!< original image
                          Scalar ***aaafDest, //!< store filtered image here (fluctuation magnitude)
                          Scalar const *const *const *aaafMask, //!< optional: if not nullptr then ignore voxel ix,iy,iz if aaafMask[iz][iy][ix]==0
                          const Scalar radius[3],  //!< radius (=sigma/√3) of neighborhooed over which we search in x,y,z directions (ellipsoidal shaped search area)
                          Scalar template_background_exponent=2, //!< exponent controlling sharpness of the (generalized) Gaussian (slow if != 2)
                          Scalar filter_truncate_ratio=2.5, //!< width over which we search is this many times larger than the gaussian width parameter (sigma)
                          Scalar filter_truncate_threshold=0.02, //!< alternatively, specify how much the Gaussian must decay before we truncate it
                          bool normalize = true, //!< normalize the result?
                          ostream *pReportProgress = nullptr //!< report progress to the user?
                          )
{
  if (filter_truncate_ratio < 0.0)
    // Choose the filter domain window based on the "truncate_threshold"
    //    filter_truncate_threshold = exp(-(filter_truncate_ratio)^exponent);
    //    -> filter__truncate_ratio^exponent = -log(filter_truncate_threshold)
    filter_truncate_ratio = pow(-log(filter_truncate_threshold),
                                1.0 / template_background_exponent);

  LocalFluctuationsByRadius(image_size,
                            aaafSource,
                            aaafDest,
                            aaafMask,
                            radius,
                            template_background_exponent,
                            filter_truncate_ratio,
                            normalize,
                            pReportProgress);

} //LocalFluctuations()




} //namespace visfd




#endif //#ifndef _FILTER_3D_VARIANTS_HPP
