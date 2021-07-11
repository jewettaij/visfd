///   @file filter2d.hpp
///   @brief classes and functions that apply filters to 2D arrays
///   @author Andrew Jewett
///   @date 2018-2-26

#ifndef _FILTER2D_HPP
#define _FILTER2D_HPP

#include <cstring>
#include <ostream>
#include <cassert>
#include <cmath>
#include <limits>
using namespace std;
#include <alloc2d.hpp>



namespace visfd {


/// @class Filter2D
/// @brief A simple class for general linear (convolutional) filters in 2D
///
/// @note  In practice, this class is not used often because most operations
///        we need can be performed with separable filters which are faster.

template<typename Scalar, typename Integer>

class Filter2D {
public:
  Scalar *afH;
  Scalar **aafH;
  Integer halfwidth[2]; //num pixels from filter center to edge in x,y directions
  Integer array_size[2]; //size of the array in x,y directions (in voxels)



  /// @brief  Apply the filter to a 2D image (aafSource[][]).
  /// Save the results in the "aafDest" array.  (A "mask" is optional.)
  /// All arrays are 2D and assumed to be the same size.
  ///     
  /// @code
  /// If (normalize == false), the filter computes the convolution of h and f:
  ///        ___
  ///        \
  /// g[i] = /__  h[j] * f[i-j] * mask[i-j]
  ///         j
  ///     (sum over the width of the filter)
  ///
  /// ...where "i", "j", and "h", "f", and "mask" are defind below.
  /// Otherwise, if (normalize == true) it computes:
  ///        ___                               /  ___
  ///        \                                /   \
  /// g[i] = /__  h[j] * f[i-j] * mask[i-j]  /    /__  h[j] * mask[i-j]
  ///         j                             /      j
  ///
  /// where: f[i] is the original image at position i (ix, iy)
  ///         i = shorthand for (ix,iy) = a location in the filtered image
  ///         j = shorthand for (jx,jy) is summed over the entries in the filter
  ///       i-j = shorthand for (ix-jx, iy-jy)
  ///        g[i] is the data after the filter has been applied
  ///        h[j] is the filter
  ///        mask[i] selects the pixels we care about (usually either 0 or 1)
  ///          (If not specified, assume it is 1 inside the image, 0 outside.)
  /// @endcode
  ///
  /// @param size_source contains size of the source image (in the x,y directions)
  /// @param aafSource[][] is the source array (source image) <==> "f[i]"
  /// @param aafDest[][] will store the image after filtering <==> "g[i]"
  /// @param aafMask[][]==0 whenever we want to ignore entries in afSource[]. Optional.
  /// @param normalize  This boolean parameter = true if you want to divide g[i] by the sum of the weights considered. Optional.
  /// (useful if the sum of your filter elements, h[j], is 1, and if the sum was
  ///  not complete because some entries lie outside the mask or the boundary.)

  void Apply(const Integer size_source[2],
             Scalar const *const *aafSource,
             Scalar **aafDest,
             Scalar const *const *aafMask = nullptr,
             bool normalize = false) const
  {
    Scalar *afDenominator = nullptr;
    Scalar **aafDenominator = nullptr;
    if (normalize)
      Alloc2D(size_source, &afDenominator, &aafDenominator);

    Apply(size_source,
          aafSource,
          aafDest,
          aafMask,
          aafDenominator);

    if (normalize) {
      for (Integer iy=0; iy < size_source[1]; iy++)
        for (Integer ix=0; ix < size_source[0]; ix++)
          if (aafDenominator[iy][ix] > 0.0)
            aafDest[iy][ix] /= aafDenominator[iy][ix];
      Dealloc2D(size_source, &afDenominator, &aafDenominator);
    }
  }



  /// @brief  Apply the filter to a 2D image (aafSource[][]).
  ///         This version is identical to the other version of Apply()
  ///         except that this version both d[i] and g[i] whenever
  ///         you supply a non-null afDenominator[] argument (see below).
  ///         It also does not normalize the result (by dividing g[i] / d[i]].
  ///     
  /// @code
  /// This function computes the convolution of h and f, defind as:
  ///        ___
  ///        \
  /// g[i] = /__  h[j] * f[i-j] * mask[i-j]
  ///         j
  ///     (sum over the width of the filter)
  ///
  /// ...where "i", "j", and "h", "f", and "mask" are defind below.
  /// If aafDenominator!=nullptr, then it also computes the "denominator":
  ///        ___
  ///        \
  /// d[i] = /__  h[j] * mask[i-j]
  ///         j
  ///
  /// where: f[i] is the original array of (source) data at position i (ix,iy)
  ///         i = shorthand for (ix,iy) = a location in the filtered image
  ///         j = shorthand for (jx,jy) is summed over the entries in the filter
  ///       i-j = shorthand for (ix-jx, iy-jy)
  ///        g[i] is the data after the filter has been applied
  ///        h[j] is the filter
  ///        d[i] is the "denominator" = sum of all the filter weights h[j]...
  ///             ...considered when calculating g[i] (weighted by the mask)
  ///        mask[i] selects the pixels we care about (usually either 0 or 1)
  ///          (If not specified, assume it is 1 inside the image, 0 outside.)
  /// @endcode
  ///
  /// @param size_source contains size of the source image (in the x,y directions)
  /// @param aafSource[][] is the source array (source image) <==> "f[i]"
  /// @param aafDest[][] will store the image after filtering <==> "g[i]"
  /// @param aafMask[][]==0 whenever we want to ignore entries in afSource[][]. Optional.
  /// @param aafDenominator[][] will store d[i] if you supply a non-null pointer

  void Apply(const Integer size_source[2],
             Scalar const *const *aafSource,
             Scalar **aafDest,
             Scalar const *const *aafMask = nullptr,
             Scalar **aafDenominator = nullptr) const
             //bool precompute_mask_times_source = true) const
  {

    #pragma omp parallel for collapse(2)
    for (Integer iy=0; iy<size_source[1]; iy++) {

      for (Integer ix=0; ix<size_source[0]; ix++) {

        if ((aafMask) && (aafMask[iy][ix] == 0.0)) {
          aafDest[iy][ix] = 0.0;
          if (aafDenominator[iy][ix])
            aafDenominator[iy][ix] = 0.0;
          continue;
        }
          
        Scalar g = 0.0;
        Scalar denominator = 0.0;

        for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {

          Integer iy_jy = iy-jy;
          if ((iy_jy < 0) || (size_source[1] <= iy_jy))
            continue;

          for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {

            Integer ix_jx = ix-jx;
            if ((ix_jx < 0) || (size_source[0] <= ix_jx))
              continue;

            Scalar filter_val = aafH[jy][jx];

            if (aafMask) {
              Scalar mask_val = aafMask[iy_jy][ix_jx];
              if (mask_val == 0.0)
                continue;
              filter_val *= mask_val;
            }
            //Note: The "filter_val" also is needed to calculate
            //      the denominator used in normalization.
            //      It is unusual to use a mask unless you intend
            //      to normalize the result later, but I don't enforce this

            Scalar delta_g = 
              filter_val * aafSource[iy_jy][ix_jx];

            g += delta_g;

            if (aafDenominator)
              denominator += filter_val;
          }
        }

        if (aafDenominator)
          aafDenominator[ix][iy] = denominator;

        aafDest[iy][ix] = g;
      }
    }
  } //Apply()


  void Normalize() {
    // Make sure the sum of the filter weights is 1
    Scalar total = 0.0;
    for (Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
      for (Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
        total += aafH[iy][ix];
    for (Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
      for (Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
        aafH[iy][ix] /= total;
  }


  void Alloc(Integer const set_halfwidth[2]) {
    //Integer array_size[2];
    for(Integer d=0; d < 2; d++) {
      halfwidth[d] = set_halfwidth[d];
      array_size[d] = 1 + 2*halfwidth[d];
    }
    Alloc2D(array_size, &afH, &aafH);
    for (Integer iy = 0; iy < array_size[1]; iy++)
      for (Integer ix = 0; ix < array_size[0]; ix++)
        aafH[iy][ix] = -1.0e38; //(if uninitiliazed memory read, we will know)

    //shift pointers to enable indexing from i = -halfwidth .. +halfwidth
    for (Integer iy = 0; iy < array_size[1]; iy++)
      aafH[iy] += halfwidth[0];
    aafH += halfwidth[1];
  }


  void Dealloc() {
    if (! aafH) {
      assert(! afH);
      Init();
      return;
    }
    //Integer array_size[2];
    for(Integer d=0; d < 2; d++) {
      array_size[d] = 1 + 2*halfwidth[d];
    }
    //shift pointers back to normal
    for (Integer iy = -halfwidth[1]; iy <= halfwidth[1]; iy++)
      aafH[iy] -= halfwidth[0];
    aafH -= halfwidth[1];
    //then deallocate
    Dealloc2D(array_size, &afH, &aafH);
    for(Integer d=0; d < 2; d++) {
      array_size[d] = -1;
      halfwidth[d] = -1;
    }
  }


  void Resize(Integer const set_halfwidth[2]) {
    //Integer array_size[2];
    if (aafH) {
      assert(afH);
      for(Integer d=0; d < 2; d++)
        array_size[d] = 1 + 2*halfwidth[d];
      Dealloc2D(array_size, &afH, &aafH);
    }
    for(Integer d=0; d < 2; d++) {
      halfwidth[d] = set_halfwidth[d];
      array_size[d] = 1 + 2*halfwidth[d];
    }
    Alloc2D(array_size, &afH, &aafH);
  }


  void Init() {
    halfwidth[0] = -1;
    halfwidth[1] = -1;
    array_size[0] = -1;
    array_size[1] = -1;
    afH = nullptr;
    aafH = nullptr;
  }


  Filter2D(Integer const set_halfwidth[2]) {
    Init();
    Resize(set_halfwidth);
  }


  Filter2D() {
    Init();
  }


  // destructor, copy and move constructor, swap, and assignment operator

  ~Filter2D() {
    Dealloc();
  }


  Filter2D(const Filter2D<Scalar, Integer>& source) {
    Init();
    Resize(source.halfwidth); // allocates and initializes afH and aafH
    //for(Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
    //  for(Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
    //    aafH[iy][ix] = source.aafH[iy][ix];
    // -- Use std:copy() instead: --
    std::copy(source.afH,
              source.afH + (array_size[0] * array_size[1]),
              afH);
  }


  void swap(Filter2D<Scalar, Integer> &other) {
    std::swap(afH, other.afH);
    std::swap(aafH, other.aafH);
    std::swap(halfwidth, other.halfwidth);
    std::swap(array_size, other.array_size);
  }

  // Move constructor (C++11)
  Filter2D(Filter2D<Scalar, Integer>&& other) {
    Init();
    this->swap(other);
  }

  // Using the "copy-swap" idiom for the assignment operator
  Filter2D<Scalar, Integer>&
    operator = (Filter2D<Scalar, Integer> source) {
    this->swap(source);
    return *this;
  }

}; // class Filter2D




/// @brief
/// Create a 2D filter and fill it with a "generalized Gaussian" function:
///
/// @code
///    h_xy(r) = A*exp(-r^m)
/// where   r  = sqrt((x/σ_x)^2 + (y/σ_y)^2)
///   and   A  is determined by normalization of the discrete sum
/// @endcode
///
/// @note "A" is equal to the value stored in the middle of the array,
///       The caller can determine what "A" is by looking at this value.

template<typename Scalar>

Filter2D<Scalar, int>
GenFilterGenGauss2D(const Scalar width[2],  //"σ_x", "σ_y" parameters
                    Scalar m_exp,           //"m" exponent parameter
                    const int halfwidth[2], //size of 2D array to store filter
                    Scalar *pA=nullptr      //optional:report A coeff to user
                    )
{
  Scalar truncate_threshold = 1.0;
  for (int d=0; d<2; d++) {
    Scalar h = ((width[d]>0) ? exp(-pow(halfwidth[d]/width[d], m_exp)) : 1.0);
    if (h < truncate_threshold)
      truncate_threshold = h;
  }

  Filter2D<Scalar, int> filter(halfwidth);
  Scalar total = 0;
  for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
    for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
      Scalar r = sqrt(SQR(ix/width[0]) + SQR(iy/width[1]));
      Scalar h = ((r>0) ? exp(-pow(r, m_exp)) : 1.0);
      if (std::abs(h) < truncate_threshold)
        h = 0.0; // this eliminates corner entries which fall below threshold
                 // (and eliminates anisotropic artifacts due to these corners)
                 // There's no reason to keep any entries less than min value.

      filter.aafH[iy][ix] = h;
                 
      total += h;
    }
  }

  // normalize:
  for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
    for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
      filter.aafH[iy][ix] /= total;

  return filter;
} //GenFilterGenGauss2D(width, m_exp, halfwidth)






/// @brief
/// Create a 2D filter and fill it with a "generalized Gaussian" function.
///
/// @code
///    h_xy(r) = A*exp(-r^m)
/// where   r  = sqrt((x/σ_x)^2 + (y/σ_y)^2)
///   and   A  is determined by normalization of the discrete sum
/// @endcode
/// @note "A" is equal to the value stored in the middle of the array,
///       The caller can determine what "A" is by looking at this value.
/// @note The width of the filter window in each direction equals the "width"
///       in that direction multiplied by "filter_cutoff_ratio"
template<typename Scalar>
Filter2D<Scalar, int>
GenFilterGenGauss2D(const Scalar width[2], //!< "σ_x", "σ_y" parameters
                    Scalar m_exp,          //!< "m" parameter in formula
                    Scalar filter_cutoff_ratio=2.5, //!< how many sigma (σ) before truncating?
                    Scalar *pA=nullptr     //!< optional:report A coeff to user?
                    )
{
  // choose the width of the filter window based on the filter_cutoff_ratio
  int halfwidth[2];
  int ix = 0;

  for (int d=0; d<2; d++)
    halfwidth[d] = floor(width[d]*filter_cutoff_ratio);

  return GenFilterGenGauss2D(width,
                             m_exp,
                             halfwidth,
                             pA);

} //GenFilterGenGauss2D(width, m_exp, filter_cutoff_ratio)



} //namespace visfd


#endif //#ifndef _FILTER2D_HPP
