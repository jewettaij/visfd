///   @file filter1d.hpp
///   @brief classes and functions that apply filters to 1D arrays
///   @author Andrew Jewett
///   @date 2018-2-26

#ifndef _FILTER1D_HPP
#define _FILTER1D_HPP

#include <cstring>
#include <ostream>
#include <cassert>
#include <cmath>
#include <limits>
using namespace std;


namespace visfd {





template<typename Scalar, typename Integer>

class Filter1D {
public:
  Scalar *afH;
  Integer halfwidth; //distance from the filter center to its edge, in pixels
  Integer array_size; //size of the array in x,y directions (in voxels)


  /// @brief  Calculate the convolution of f[] (afSource) and h[] (afH):
  /// @code
  ///        ___
  ///        \
  /// g[i] = /__  h[j] * f[i-j] * Theta[i-j]
  ///         j
  ///     (sum over the width of the filter)
  ///     Theta[i] = 1 if inside the array (0<=i<size_source), 0 otherwise
  /// @endcode
  /// Save the results ("g[i]") in the "afDest" array.
  /// (Sparse input optimizations: Array entries far away from non-zero f[i] 
  ///  values are skipped and set to 0, without performing any filtering. This
  ///  only adds a small amount of overhead when filtering non-sparse arrays.)

  void Apply(Integer size_source, 
             Scalar const *afSource,
             Scalar *afDest)
  {
    assert(afDest != afSource);

    // -- Sparse input optimaztion: --
    // Scan along the entries of afSource[] array, using a for loop of the form
    //   for (I=0, I<size_source; I++)
    // Let "m" denote the number entries visited since we last encountered
    // a non-zero entry in the afSource[] array.

    Integer m = array_size; // initialize with a large number

    // Once "m" reaches "array_size"(a precomputed value equaling 2*halfwidth+1)
    // THEN we know that the array entry located at position "i", defined as:
    //   i = I - halfwidth
    // is at least "halfwidth" entries away from the nearest non-zero afSource[]
    // array entry.  In that case we don't need to apply the fiter at this
    // location (i), because any array entries which could effect the result
    // at (i) lie beyond the range of the filter.
   
    // First look-ahead "halfwidth" entries to calculate "m" beforehand:
    Integer init_width = halfwidth;
    if (init_width > size_source)
      init_width = size_source;
    for (Integer I=0; I<init_width; I++) {
      if (afSource[I] != 0.0)
        m = 0;
      else
        m++;
    }

    // Then loop over the entries in the table, updating "m" and performing
    // the filtering operation (sum) only when necessary.
    Integer I = halfwidth;
    for (Integer i=0; i<size_source; i++) {

      // update m
      if ((I < size_source) && (afSource[I] != 0.0))
        m = 0;
      else
        m++;
      I++;
      if (m >= array_size) {
        afDest[i] = 0.0;
        continue;
      }
      Scalar g = 0.0;
      for (Integer j=-halfwidth; j<=halfwidth; j++) {
        Integer i_j = i-j;
        if ((i_j < 0) || (size_source <= i_j))
          continue;
        g += afH[j] * afSource[i_j];
      }
      afDest[i] = g;
    }
  } //Apply()




  /// @brief  Apply the filter to data in the original source array ("afSource")
  /// Save the results in the "aafDest" array.
  ///
  /// @code
  /// If (normalize == false), the filter computes the convolution of h and f:
  ///        ___
  ///        \
  /// g[i] = /__  h[j] * f[i-j] * mask[i-j]
  ///         j
  ///     (sum over the width of the filter)
  /// Otherwise, if (normalize == true) it computes:
  ///        ___                               /  ___
  ///        \                                /   \
  /// g[i] = /__  h[j] * f[i-j] * mask[i-j]  /    /__  h[j] * mask[i-j]
  ///         j                             /      j
  ///
  /// where: f[i] is the original array of (source) data at position i
  ///        g[i] is the data after the filter has been applied
  ///        h[j] is the filter
  ///        mask[i] selects the entries we care about (usually either 0 or 1)
  ///          (If not specified, assume it is 1 inside the array, 0 outside.)
  /// @endcode
  ///
  /// @param size_source contains size of the source image (in the x,y directions)
  /// @param afSource[] is the original source data <==> h(j)
  /// @param afDest[] will store the result after filtering <==> g(j)
  /// @param afMask[]==0 whenever we want to ignore entries in afSource[]. Optional.
  /// @param normalize  This boolean parameter = true if you want to
  ///  divide g(i) by the sum of the weights considered. Optional.
  /// (useful if the sum of your filter elements, h(j), is 1, and if the sum was
  ///  not complete because some entries lie outside the mask or the boundary.)


  void Apply(Integer size_source, 
             Scalar const *afSource,
             Scalar *afDest,
             Scalar const *afMask,
             bool normalize) const
  {
    Scalar *afDenominator = nullptr;
    if (normalize)
      afDenominator = new Scalar[size_source];

    Apply(size_source, 
          afSource,
          afDest,
          afMask,
          afDenominator);

    if (normalize) {
      for (Integer i=0; i < size_source; i++)
        if (afDenominator[i] > 0.0)
          afDest[i] /= afDenominator[i];
      delete [] afDenominator;
    }
  }


  /// @brief  Apply the filter to data in the original source array ("afSource")
  ///         This version is identical to the other version of Apply() except
  ///         that this version returns g[i] AND the normalization d[i] whenever
  ///         you supply a non-null afDenominator[] argument (see below).
  ///         This version does NOT normalize the result (by dividing g[i]/d[i])
  /// @code
  /// This function computes both g[i] and d[i]  (d[i] is optional)  where:
  ///        ___
  ///        \
  /// g[i] = /__  h[j] * mask[i-j] * f[i-j]
  ///         j
  ///        ___
  ///        \
  /// d[i] = /__  h[j] * mask[i-j]
  ///         j
  ///     (sum over the width of the filter)
  ///
  /// where: f(i) is the original array of (source) data at position i
  ///        g(i) is the data after the filter has been applied
  ///        h(j) is the filter
  ///        mask[i] selects the entries we care about (usually either 0 or 1)
  ///          (If not specified, assume it is 1 inside the array, 0 outside.)
  ///        d[i] is the "denominator" = sum of all the filter weights h[j]...
  ///             ...considered when calculating g[i] (weighted by the mask)
  /// @endcode
  ///       
  ///  Note on mask functions:
  ///     When summing over j, we ignore contributions from positions where
  ///     mask(i-j) is zero.  We don't count them in the average.
  ///
  /// @param size_source contains size of the source image (in the x,y directions)
  /// @param afSource[] is the original source data <==> h(j)
  /// @param afDest[] will store the result after filtering <==> g(j)
  /// @param afMask[]==0 whenever we want to ignore entries in afSource[]. Optional.
  /// @param afDenominator[] will store d(i) if you supply a non-null pointer.


  void Apply(Integer const size_source, 
             Scalar const *afSource,
             Scalar *afDest,
             Scalar const *afMask,
             Scalar *afDenominator = nullptr) const
  {
    assert(afDest != afSource);
    assert(afDest != afMask);

    // -- Sparse input optimaztion: --
    // Scan along the entries of afMask[] array, using a for loop of the form
    //   for (I=0, I<size_source; I++)
    // Let "m" denote the number entries visited since we last encountered
    // a non-zero entry in the afMask[] array.

    Integer m = array_size; // initialize with a large number

    // Once "m" reaches "array_size"(a precomputed value equaling 2*halfwidth+1)
    // THEN we know that the array entry located at position "i", defined as:
    //   i = I - halfwidth
    // is at least "halfwidth" entries away from the nearest non-zero afMask[]
    // array entry.  In that case we don't need to apply the fiter at this
    // location (i), because any array entries which could effect the result
    // at (i) lie beyond the range of the filter.
    
    // First look-ahead "halfwidth" entries to calculate "m" beforehand:
    Integer init_width = halfwidth;
    if (init_width > size_source)
      init_width = size_source;
    for (Integer I=0; I<init_width; I++) {
      if ((afMask == nullptr) || (afMask[I] != 0.0))
        m = 0;
      else
        m++;
    }

    // Then loop over the entries in the table, updating "m" and performing
    // the filtering operation (sum) only when necessary.
    Integer I =  halfwidth;
    for (Integer i=0; i<size_source; i++) {

      // update m
      if ((I < size_source) && ((afMask == nullptr) || (afMask[I] != 0.0)))
        m = 0;
      else
        m++;
      I++;
      if (m >= array_size) {
        afDest[i] = 0.0;
        if (afDenominator)
          afDenominator[i] = 0.0;
        continue;
      }

      // If we got this far, then there is a nearby non-zero entry in the
      // afMask[] array, so we proceed to apply the filter at this location, i

      Scalar g = 0.0;           // g[i], the result after filtering
      Scalar denominator = 0.0; // "d[i]" in the comments above

      // Inner loop:

      for (Integer j=-halfwidth; j<=halfwidth; j++) {

        Integer i_j = i-j;

        if ((i_j < 0) || (size_source <= i_j))
          continue;

        Scalar filter_val = afH[j];
        if (afMask)
          filter_val *= afMask[i_j];
          //Note: The "filter_val" also is needed to calculate
          //      the denominator used in normalization.
          //      It is unusual to use a mask unless you intend
          //      to normalize the result later, but I don't enforce this

        Scalar delta_g = filter_val * afSource[i_j];

        g += delta_g;

        if (afDenominator)
          denominator += filter_val;
      }

      if (afDenominator)
        afDenominator[i] = denominator;

      afDest[i] = g;
    } //for (Integer i=0; i<size_source; i++)

  } //Apply()


  
  void Normalize() {
    // Make sure the sum of the filter weights is 1
    Scalar total = 0.0;
    for (Integer i=-halfwidth; i<=halfwidth; i++)
      total += afH[i];
    for (Integer i=-halfwidth; i<=halfwidth; i++)
      afH[i] /= total;
  }


  void Init() {
    halfwidth = -1;
    afH = nullptr;
  }


  void Alloc(Integer set_halfwidth) {
    halfwidth = set_halfwidth;
    array_size = 1 + 2*halfwidth;
    afH = new Scalar [array_size];
    for (Integer i = 0; i < array_size; i++)
      afH[i] = -1.0e38; //(if uninitiliazed memory read, we will know)

    //shift pointers to enable indexing from i = -halfwidth .. +halfwidth
    afH += halfwidth;
  }


  void Dealloc() {
    if (afH) {
      //indexing starts at -halfwidth shift pointers back to normal
      afH -= halfwidth;
      delete [] afH;
    }
    halfwidth = -1;
    array_size = -1;
  }


  void Resize(Integer set_halfwidth) {
    Dealloc();
    Alloc(set_halfwidth);
  }


  Filter1D() {
    Init();
  }


  Filter1D(Integer halfwidth) {
    Init();
    Resize(halfwidth);
  }


  Filter1D(const Filter1D<Scalar, Integer>& source) {
    Init();
    Resize(source.halfwidth); // allocates and initializes afH
    //for(Integer i=-halfwidth; i<=halfwidth; i++)
    //  afH[i] = source.afH[i];
    // -- Use std:copy() instead: --
    std::copy(source.afH, source.afH + array_size, afH);
  }


  ~Filter1D() {
    Dealloc();
  }

  void swap(Filter1D<Scalar, Integer> &other) {
    std::swap(afH, other.afH);
    std::swap(halfwidth, other.halfwidth);
    std::swap(array_size, other.array_size);
  }

  Filter1D<Scalar, Integer>&
    operator = (Filter1D<Scalar, Integer> source) {
    this->swap(source);
    return *this;
  }

}; // class Filter1D




//void swap(Filter1D<Scalar, Integer> &a, Filter1D<Scalar, Integer> &b) {
//  a.swap(b);
//}




/// @brief
/// This function generates a 1-D filter and fills its array with values
/// corresponding to a normalized Gaussian evaluated at evenly spaced intervals.
/// The caller must specify the "σ" parameter (width of the Gaussian,
/// in units of pixels/voxels), in addition to the "halfwidth" parameter, which
/// indicates the number of entries in the array (in units of pixels/voxels).

template<typename Scalar>

Filter1D<Scalar, int>
GenFilterGauss1D(Scalar sigma,  //!< The "σ" paramgeter in the Gaussian
                 int halfwidth,  //!< number of entries in the filter array / 2
                 ostream *pReportProgress = nullptr //!< report filter details to the user?  (WARNING: This argument is currently ignored.)
                 )
{
  Filter1D<Scalar, int> filter(halfwidth);

  long double sum = 0.0;
  for (int i=-halfwidth; i<=halfwidth; i++) {
    if (sigma == 0.0) //(When sigma==0, use a Kronecker delta function)
      filter.afH[i] = ((i == 0) ? 1.0 : 0.0);
    else {

      // Old method, using an ordinary continuous Gaussian kernel:
      //filter.afH[i] = std::exp(-(i*i)/(2.0*sigma*sigma));
      //
      // New method, using a discrete Gaussian kernel:
      // https://en.wikipedia.org/wiki/Scale_space_implementation#The_discrete_Gaussian_kernel

      // Use long double's instead of floats during the intermediate
      // stages of the calculation in order to increase numerical precision.
      long double S = sigma;
      long double I = i;

      if ((S <= 10.0) && (abs(I) <= 20.0)) {
        // For small S and I, use a discrete Gaussian kernel
        long double h = exp(-S*S) * std::cyl_bessel_i(std::abs(I), S*S);
        filter.afH[i] = h;
      }
      else {
        // If either I or S > 10.0, then the bessel-function formula
        // may not be numerically stable.  In these cases, it is probably
        // it's numerically more accurate to approximate the result
        // using the formula for an ordinary (continuous) Gaussian.
        filter.afH[i] = std::exp(-(I*I)/(2.0*S*S)) / std::sqrt(2*S*S*M_PI);
      }


      // (We will make sure this is normalized later.)
    }
    sum += filter.afH[i];
  }

  // normalize:
  for (int i=-halfwidth; i<=halfwidth; i++)
    filter.afH[i] /= sum;

  return filter;
} //GenFilterGauss1D(sigma, halfwidth)



} //namespace visfd



#endif //#ifndef _FILTER1D_HPP
