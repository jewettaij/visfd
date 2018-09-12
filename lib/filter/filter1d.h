#ifndef _FILTER1D_H
#define _FILTER1D_H

#include <cstring>
using namespace std;


template<class RealNum >
inline RealNum ABS(RealNum x) { return ((x<0.0) ? -x: x); }

template<class RealNum >
inline RealNum SQR(RealNum x) { return x*x; }

template<class RealNum >
inline RealNum MAX(RealNum x, RealNum y) { return ((x<y) ? y : x); }



template<class RealNum, class Integer>

class Filter1D {
public:
  RealNum *afH;
  Integer halfwidth; //distance from the filter center to its edge, in pixels
  Integer array_size; //size of the array in x,y directions (in voxels)


  /// @brief  Apply the filter to data in the original source array ("afSource")
  /// Save the results in the "aafDest" array.  (A "mask" is optional.)
  ///
  /// @code
  /// If afMask == NULL, then filter computes g(i):
  ///        ___
  ///        \
  /// g(i) = /__  h(j) * f(i-j)
  ///         j
  /// Otherwise, if afMask!=NULL and normalize == true, it computes
  ///        ___                               /  ___
  ///        \                                /   \
  /// g(i) = /__  h(j) * f(i-j) * mask(i-j)  /    /__  h(j) * mask(i-j)
  ///         j                             /      j
  ///
  /// where: f(i) is the original array of (source) data at position i
  ///        g(i) is the data after the filter has been applied
  ///        h(j) is the filter
  ///        mask(i) selects the entries we care about (usually either 0 or 1)
  /// @endcode
  ///
  /// @param size_source contains size of the source image (in the x,y directions)
  /// @param afSource[] is the original source data <==> h(j)
  /// @param afDest[] will store the result after filtering <==> g(j)
  /// @param afMask[]==0 whenever we want to ignore entries in afSource[]. Optional.
  /// @param normalize  This boolean parameter = true if you want to divide g(i) by the sum of the weights considered. Optional.
  /// (useful if the sum of your filter elements, h(j), is 1, and if the sum was
  ///  not complete because some entries lie outside the mask or the boundary.)

  void Apply(Integer const size_source, 
             RealNum const *afSource,
             RealNum *afDest,
             RealNum const *afMask = NULL,
             bool normalize = false) const
  {
    RealNum *afDenominator = NULL;
    if (normalize)
      afDenominator = new RealNum[size_source];

    Apply(size_source, 
          afSource,
          afDest,
          afMask,
          afDenominator);

    if (normalize) {
      for (int i=0; i < size_source; i++)
        if (afDenominator[i] > 0.0)
          afDest[i] /= afDenominator[i];
      delete [] afDenominator;
    }
  }


  /// @brief  Apply the filter to data in the original source array ("afSource")
  ///         This version is identical to the other version of Apply()
  ///         except that this version both d(i) and g(i) whenever
  ///         you supply a non-NULL afDenominator[] argument (see below).
  ///         It also does not normalize the result (by dividing g(i) / d(i)).
  /// @code
  /// If afMask == NULL, then filter computes g(i):
  ///        ___
  ///        \
  /// g(i) = /__  h(j) * f(i-j)
  ///         j
  /// Otherwise, if afMask!=NULL and afDenominator!=NULL, it computes g(i), d(i)
  ///        ___
  ///        \
  /// g(i) = /__  h(j) * f(i-j) * mask(i-j)
  ///         j
  ///        ___
  ///        \
  /// d(i) = /__  h(j) * mask(i-j)
  ///         j
  ///
  /// where: f(i) is the original array of (source) data at position i
  ///        g(i) is the data after the filter has been applied
  ///        h(j) is the filter
  ///        mask(i) is usually either 0 or 1
  ///        d(i) is the "denominator" = sum of the filter weights after masking
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
  /// @param afDenominator[] will store d(i) if you supply a non-NULL pointer.


  void Apply(Integer const size_source, 
             RealNum const *afSource,
             RealNum *afDest,
             RealNum const *afMask = NULL,
             RealNum *afDenominator = NULL) const
  {
    for (Integer ix=0; ix<size_source; ix++) {

      if ((afMask) && (afMask[ix] == 0.0)) {
        afDest[ix] = 0.0;
        continue;
      }
          
      RealNum g = 0.0;
      RealNum denominator = 0.0;

      for (Integer jx=-halfwidth; jx<=halfwidth; jx++) {

        int ix_jx = ix-jx;

        if ((ix_jx < 0) || (size_source <= ix_jx))
          continue;

        RealNum filter_val = afH[jx];
        if (afMask)
          filter_val *= afMask[ix_jx];
          //Note: The "filter_val" also is needed to calculate
          //      the denominator used in normalization.
          //      It is unusual to use a mask unless you intend
          //      to normalize the result later, but I don't enforce this

        RealNum delta_g = filter_val * afSource[ix_jx];

        g += delta_g;

        if (afDenominator)
          denominator += filter_val;
      }

      if (afDenominator)
        afDenominator[ix] = denominator;

      afDest[ix] = g;
    }
  } //Apply()


  
  void Normalize() {
    // Make sure the sum of the filter weights is 1
    RealNum total = 0.0;
    for (Integer ix=-halfwidth; ix<=halfwidth; ix++)
      total += afH[ix];
    for (Integer ix=-halfwidth; ix<=halfwidth; ix++)
      afH[ix] /= total;
  }


  void Init() {
    halfwidth = -1;
    afH = NULL;
  }


  void Alloc(Integer set_halfwidth) {
    halfwidth = set_halfwidth;
    array_size = 1 + 2*halfwidth;
    afH = new RealNum [array_size];
    for (int i = 0; i < array_size; i++)
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


  inline Filter1D(const Filter1D<RealNum, Integer>& source) {
    Init();
    Resize(source.halfwidth); // allocates and initializes afH
    //for(Int ix=-halfwidth; ix<=halfwidth; ix++)
    //  afH[ix] = source.afH[ix];
    // -- Use memcpy() instead: --
    //memcpy(afH,
    //       source.afH,
    //       array_size * sizeof(RealNum));
    // -- Use std:copy() instead: --
    std::copy(source.afH, source.afH + array_size, afH);
  }


  ~Filter1D() {
    Dealloc();
  }

  inline void swap(Filter1D<RealNum, Integer> &other) {
    std::swap(afH, other.afH);
    std::swap(halfwidth, other.halfwidth);
    std::swap(array_size, other.array_size);
  }

  inline Filter1D<RealNum, Integer>&
    operator = (Filter1D<RealNum, Integer> source) {
    this->swap(source);
    return *this;
  }


}; // class Filter1D



//void swap(Filter1D<RealNum, Integer> &a, Filter1D<RealNum, Integer> &b) {
//  a.swap(b);
//}



template<class RealNum>
// GenFilterGauss1D generates a 1-D filter and fills its array with values
// corresponding to a normalized Gaussian evaluated at evenly spaced intervals.
// The caller must specify the "σ" parameter (width of the Gaussian,
// in units of pixels/voxels), in addition to the "halfwidth" parameter, which
// indicates the number of entries in the array (in units of pixels/voxels).
Filter1D<RealNum, int>
GenFilterGauss1D(RealNum sigma,  // The "σ" paramgeter in the Gaussian
                 int halfwidth,  // number of entries in the filter array / 2
                 ostream *pReportProgress = NULL)
{
  Filter1D<RealNum, int> filter(halfwidth);

  RealNum sum = 0.0;
  for (int i=-halfwidth; i<=halfwidth; i++) {
    if (sigma == 0.0) //(When sigma==0, use a Kronecker delta function)
      filter.afH[i] = ((i == 0) ? 1.0 : 0.0);
    else
      filter.afH[i] = exp(-(i*i)/(2.0*sigma*sigma)); // will normalize later
    sum += filter.afH[i];
  }

  // normalize:
  for (int i=-halfwidth; i<=halfwidth; i++)
    filter.afH[i] /= sum;
  return filter;
} //GenFilterGauss1D(sigma, halfwidth)



#endif //#ifndef _FILTER1D_H
