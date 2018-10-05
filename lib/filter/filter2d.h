#ifndef _FILTER2D_H
#define _FILTER2D_H

#include <cstring>
#include <ostream>
using namespace std;
#include <alloc2d.h>





template<class RealNum, class Integer>

class Filter2D {
public:
  RealNum *afH;
  RealNum **aafH;
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
  /// g[i] = /__  h[j] * f[i-j] * Theta[i-j]
  ///         j
  ///     (sum over the width of the filter)
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
  ///          (If not supplied, assumed it is 1 everywhere inside, 0 outside)
  ///        Theta[i] = 1 if i is inside the image boundaries, 0 otherwise.
  ///          (In other words, we only consider pixels within the image.)
  /// @endcode
  ///
  /// @param size_source contains size of the source image (in the x,y directions)
  /// @param aafSource[][] is the source array (source image) <==> "f[i]"
  /// @param aafDest[][] will store the image after filtering <==> "g[i]"
  /// @param aafMask[][]==0 whenever we want to ignore entries in afSource[]. Optional.
  /// @param normalize  This boolean parameter = true if you want to divide g[i] by the sum of the weights considered. Optional.
  /// (useful if the sum of your filter elements, h[j], is 1, and if the sum was
  ///  not complete because some entries lie outside the mask or the boundary.)

  void Apply(Integer const size_source[2],
             RealNum const *const *aafSource,
             RealNum **aafDest,
             RealNum const *const *aafMask = NULL,
             bool normalize = false) const
  {
    RealNum *afDenominator = NULL;
    RealNum **aafDenominator = NULL;
    if (normalize)
      Alloc2D(size_source, &afDenominator, &aafDenominator);

    Apply(size_source,
          aafSource,
          aafDest,
          aafMask,
          aafDenominator);

    if (normalize) {
      for (int iy=0; iy < size_source[1]; iy++)
        for (int ix=0; ix < size_source[0]; ix++)
          if (aafDenominator[iy][ix] > 0.0)
            aafDest[iy][ix] /= aafDenominator[iy][ix];
      Dealloc2D(size_source, &afDenominator, &aafDenominator);
    }
  }



  /// @brief  Apply the filter to a 2D image (aafSource[][]).
  ///         This version is identical to the other version of Apply()
  ///         except that this version both d[i] and g[i] whenever
  ///         you supply a non-NULL afDenominator[] argument (see below).
  ///         It also does not normalize the result (by dividing g[i] / d[i]].
  ///     
  /// @code
  /// If afMask == NULL, then filter computes g[i] and d[i]:
  ///        ___
  ///        \
  /// g[i] = /__  h[j] * f[i-j] * Theta[i-j]
  ///         j
  ///        ___
  ///        \
  /// d[i] = /__  h[j] * Theta[i-j]
  ///         j
  ///     (sum over the width of the filter)
  /// Otherwise, if afMask!=NULL and afDenominator!=NULL, it computes:
  ///        ___
  ///        \
  /// g[i] = /__  h[j] * f[i-j] * mask[i-j]
  ///         j
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
  ///        d[i] is the "denominator" = sum of the filter weights considered
  ///        mask[i] selects the pixels we care about (usually either 0 or 1)
  ///          (If not supplied, assumed it is 1 everywhere inside, 0 outside)
  ///        Theta[i] = 1 if i is inside the image boundaries, 0 otherwise.
  ///          (In other words, we only consider pixels within the image.)
  /// @endcode
  ///
  /// @param size_source contains size of the source image (in the x,y directions)
  /// @param aafSource[][] is the source array (source image) <==> "f[i]"
  /// @param aafDest[][] will store the image after filtering <==> "g[i]"
  /// @param aafMask[][]==0 whenever we want to ignore entries in afSource[][]. Optional.
  /// @param aafDenominator[][] will store d[i] if you supply a non-NULL pointer

  void Apply(Integer const size_source[2],
             RealNum const *const *aafSource,
             RealNum **aafDest,
             RealNum const *const *aafMask = NULL,
             RealNum **aafDenominator = NULL) const
             //bool precompute_mask_times_source = true) const
  {

    #pragma omp parallel for collapse(2)
    for (Integer iy=0; iy<size_source[1]; iy++) {

      for (Integer ix=0; ix<size_source[0]; ix++) {

        if ((aafMask) && (aafMask[iy][ix] == 0.0)) {
          aafDest[iy][ix] = 0.0;
          continue;
        }
          
        RealNum g = 0.0;
        RealNum denominator = 0.0;

        for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {

          int iy_jy = iy-jy;
          if ((iy_jy < 0) || (size_source[1] <= iy_jy))
            continue;

          for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {

            int ix_jx = ix-jx;
            if ((ix_jx < 0) || (size_source[0] <= ix_jx))
              continue;

            RealNum filter_val = aafH[jy][jx];

            if (aafMask)
              filter_val *= aafMask[iy_jy][ix_jx];
              //Note: The "filter_val" also is needed to calculate
              //      the denominator used in normalization.
              //      It is unusual to use a mask unless you intend
              //      to normalize the result later, but I don't enforce this

            RealNum delta_g = 
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
    RealNum total = 0.0;
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
    for (int iy = 0; iy < array_size[1]; iy++)
      for (int ix = 0; ix < array_size[0]; ix++)
        aafH[iy][ix] = -1.0e38; //(if uninitiliazed memory read, we will know)

    //shift pointers to enable indexing from i = -halfwidth .. +halfwidth
    aafH += halfwidth[1];
    for (int iy = 0; iy < array_size[1]; iy++)
      aafH[iy] += halfwidth[0];
  }


  void Dealloc() {
    //Integer array_size[2];
    for(Integer d=0; d < 2; d++) {
      array_size[d] = 1 + 2*halfwidth[d];
      halfwidth[d] = -1;
    }
    //shift pointers back to normal
    aafH -= halfwidth[1];
    for (int iy = 0; iy < array_size[1]; iy++)
      aafH[iy] -= halfwidth[0];
    //then deallocate
    Dealloc2D(array_size, &afH, &aafH);
    for(Integer d=0; d < 2; d++)
      array_size[d] = -1;
  }


  void Resize(Integer const set_halfwidth[2]) {
    //Integer array_size[2];
    if (afH && aafH) {
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


  inline Filter2D(const Filter2D<RealNum, Integer>& source) {
    Init();
    Resize(source.halfwidth); // allocates and initializes afH and aafH
    //for(Int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
    //  for(Int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
    //    aafH[iy][ix] = source.aafH[iy][ix];
    // -- Use memcpy() instead: --
    //memcpy(afH,
    //       source.afH,
    //       (array_size[0] * array_size[1])
    //       *sizeof(RealNum));
    // -- Use std:copy() instead: --
    std::copy(source.afH,
              source.afH + (array_size[0] * array_size[1]),
              afH);
  }


  void Init() {
    halfwidth[0] = -1;
    halfwidth[1] = -1;
    afH = NULL;
    aafH = NULL;
  }


  Filter2D(Integer const set_halfwidth[2]) {
    Init();
    Resize(set_halfwidth);
  }


  Filter2D() {
    Init();
  }


  ~Filter2D() {
    Dealloc();
  }


  inline void swap(Filter2D<RealNum, Integer> &other) {
    std::swap(afH, other.afH);
    std::swap(aafH, other.aafH);
    std::swap(halfwidth, other.halfwidth);
    std::swap(array_size, other.array_size);
  }


  inline Filter2D<RealNum, Integer>&
    operator = (Filter2D<RealNum, Integer> source) {
    this->swap(source);
    return *this;
  }

}; // class Filter2D





template<class RealNum>
// Create a 2D filter and fill it with a "generalized Gaussian" function:
//    h_xy(r) = A*exp(-r^m)
// where   r  = sqrt((x/σ_x)^2 + (y/σ_y)^2)
//   and   A  is determined by normalization of the discrete sum
// Note: "A" is equal to the value stored in the middle of the array,
//       The caller can determine what "A" is by looking at this value.
Filter2D<RealNum, int>
GenFilterGenGauss2D(RealNum width[2],    //"σ_x", "σ_y" parameters
                    RealNum m_exp,       //"m" exponent parameter
                    int halfwidth[2],
                    RealNum *pA=NULL,    //optional:report A coeff to user
                    ostream *pReportProgress = NULL)
{
  RealNum window_threshold = 1.0;
  for (int d=0; d<2; d++) {
    RealNum h = ((width[d]>0) ? exp(-pow(halfwidth[d]/width[d], m_exp)) : 1.0);
    if (h < window_threshold)
      window_threshold = h;
  }

  Filter2D<RealNum, int> filter(halfwidth);
  RealNum total = 0;
  for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
    for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
      RealNum r = sqrt(SQR(ix/width[0]) + SQR(iy/width[1]));
      RealNum h = ((r>0) ? exp(-pow(r, m_exp)) : 1.0);
      if (ABS(h) < window_threshold)
        h = 0.0; // this eliminates corner entries which fall below threshold
                 // (and eliminates anisotropic artifacts due to these corners)
                 // There's no reason to keep any entries less than min value.

      filter.aafH[iy][ix] = h;
                 
      total += h;
    }
  }

  // normalize:
  for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
    for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
      filter.aafH[iy][ix] /= total;
                 
      //FOR DEBUGGING REMOVE EVENTUALLY
      if (pReportProgress)
        *pReportProgress << "GenGauss2D:" //<< window_threshold
                         <<" aafH["<<iy<<"]["<<ix<<"] = "
                         << filter.aafH[iy][ix] << endl;
    }
  }
  return filter;
} //GenFilterGenGauss2D(width, m_exp, halfwidth)






template<class RealNum>
Filter2D<RealNum, int>
GenFilterGenGauss2D(RealNum width[2],            //"s_x", "s_y" parameters
                    RealNum m_exp,               //"m" parameter in formula
                    RealNum filter_cutoff_ratio,
                    RealNum *pA=NULL,    //optional:report A coeff to user
                    ostream *pReportProgress = NULL)
{
  // choose the width of the filter window based on the filter_cutoff_ratio
  int halfwidth[2];
  int ix = 0;

  for (int d=0; d<2; d++) {
    halfwidth[d] = floor(width[d]*filter_cutoff_ratio);
  }
  return GenFilterGenGauss2D(width,
                             m_exp,
                             halfwidth,
                             pA,
                             pReportProgress);
} //GenFilterGenGauss2D(width, m_exp, filter_cutoff_ratio)





template<class RealNum>
// Create a 2D filter and fill it with a difference of (generalized) Gaussians:
// This version requires that the caller has already created individual
// filters for the two gaussians.
// All this function does is subtract one filter from the other (and rescale).
Filter2D<RealNum, int> 
_GenFilterDogg2D(RealNum width_a[2],  //"a" parameter in formula
                 RealNum width_b[2],  //"b" parameter in formula
                 RealNum m_exp,  //"m" parameter in formula
                 RealNum n_exp,  //"n" parameter in formula
                 Filter2D<RealNum, int>& filterXY_A, //filters for the two
                 Filter2D<RealNum, int>& filterXY_B, //gaussians
                 RealNum *pA=NULL, //optional:report A,B coeffs to user
                 RealNum *pB=NULL, //optional:report A,B coeffs to user
                 ostream *pReportProgress = NULL)
{

  RealNum A, B;
  //A, B = height of the central peak
  A = filterXY_A.aafH[0][0];
  B = filterXY_B.aafH[0][0];

  // The "difference of gaussians" filter is the difference between
  // these two (generalized) gaussian filters.
  int halfwidth[2];
  halfwidth[0] = MAX(filterXY_A.halfwidth[0], filterXY_B.halfwidth[0]);
  halfwidth[1] = MAX(filterXY_A.halfwidth[1], filterXY_B.halfwidth[1]);
  Filter2D<RealNum, int> filter(halfwidth);

  if (pReportProgress)
    *pReportProgress  << "Array of 2D filter entries:" << endl;

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
                         


      //FOR DEBUGGING REMOVE EVENTUALLY
      if (pReportProgress)
        *pReportProgress << "GenDogg2D: aafH["<<iy<<"]["<<ix<<"] = "
                         << filter.aafH[iy][ix] << endl;

      // *pReportProgress  << aafH[iy][ix];
      //if (ix == 0)
      //  *pReportProgress << "\n";
      //else
      //  *pReportProgress << " ";

    } // for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
  } // for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {

  if (pReportProgress)
    *pReportProgress << "\n";

  // COMMENTING OUT 1/(A-B)
  //if (pA && pB) {
  //  *pA = A/(A-B); // Rescale A and B numbers returned to the caller
  //  *pB = B/(A-B); // (because we divided the array entries by (A-B) earlier)
  //}

  return filter;

} //_GenFilterDogg2D()



template<class RealNum>
// Create a 2D filter and fill it with a difference of (generalized) Gaussians:
Filter2D<RealNum, int> 
GenFilterDogg2D(RealNum width_a[2],  //"a" parameter in formula
                RealNum width_b[2],  //"b" parameter in formula
                RealNum m_exp,       //"m" parameter in formula
                RealNum n_exp,       //"n" parameter in formula
                int halfwidth[2],
                RealNum *pA = NULL,  //optional:report A,B coeffs to user
                RealNum *pB = NULL,  //optional:report A,B coeffs to user
                ostream *pReportProgress = NULL)
{
  Filter2D<RealNum, int> filterXY_A =
    GenFilterGenGauss2D(width_a,      //"a_x", "a_y" gaussian width parameters
                        m_exp,        //"n" exponent parameter
                        halfwidth);
                        //pReportProgress);

  Filter2D<RealNum, int> filterXY_B =
    GenFilterGenGauss2D(width_b,      //"b_x", "b_y" gaussian width parameters
                        n_exp,        //"n" exponent parameter
                        halfwidth);
                        //pReportProgress);

  return _GenFilterDogg2D(width_a,
                          width_b,
                          m_exp,
                          n_exp,
                          filterXY_A, filterXY_B,
                          pA,
                          pB,
                          pReportProgress);
} //GenFilterDogg2D(...halfwidth...)



#endif //#ifndef _FILTER2D_H
