///   @file filter3d.h
///   @brief a collection of common filters for 3D arrays
///   @author Andrew Jewett
///   @date 2018-9-14

#ifndef _FILTER3D_H
#define _FILTER3D_H

#include <cstring>
#include <cmath>
#include <ostream>
#include <vector>
#include <tuple>
#include <cassert>
using namespace std;
#include <alloc3d.h>
#include <filter1d.h>  // defines "Filter1D" (used in "ApplyGauss()")
#include <filter3d_utils.h> // defines "AverageArr()", and similar functions...
#include <eigen3_simple.h>






/// @brief A simple class for general linear convolutional filters in 3D
///
/// @note  In practice, this class is not used often because separable filters
///        based on Gaussians are much much faster.
///        -A 2018-9-11

template<class Scalar, class Integer>

class Filter3D {

public:
  Scalar *afH;         //!< contiguous block of memory storing the filter array
  Scalar ***aaafH;     //!< the same array which can be indexed using [i][j][k] notation
  Integer halfwidth[3]; //!< num pixels from filter center to edge in x,y,z directions
  Integer array_size[3]; //!<size of the array in x,y,z directions (2*halfwidth+1)

  /// @brief  Apply the filter to a 3D image (aaafSource[][][]).
  /// Save the results in the "aafDest" array.  (A "mask" is optional.)
  /// All arrays are 3D and assumed to be the same size.
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
  /// where: f[i] is the original image at position i (ix,iy,iz)
  ///      i = shorthand for (ix,iy,iz) = a location in the filtered image
  ///      j = shorthand for (jx,jy,jz) is summed over the entries in the filter
  ///    i-j = shorthand for (ix-jx, iy-jy, iz-jz)
  ///        g[i] is the image after the filter has been applied
  ///        h[j] is the filter
  ///        mask[i] selects the pixels we care about (usually either 0 or 1)
  ///          (If not supplied, assumed it is 1 everywhere inside, 0 outside)
  ///        Theta[i] = 1 if i is inside the image boundaries, 0 otherwise.
  ///          (In other words, we only consider voxels within the image.)
  /// @endcode
  ///
  /// @param size_source contains size of the source image (in the x,y,z directions)
  /// @param aaafSource[][][] is the source array (source image) <==> "f[i]"
  /// @param aaafDest[][][] will store the image after filtering <==> "g[i]"
  /// @param aaafMask[][][]==0 whenever we want to ignore entries in afSource[]. Optional.
  /// @param normalize  This boolean parameter = true if you want to divide g(i) by the sum of the weights considered. Optional.
  /// (useful if the sum of your filter elements, h(j), is 1, and if the sum was
  ///  not complete because some entries lie outside the mask or the boundary.)

  void Apply(Integer const size_source[3],
             Scalar const *const *const *aaafSource,
             Scalar ***aaafDest,
             Scalar const *const *const *aaafMask = NULL,
             bool normalize = false,
             ostream *pReportProgress = NULL) const
  {
    Scalar *afDenominator = NULL;
    Scalar ***aaafDenominator = NULL;
    if (normalize)
      Alloc3D(size_source, &afDenominator, &aaafDenominator);

    Apply(size_source,
          aaafSource,
          aaafDest,
          aaafMask,
          aaafDenominator,
          pReportProgress);

    if (normalize) {
      for (Integer iz=0; iz < size_source[2]; iz++)
        for (Integer iy=0; iy < size_source[1]; iy++)
          for (Integer ix=0; ix < size_source[0]; ix++)
            if (aaafDenominator[iz][iy][ix] > 0.0)
              aaafDest[iz][iy][ix] /= aaafDenominator[iz][iy][ix];
      Dealloc3D(size_source, &afDenominator, &aaafDenominator);
    }
  }

  /// @brief  Apply the filter to a 3D image (aaafSource[][][]).
  ///         This version is identical to the other version of Apply()
  ///         except that this version returns both d(i) and g(i) whenever
  ///         you supply a non-NULL afDenominator[] argument (see below).
  ///         It also does NOT normalize the result (by dividing g(i) / d(i)).
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
  /// where: f[i] is the original array of (source) data at position i
  ///      i = shorthand for (ix,iy,iz) = a location in the filtered image
  ///      j = shorthand for (jx,jy,jz) is summed over the entries in the filter
  ///    i-j = shorthand for (ix-jx, iy-jy, iz-jz)
  ///        g[i] is the data after the filter has been applied
  ///        h[j] is the filter
  ///        d[i] is the "denominator" = sum of the filter weights considered
  ///        mask[i] selects the pixels we care about (usually either 0 or 1)
  ///          (If not supplied, assumed it is 1 everywhere inside, 0 outside]
  ///        Theta[i] = 1 if i is inside the image boundaries, 0 otherwise.
  ///          (In other words, we only consider voxels within the image.)
  /// @endcode
  ///
  /// @param size_source contains size of the source image (in the x,y,z directions)
  /// @param aaafSource[][][] is the source array (source image) <==> "f[i]"
  /// @param aaafDest[][][] will store the image after filtering <==> "g[i]"
  /// @param aaafMask[][][]==0 whenever we want to ignore entries in afSource[][]. Optional.
  /// @param aaafDenominator[][][] will store d[i] if you supply a non-NULL pointer

  void Apply(Integer const size_source[3],
             Scalar const *const *const *aaafSource,
             Scalar ***aaafDest,
             Scalar const *const *const *aaafMask = NULL,
             Scalar ***aaafDenominator = NULL,
             ostream *pReportProgress = NULL) const
  {


    if (pReportProgress)
      *pReportProgress << "  progress: processing plane#" << endl;

    // The mask should be 1 everywhere we want to consider, and 0 elsewhere.
    // Multiplying the density in the tomogram by the mask removes some of 
    // the voxels from consideration later on when we do the filtering.
    // (Later, we will adjust the weight of the average we compute when we
    //  apply the filter in order to account for the voxels we deleted now.)

    for (Integer iz=0; iz<size_source[2]; iz++) {

      if (pReportProgress)
        *pReportProgress << "  " << iz+1 << " / " << size_source[2] << "\n";

      #pragma omp parallel for collapse(2)
      for (Integer iy=0; iy<size_source[1]; iy++) {

        for (Integer ix=0; ix<size_source[0]; ix++) {

          // Calculate the effect of the filter on
          // the voxel located at position ix,iy,iz

          if ((aaafMask) && (aaafMask[iz][iy][ix] == 0.0)) {
            aaafDest[iz][iy][ix] = 0.0;
            continue;
          }

          aaafDest[iz][iy][ix] =
            ApplyToVoxel(ix, iy, iz,
                         size_source,
                         aaafSource,
                         aaafMask,
                         (aaafDenominator
                          ? &(aaafDenominator[iz][iy][ix])
                          : NULL));
        }
      }
    }
  } // Apply()



  inline Filter3D(const Filter3D<Scalar, Integer>& source) {
    Init();
    Resize(source.halfwidth); // allocates and initializes afH and aaafH
    //for(Integer iz=-halfwidth[2]; iz<=halfwidth[2]; iz++)
    //  for(Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
    //    for(Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
    //      aaafH[iz][iy][ix] = source.aaafH[iz][iy][ix];
    // -- Use memcpy() instead: --
    //memcpy(afH,
    //       source.afH,
    //       array_size[0] * array_size[1] * array_size[2],
    //       *sizeof(Scalar));
    // -- Use std:copy() instead: --
    std::copy(source.afH,
              source.afH + (array_size[0] * array_size[1] * array_size[2]),
              afH);
  }


  Filter3D(Integer const set_halfwidth[3]) {
    Init();
    Resize(set_halfwidth);
  }


  Filter3D() {
    Init();
  }


  ~Filter3D() {
    Dealloc();
  }


  inline void swap(Filter3D<Scalar, Integer> &other) {
    std::swap(afH, other.afH);
    std::swap(aaafH, other.aaafH);
    std::swap(halfwidth, other.halfwidth);
    std::swap(array_size, other.array_size);
  }


  inline Filter3D<Scalar, Integer>&
    operator = (Filter3D<Scalar, Integer> source) {
    this->swap(source);
    return *this;
  }


  /// @ brief   Make sure the sum of the filter weights (in aaafH) is 1
  void Normalize() {
    Scalar total = 0.0;
    for (Integer iz=-halfwidth[2]; iz<=halfwidth[2]; iz++)
      for (Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
        for (Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
          total += aaafH[iz][iy][ix];
            
    assert(total > 0.0);
    for (Integer iz=-halfwidth[2]; iz<=halfwidth[2]; iz++)
      for (Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
        for (Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
          aaafH[iz][iy][ix] /= total;
  }

  /// @brief Calculate the (weighted) average value of the filter array aaafH
  /// @param aaafW optional weights used when calculating the averages
  /// @return the (weighted) average value of the filter array aaafH
  Scalar Average(Scalar const *const *const *aaafW=NULL) const {
    return AverageArr(array_size, aaafH, aaafW);
  }

  /// @brief Calculate the (weighted) average squared values in the filter array, aaafH
  /// @param aaafW optional weights used when calculating the averages
  /// @return the (weighted) average squared values in the filter array aaafH
  Scalar AverageSqr(Scalar const *const *const *aaafW=NULL) const {
    return _AveSqrArr(array_size, aaafH, aaafW);
  }

  /// @brief Calculate the (weighted) standard deviation of the filter array values
  /// @param aaafW optional weights used when calculating the standard deviation
  /// @return the (weighted) standard deviation of the filter array aaafH
  Scalar StdDev(Scalar const *const *const *aaafW=NULL) const {
    return StdDevArr(array_size, aaafH, aaafW);
  }

  /// @brief Calculate the (weighted) sum of the filter array values, aaafH
  /// @param aaafW optional weights used when calculating the sum
  /// @return the (weighted) sum of the filter array values
  Scalar Sum(Scalar const *const *const *aaafW=NULL) const {
    return _SumArr(array_size, aaafH, aaafW);
  }

  /// @brief Calculate the (weighted) sum of the squared filter array values
  /// @param aaafW optional weights used when calculating the sum
  /// @return the (weighted) sum of the squared filter array values
  Scalar SumSqr(Scalar const *const *const *aaafW=NULL) const {
    return _SumSqrArr(array_size, aaafH, aaafW);
  }

  /// @brief Add a number to all of the filter array values, aaafH
  /// @param offset  the number to add
  void AddScalar(Scalar offset) {
    AddScalarArr(offset, array_size, aaafH);
  }

  /// @brief multiply all of the filter array values (aaafH) by a number
  /// @param offset  the number to multiply
  void MultiplyScalar(Scalar scale) {
    MultiplyScalarArr(scale, array_size, aaafH);
  }



 private:
  /// @brief Apply a filter to the source image (aaafSource) at a particular
  ///        voxel location.  This function was not intended for public use.
  /// @param ix  the voxel's position in that 3D image
  /// @param iy  the voxel's position in that 3D image
  /// @param iz  the voxel's position in that 3D image
  /// @param size_source contains size of the source image (in the x,y,z directions)
  /// @param aaafSource[][][] is the source image
  /// @param optional. aaafMask[i][j][k]==0 for voxels we want to exclude from consideration
  /// @param pDenominator=if you want to store the sum of the weights considered, pass a pointer to a number
  /// (useful if the sum was not complete due to some voxels being masked out,
  ///  or because the filter extends beyond the boundaries of the image)
  Scalar ApplyToVoxel(Integer ix,
                       Integer iy,
                       Integer iz,
                       Integer const size_source[3],
                       Scalar const *const *const *aaafSource,
                       Scalar const *const *const *aaafMask = NULL,
                       Scalar *pDenominator = NULL) const

  {
    Scalar g = 0.0;
    Scalar denominator = 0.0;

    for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++) {
      Integer iz_jz = iz-jz;
      if ((iz_jz < 0) || (size_source[2] <= iz_jz))
        continue;

      for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {
        Integer iy_jy = iy-jy;
        if ((iy_jy < 0) || (size_source[1] <= iy_jy))
          continue;

        for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {
          Integer ix_jx = ix-jx;
          if ((ix_jx < 0) || (size_source[0] <= ix_jx))
            continue;

          Scalar filter_val = aaafH[jz][jy][jx];

          if (aaafMask) {
            Scalar mask_val = aaafMask[iz_jz][iy_jy][ix_jx];
            if (mask_val = 0.0)
              continue;
            filter_val *= mask_val;
          }
          //Note: The "filter_val" also is needed to calculate
          //      the denominator used in normalization.
          //      It is unusual to use a mask unless you intend
          //      to normalize the result later, but I don't enforce this

          Scalar delta_g = 
            filter_val * aaafSource[iz_jz][iy_jy][ix_jx];

          g += delta_g;

          if (pDenominator)
            denominator += filter_val;
        }
      }
    }

    if (pDenominator)
      *pDenominator = denominator;

    return g;
  } // ApplyToVoxel()



  /// @brief allocate space for the filter array
  void Alloc(Integer const set_halfwidth[3]) {
    //Integer array_size[3];
    for(int d=0; d < 3; d++) {
      halfwidth[d] = set_halfwidth[d];
      array_size[d] = 1 + 2*halfwidth[d];
    }
    Alloc3D(array_size, &afH, &aaafH);
    for (Integer iz = 0; iz < array_size[2]; iz++)
      for (Integer iy = 0; iy < array_size[1]; iy++)
        for (Integer ix = 0; ix < array_size[0]; ix++)
          aaafH[iz][iy][ix] = -1.0e38; //(if uninitiliazed memory read, we will know)

    //shift pointers to enable indexing from i = -halfwidth .. +halfwidth
    for (Integer iz = 0; iz < array_size[2]; iz++) {
      for (Integer iy = 0; iy < array_size[1]; iy++) {
        aaafH[iz][iy] += halfwidth[0];
      }
      aaafH[iz] += halfwidth[1];
    }
    aaafH += halfwidth[2];
  }


  /// @brief allocate space used by the filter array
  void Dealloc() {
    if (! aaafH) {
      assert(! afH);
      Init();
      return;
    }
    //Integer array_size[3];
    for(int d=0; d < 3; d++) {
      array_size[d] = 1 + 2*halfwidth[d];
    }
    //shift pointers back to normal
    for (Integer iz = -halfwidth[2]; iz <= halfwidth[2]; iz++) {
      for (Integer iy = -halfwidth[1]; iy <= halfwidth[1]; iy++) {
        aaafH[iz][iy] -= halfwidth[0];
      }
      aaafH[iz] -= halfwidth[1];
    }
    aaafH -= halfwidth[2];
    //then deallocate
    Dealloc3D(array_size, &afH, &aaafH);
    for(int d=0; d < 3; d++) {
      array_size[d] = -1;
      halfwidth[d] = -1;
    }
  }


  void Resize(Integer const set_halfwidth[3]) {
    Dealloc();
    Alloc(set_halfwidth);
  }


  void Init() {
    halfwidth[0] = -1;
    halfwidth[1] = -1;
    halfwidth[2] = -1;
    array_size[0] = -1;
    array_size[1] = -1;
    array_size[2] = -1;
    afH = NULL;
    aaafH = NULL;
  }
  



  #ifndef DISABLE_TEMPLATE_MATCHING
private:
  /// @brief  Compute the sum weighted difference between the local entries
  ///         in an array centered around position ix,iy,iz, 
  ///         and another array we call the "template" (which is assumed to
  ///         be stored in the "aaafH" filter array member for this class).
  ///         It is assumed that the weighted average of both the
  ///          source and template arrays is 0 before calling this function.
  /// @return the weighted sum of the error (raised to the power of
  ///         err_exponent, which is 2 by default) summed over all entries
  ///         in the template.  Weights are typically chosen so that they
  ///         decay to 0 near the boundaries of the aaafH array.
  ///         This function is not yet intended for public use.
  Scalar
  _TemplateError(Integer ix, //!< voxel's position in the x,y,z directions
                 Integer iy, //!< voxel's position in the x,y,z directions
                 Integer iz, //!< voxel's position in the x,y,z directions
                 Integer const size_source[3], //!< size of the source array
                 Scalar const *const *const *aaafSource,   //!< the array 
                 Scalar scale, //!< how much to scale template intensities
                 Scalar const *const *const *aaafW, //!< weights used in all summations
                 Scalar err_exponent=2, //!< exponent used for calculating error
                 Scalar const *const *const *aaafMask = NULL, //!< optional array storing voxels we should ignore (if 0)
                 ostream *pReportProgress = NULL //!< optional ostream for printing out progress of the calculation
                 ) const 
  {
      Scalar g = 0.0;
      Scalar denominator = 0.0;

      for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++) {
      Integer iz_jz = iz-jz;
      if ((iz_jz < 0) || (size_source[2] <= iz_jz))
        continue;

      for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {
        Integer iy_jy = iy-jy;
        if ((iy_jy < 0) || (size_source[1] <= iy_jy))
          continue;

        for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {
          Integer ix_jx = ix-jx;
          if ((ix_jx < 0) || (size_source[0] <= ix_jx))
            continue;

          Scalar delta_g = 
            (scale * aaafH[jz][jy][jx]
             -
             aaafSource[iz_jz][iy_jy][ix_jx]);

          if (err_exponent == 2.0)
            delta_g *= delta_g;
          else
            delta_g = pow(delta_g, err_exponent);

          //if (! precompute_mask_times_source)
          delta_g *= aaafW[jz][jy][jx];


          if (aaafMask)     //<--rarely used. I may delete "aaafMask" later
            delta_g *= aaafMask[iz_jz][iy_jy][ix_jx];

          g += delta_g;

          //if (normalize) {
          if (aaafMask)     //<--rarely used. I may delete "aaafMask" later
            denominator +=
              aaafMask[iz_jz][iy_jy][ix_jx]
              *
              aaafW[jz][jy][jx];

          else
            //denominator += 1.0;
            denominator += aaafW[jz][jy][jx];

          //} //if (normalize)

        }
      }
    }

    //if (normalize) {
    if (denominator > 0.0)
      g /= denominator;
    else
      g = 1.0;

    return g;
  } // _TemplateError()





  /// @brief Calculate the TemplateError() everywhere in the source array
  ///        (aaafSource).  Store the result in aaafDest.
  ///         This function was not intended for public use.

  void
  _ScanTemplateError(Integer const size_source[3], //!< size of the source array
                     Scalar const *const *const *aaafSource, //!< source array
                     Scalar ***aaafDest, //!< store the template error results here
                     Scalar ***aaafC, //!< how much to scale template intensities
                     Scalar const *const *const *aaafW, //!< weights used in all summations
                     Scalar err_exponent=2, //!< exponent used for calculating error
                     Scalar const *const *const *aaafMask = NULL, //!< optional: indicate which entries should be ignored
                     ostream *pReportProgress = NULL //!< optional: print out the progress of the calculation
                     ) const
  {
    if (pReportProgress)
      *pReportProgress << "  progress: processing plane#" << endl;

    for (Integer iz=0; iz<size_source[2]; iz++) {

      if (pReportProgress)
        *pReportProgress << "  " << iz+1 << " / " << size_source[2] << "\n";

      for (Integer iy=0; iy<size_source[1]; iy++) {

        for (Integer ix=0; ix<size_source[0]; ix++) {

          // Calculate the effect of the filter on
          // the voxel located at position ix,iy,iz

          if ((aaafMask) && (aaafMask[iz][iy][ix] == 0.0)) {
            aaafDest[iz][iy][ix] = 1.0;
            continue;
          }
          
          aaafDest[iz][iy][ix] =
            TemplateError(ix, iy, iz,
                          size_source,
                          aaafSource,
                          aaafC[iz][iy][ix],
                          aaafW,
                          err_exponent,
                          NULL); //aaafMask);
                          //normalize);
        }
      }
    }
  } //_ScanTemplateError()

  #endif //#ifndef DISABLE_TEMPLATE_MATCHING

}; // class Filter3D






/// @brief Create a 3D filter and fill it with a "generalized Gaussian" function
///   @verbatim h(x,y,z) = A*exp(-r^m)  @endverbatim
/// where  @verbatim r  = sqrt((x/σ_x)^2 + (y/σ_y)^2 + (z/σ_z)^2) @endverbatim
///   and  "A" is determined by normalization of the discrete sum
/// @note "A" is equal to the value stored in the middle of the array,
///       The caller can determine what "A" is by looking at this value.

template<class Scalar>
Filter3D<Scalar, int>
GenFilterGenGauss3D(Scalar width[3],    //!< "σ_x", "σ_y", "σ_z" parameters
                    Scalar m_exp,       //!< "m" exponent parameter
                    int truncate_halfwidth[3], //!< size of filter window
                    Scalar *pA=NULL,    //!< optional:report A coeff to user
                    ostream *pReportEquation=NULL//!< optional:report equation used to the user
                    )
{
  Scalar window_threshold = 1.0;
  for (int d=0; d<3; d++) {
    Scalar h = ((width[d]>0)
                 ? exp(-pow(truncate_halfwidth[d]/width[d], m_exp))
                 : 1.0);
    if (h < window_threshold)
      window_threshold = h;
  }
  Filter3D<Scalar, int> filter(truncate_halfwidth);
  Scalar total = 0;
  for (int iz=-filter.halfwidth[2]; iz<=filter.halfwidth[2]; iz++) {
    for (int iy=-filter.halfwidth[1]; iy<=filter.halfwidth[1]; iy++) {
      for (int ix=-filter.halfwidth[0]; ix<=filter.halfwidth[0]; ix++) {
        Scalar r = sqrt(SQR(ix/width[0])+SQR(iy/width[1])+SQR(iz/width[2]));
        Scalar h = ((r>0) ? exp(-pow(r, m_exp)) : 1.0);
        if (ABS(h) < window_threshold)
          h = 0.0; //This eliminates corner entries which fall below threshold
                   //(and eliminates anisotropic artifacts due to these corners)
                   //There's no reason to keep any entries less than min value.
        filter.aaafH[iz][iy][ix] = h;
                    
        total += h;
      }
    }
  }
  // normalize:
  for (int iz=-filter.halfwidth[2]; iz<=filter.halfwidth[2]; iz++) {
    for (int iy=-filter.halfwidth[1]; iy<=filter.halfwidth[1]; iy++) {
      for (int ix=-filter.halfwidth[0]; ix<=filter.halfwidth[0]; ix++) {

        filter.aaafH[iz][iy][ix] /= total;

        //FOR DEBUGGING REMOVE EVENTUALLY:
        if (pReportEquation)
          *pReportEquation << "GenGauss3D:" //<< window_threshold
                           << " aaafH["<<iz<<"]["<<iy<<"]["<<ix<<"] = "
                           << filter.aaafH[iz][iy][ix] << endl;
      }
    }
  }

  // The coefficient in front of the Gaussian ("A")
  // equals the height of its central peak,
  // which is located in the middle of the array
  Scalar A = filter.aaafH[0][0][0];


  if (pA) {
    *pA = A;
  }

  if (pReportEquation) {
    *pReportEquation << " Filter Used:\n"
      " h(x,y,z)   = A*exp(-((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2)^(m/2))\n"
      "  ... where      A = " << A << "\n"
      "                 m = " << m_exp << "\n" 
      "   (a_x, a_y, a_z) = "
                     << "(" << width[0]
                     << " " << width[1]
                     << " " << width[2] << ")\n";
    *pReportEquation << " You can plot a slice of this function\n"
                     << "     in the X direction using:\n"
      " draw_filter_1D.py -gauss " << A << " " << width[0] << endl;
    *pReportEquation << " and in the Y direction using:\n"
      " draw_filter_1D.py -gauss " << A << " " << width[1] << endl;
    *pReportEquation << " and in the Z direction using:\n"
      " draw_filter_1D.py -gauss " << A << " " << width[2] << endl;
    *pReportEquation << " You can plot a slice of this function\n"
                     << "     in the X direction using:\n"
      " draw_filter_1D.py -ggauss " << A
                     << " " << width[0]
                     << " " << m_exp << endl;
    *pReportEquation << " and in the Y direction using:\n"
      " draw_filter_1D.py -ggauss " << A
                     << " " << width[1]
                     << " " << m_exp << endl;
    *pReportEquation << " and in the Z direction using:\n"
      " draw_filter_1D.py -ggauss " << A
                     << " " << width[2]
                     << " " << m_exp << endl;
  }
  return filter;
} //GenFilterGenGauss3D(width, m_exp, truncate_halfwidth)







/// @brief Create a 3D filter and fill it with a "generalized Gaussian" function
///   @verbatim h(x,y,z) = A*exp(-r^m)  @endverbatim
/// where  @verbatim r  = sqrt((x/σ_x)^2 + (y/σ_y)^2 + (z/σ_z)^2) @endverbatim
///   and  "A" is determined by normalization of the discrete sum
/// @note "A" is equal to the value stored in the middle of the array (aaafH),
///       The caller can determine what "A" is by looking at aaafH there

template<class Scalar>
Filter3D<Scalar, int>
GenFilterGenGauss3D(Scalar width[3],    //!< "σ_x", "σ_y", "σ_z" parameters
                    Scalar m_exp,       //!< "m" parameter in formula
                    Scalar filter_cutoff_ratio=2.5, //!< how many sigma (σ) before truncating?
                    Scalar *pA=NULL,    //!< optional:report A coeff to user
                    ostream *pReportEquation = NULL //!< optional:report equation used to the user
                    )
{
  // choose the width of the filter window based on the filter_cutoff_ratio
  int truncate_halfwidth[3];
  int ix = 0;

  for (int d=0; d<3; d++) {
    truncate_halfwidth[d] = floor(width[d]*filter_cutoff_ratio);
  }
  return GenFilterGenGauss3D(width,
                             m_exp,
                             truncate_halfwidth,
                             pA,
                             pReportEquation);
} //GenFilterGenGauss3D(width, m_exp, filter_cutoff_ratio)






/// @brief  Create a 3D filter and fill it with a difference of (generalized) Gaussians 
///
/// This version requires that the caller has already created individual
/// filters for the two gaussians.
/// All this function does is subtract one filter from the other (and rescale).
/// This function was not intended for public use.

template<class Scalar>
static
Filter3D<Scalar, int> 
_GenFilterDogg3D(Scalar width_a[3],  //!< "a" parameter in formula
                 Scalar width_b[3],  //!< "b" parameter in formula
                 Scalar m_exp,  //!< "m" parameter in formula
                 Scalar n_exp,  //!< "n" parameter in formula
                 Filter3D<Scalar, int>& filter_A, //!< filters for the two
                 Filter3D<Scalar, int>& filter_B, //!< gaussians
                 Scalar *pA=NULL, //!< optional:report A,B coeffs to user
                 Scalar *pB=NULL, //!< optional:report A,B coeffs to user
                 ostream *pReportEquation = NULL //!< optional: report equation to the user
                 )
{
  Scalar A, B;
  //A, B = height of the central peak
  A = filter_A.aaafH[0][0][0];
  B = filter_B.aaafH[0][0][0];


  // The "difference of gaussians" filter is the difference between
  // these two (generalized) gaussian filters.
  int halfwidth[3];
  halfwidth[0] = MAX(filter_A.halfwidth[0], filter_B.halfwidth[0]);
  halfwidth[1] = MAX(filter_A.halfwidth[1], filter_B.halfwidth[1]);
  halfwidth[2] = MAX(filter_A.halfwidth[2], filter_B.halfwidth[2]);
  Filter3D<Scalar, int> filter(halfwidth);

  //FOR DEBUGGING REMOVE EVENTUALLY
  if (pReportEquation)
    *pReportEquation << "Array of 3D filter entries:" << endl;


  for (int iz=-halfwidth[2]; iz<=halfwidth[2]; iz++) {
    for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
      for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {

        filter.aaafH[iz][iy][ix] = 0.0;

        // The two filters may have different widths, so we have to check
        // that ix,iy and iz lie within the domain of these two filters before
        // adding or subtracting their values from the final GDOG filter.
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
    *pReportEquation << " and in the Y direction using:\n"
      " draw_filter_1D.py -dogg " << A << " " << B
                     << " " << width_a[1] << " " << width_b[1]
                     << " " << m_exp << " " << n_exp << endl;
    *pReportEquation << " and in the Z direction using:\n"
      " draw_filter_1D.py -dogg " << A << " " << B
                     << " " << width_a[2] << " " << width_b[2]
                     << " " << m_exp << " " << n_exp << endl;
  } //if (pReportEquation)
  
  return filter;
} //_GenFilterDogg3D()



/// @brief Create a 3D filter and fill it with a difference of (generalized) Gaussians
///   @verbatim h(x,y,z) = A*exp(-(r/a)^m) - B*exp(-(r/b)^n)  @endverbatim
/// where  @verbatim r = sqrt(x^2 + y^2 + z^2) @endverbatim
///   and "A" and "B" are determined by normalization of each term independently
/// (It's not clear whether this kind of filter is useful when m>2 or n>2)
template<class Scalar>
Filter3D<Scalar, int> 
GenFilterDogg3D(Scalar width_a[3],   //!< "a" parameter in formula
                Scalar width_b[3],   //!< "b" parameter in formula
                Scalar m_exp,        //!< "m" parameter in formula
                Scalar n_exp,        //!< "n" parameter in formula
                int halfwidth[3],     //!< the width of the filter
                Scalar *pA=NULL,     //!< optional:report A,B coeffs to user
                Scalar *pB=NULL,     //!< optional:report A,B coeffs to user
                ostream *pReportEquation = NULL //!< optional: print params used?
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

  return _GenFilterDogg3D(width_a,
                          width_b,
                          m_exp,
                          n_exp,
                          filter_A, filter_B,
                          pA,
                          pB,
                          pReportEquation);
} //GenFilterDogg3D(...halfwidth...)






// @brief ApplySeparable3D applies a separable filter on a 3D array.
//        It assumes separate Filter1D objects have already been created 
//        which will blur the image in each direction (x,y,z).
//        This function supports masks (which exclude voxels from consideration)
//        One major feature of this function is its ability to efficiently
//        normalize the resulting filtered image in the presence of a mask
//        (as well as near the image boundaries). Consequently, the image does
//        not fade to black near the boundaries of the image or the mask.
//        This function was not intended for public use.
template<class Scalar>
static
Scalar
ApplySeparable3D(int const image_size[3], 
                 Scalar const *const *const *aaafSource,
                 Scalar ***aaafDest,
                 Scalar const *const *const *aaafMask,
                 Filter1D<Scalar, int> aFilter[3], //preallocated 1D filters
                 bool normalize = true,
                 ostream *pReportProgress = NULL)
{
  assert(aaafSource);
  assert(aaafDest);

  // This is a "separable" filter.
  // The filter is fast because we apply the filter sequentially in the 
  // in the X, Y, Z directions (instead of applying the filter simultaneously
  // in all 3 directions, requiring a sum over all voxels with a given radius).

  // Initially copy aaafSource into aaafDest
  // (We don't want to have to allocate temporary array to 
  //  store the result of each successive filter operation. 
  //  Instead just store the most recent filter operation in aaafDest,
  //  and perform each operation on whatever's currently in aaafDest.)
  for (int iz = 0; iz < image_size[2]; iz++)
    for (int iy = 0; iy < image_size[1]; iy++)
      for (int ix = 0; ix < image_size[0]; ix++)
        aaafDest[iz][iy][ix] = aaafSource[iz][iy][ix];

  // Some filters, (such as Gaussian filters) are normalizable.
  // That means these filters are weighted averaging of nearby voxels
  // (with nonnegative weights).
  // However sometimes these nearby voxels are unavailable because they either
  // lie outside the boundaries of the image, or they lie outside the mask.
  // In that case, we can "normalize" the resulting filtered image by dividing
  // the weighted average brightness of the voxels near a particular location,
  // ...by the sum of the weights which were used to calculate that average.
  // (at that location).  The sum of those weights are called the "denominator".
  // Create an array to store the denominator.
  // First create the 3D version of the denominator array:
  Scalar ***aaafDenom = NULL;
  Scalar *afDenom = NULL;

  if (normalize) {
    if (aaafMask) {
      Alloc3D(image_size, &afDenom, &aaafDenom);
      for (int iz = 0; iz < image_size[2]; iz++)
        for (int iy = 0; iy < image_size[1]; iy++)
          for (int ix = 0; ix < image_size[0]; ix++)
            aaafDenom[iz][iy][ix] = 1.0; //(default value)
    }
  } // if (normalize) 

  int d; //direction where we are applying the filter (x<==>0, y<==>1, z<==>2)

  // First, apply the filter in the Z direction (d=2):
  d = 2;
  if (pReportProgress)
    *pReportProgress << "  progress: Applying Z filter"// Processing XY plane#"
                     << endl;
  
  #pragma omp parallel
  {
    // Don't want to create unnecessary 3D arrays to store the entire image 
    // after each successive filtering (in the x,y,z directions).
    // Instead use temporary arrays which store the image, (mask, denom, etc)
    // along the direction that we intend to apply the filter at this step.
    // Then use a simple 1-D filter on that array and copy the results back.
    Scalar *afDest_tmp   = new Scalar [image_size[d]];
    Scalar *afSource_tmp = new Scalar [image_size[d]];
    Scalar *afMask_tmp   = NULL;
    if (aaafMask)
      afMask_tmp = new Scalar [image_size[d]];
    Scalar *afDenom_tmp = NULL;
    if (normalize && aaafMask)
      afDenom_tmp = new Scalar [image_size[d]];

    #pragma omp for collapse(2)
    for (int iy = 0; iy < image_size[1]; iy++) {
      for (int ix = 0; ix < image_size[0]; ix++) {

        //CONTINUEHERE: std::bad_alloc() exception here using OpenMP 2018-10-18!


        // Copy the data we need to the temporary arrays
        for (int iz = 0; iz < image_size[2]; iz++) {
          afSource_tmp[iz] = aaafDest[iz][iy][ix];  //copy from prev aaafDest
          if (aaafMask)
            afMask_tmp[iz] = aaafMask[iz][iy][ix];
        }

        // Apply the filter to the 1-D temporary arrays which contain the source
        // image data (afSource_tmp), as well as the mask data (afMask_tmp,
        // if applicable).  The following version of the "Apply()" function 
        // can apply the filter to both of these signals.  In this case, the 
        // computational overhead required for filtering both signals is only 
        // slightly higher than the cost required for filtering only one of them
        // (Later this won't be true and we will have to filter them separately)
        aFilter[d].Apply(image_size[d],
                         afSource_tmp,
                         afDest_tmp,  //<-store filtered result here
                         afMask_tmp,
                         afDenom_tmp);//<-store sum of weights considered here

        // copy the results from the temporary filters back into the 3D arrays
        for (int iz = 0; iz < image_size[2]; iz++) {
          aaafDest[iz][iy][ix] = afDest_tmp[iz];
          if (normalize && aaafMask)
            aaafDenom[iz][iy][ix] = afDenom_tmp[iz]; //copy back into aaafDenom
            // (Note: if aaafMask==NULL then we normalize using a faster method)
        }
      } //for (int ix = 0; ix < image_size[0]; ix++)
    } //for (int iy = 0; iy < image_size[1]; iy++)

    // delete the temporary arrays
    delete [] afSource_tmp;
    delete [] afDest_tmp;
    if (afMask_tmp)
      delete [] afMask_tmp;
    if (afDenom_tmp)
      delete [] afDenom_tmp;
  } //#pragma omp parallel private(afDest_tmp, ...)



  // Then apply the filter in the Y direction (d=1):
  d = 1;
  if (pReportProgress)
    *pReportProgress << "  progress: Applying Y filter"// Processing XZ plane#"
                     << endl;

  #pragma omp parallel
  {
    // Don't want to create unnecessary 3D arrays to store the entire image 
    // after each successive filtering (in the x,y,z directions).
    // Instead use temporary arrays which store the image, (mask, denom, etc)
    // along the direction that we intend to apply the filter at this step.
    // Then use a simple 1-D filter on that array and copy the results back.
    Scalar *afDest_tmp   = new Scalar [image_size[d]];
    Scalar *afSource_tmp = new Scalar [image_size[d]];
    //Scalar *afMask_tmp   = NULL;
    //if (aaafMask)
    //  afMask_tmp = new Scalar [image_size[d]];
    Scalar *afDenom_src_tmp = NULL;
    Scalar *afDenom_tmp = NULL;
    if (normalize && aaafMask) {
      afDenom_src_tmp = new Scalar [image_size[d]];
      afDenom_tmp     = new Scalar [image_size[d]];
    }

    #pragma omp for collapse(2)
    for (int iz = 0; iz < image_size[2]; iz++) {
      for (int ix = 0; ix < image_size[0]; ix++) {

        // copy the data we need to the temporary arrays
        for (int iy = 0; iy < image_size[1]; iy++) {
          afSource_tmp[iy] = aaafDest[iz][iy][ix];  //copy from prev aaafDest
          //if (aaafMask)
          //  afMask_tmp[iy] = aaafMask[iz][iy][ix];
          if (normalize && aaafMask)
            afDenom_src_tmp[iy] = aaafDenom[iz][iy][ix];
        }

        // At this point, the convolution of the 1-D filter along
        // the Z direction on the afSource[][][] array is stored in both 
        // aaafDeest[][][] and also afSource_tmp[].  Now we want to
        // apply the 1-D filter along the Y direction to that data.
        aFilter[d].Apply(image_size[d],
                         afSource_tmp,
                         afDest_tmp); //<-store filtered result here

        if (normalize && aaafMask)
          // At this point, the convolution of the 1-D filter along
          // the Z direction on the aaafMask[][][] array is stored in both 
          // aaafDenom[][][] and also afDenom_src_tmp[].  Now we want to
          // apply the 1-D filter along the Y direction to that data
          aFilter[d].Apply(image_size[d],
                           afDenom_src_tmp, //<-weights so far (summed along z)
                           afDenom_tmp);//<-store sum of weights considered here

        // copy the results from the temporary filters back into the 3D arrays
        for (int iy = 0; iy < image_size[1]; iy++) {
          aaafDest[iz][iy][ix] = afDest_tmp[iy];
          if (normalize && aaafMask)
            //copy the weights from afDenom_tmp[] into aaafDenom[][][]
            // (Note: if aaafMask==NULL then we normalize using a faster method)
            aaafDenom[iz][iy][ix] = afDenom_tmp[iy];
        }
      } //for (int ix = 0; ix < image_size[0]; ix++)
    } //for (int iz = 0; iz < image_size[2]; iz++)

    // delete the temporary arrays
    delete [] afSource_tmp;
    delete [] afDest_tmp;
    //if (afMask_tmp)
    //  delete [] afMask_tmp;
    if (afDenom_tmp)
      delete [] afDenom_tmp;
    if (afDenom_src_tmp)
      delete [] afDenom_src_tmp;
  } //#pragma omp parallel private(afDest_tmp, ...)




  // Finally apply the filter in the X direction (d=0):
  d = 0;
  if (pReportProgress)
    *pReportProgress << "  progress: Applying X filter"// Processing YZ plane#"
                     << endl;
  
  #pragma omp parallel
  {
    // Don't want to create unnecessary 3D arrays to store the entire image 
    // after each successive filtering (in the x,y,z directions).
    // Instead use temporary arrays which store the image, (mask, denom, etc)
    // along the direction that we intend to apply the filter at this step.
    // Then use a simple 1-D filter on that array and copy the results back.
    Scalar *afDest_tmp   = new Scalar [image_size[d]];
    Scalar *afSource_tmp = new Scalar [image_size[d]];
    //Scalar *afMask_tmp   = NULL;
    //if (aaafMask)
    //  afMask_tmp = new Scalar [image_size[d]];
    Scalar *afDenom_src_tmp = NULL;
    Scalar *afDenom_tmp = NULL;
    if (normalize && aaafMask) {
      afDenom_src_tmp = new Scalar [image_size[d]];
      afDenom_tmp     = new Scalar [image_size[d]];
    }

    #pragma omp for collapse(2)
    for (int iz = 0; iz < image_size[2]; iz++) {
      for (int iy = 0; iy < image_size[1]; iy++) {

        // copy the data we need to the temporary arrays
        for (int ix = 0; ix < image_size[0]; ix++) {
          afSource_tmp[ix] = aaafDest[iz][iy][ix];  //copy from prev aaafDest
          //if (aaafMask)
          //  afMask_tmp[ix] = aaafMask[iz][iy][ix];
          if (normalize && aaafMask)
            afDenom_src_tmp[ix] = aaafDenom[iz][iy][ix];
        }

        // At this point, the convolution of the 1-D filter along both
        // the Y and Z directions on the afSource[][][] array is stored in both 
        // aaafDeest[][][] and also afSource_tmp[].  Now we want to
        // apply the 1-D filter along the X direction to that data.
        aFilter[d].Apply(image_size[d],
                         afSource_tmp,
                         afDest_tmp); //<-store filtered result here

        if (normalize && aaafMask)
          // At this point, the convolution of the 1-D filter along the
          // Y and Z directions on the aaafMask[][][] array is stored in both 
          // aaafDenom[][][] and also afDenom_src_tmp[].  Now we want to
          // apply the 1-D filter along the X direction to that data
          aFilter[d].Apply(image_size[d],
                           afDenom_src_tmp, //<-weights so far(summed along y,z)
                           afDenom_tmp);//<-store sum of weights considered here

        // copy the results from the temporary filters back into the 3D arrays
        for (int ix = 0; ix < image_size[0]; ix++) {
          aaafDest[iz][iy][ix] = afDest_tmp[ix];
          if (normalize && aaafMask)
            //copy the weights from afDenom_tmp[] into aaafDenom[][][]
            // (Note: if aaafMask==NULL then we normalize using a faster method)
            aaafDenom[iz][iy][ix] = afDenom_tmp[ix];
        }
      } //for (int iy = 0; iy < image_size[1]; iy++)
    } //for (int iz = 0; iz < image_size[2]; iz++)

    // delete the temporary arrays
    delete [] afSource_tmp;
    delete [] afDest_tmp;
    //if (afMask_tmp)
    //  delete [] afMask_tmp;
    if (afDenom_tmp)
      delete [] afDenom_tmp;
    if (afDenom_src_tmp)
      delete [] afDenom_src_tmp;
  } //#pragma omp parallel private(afDest_tmp, ...)



  if (normalize) {
    if (aaafMask) {
      assert(aaafDenom);
      for (int iz = 0; iz < image_size[2]; iz++)
        for (int iy = 0; iy < image_size[1]; iy++)
          for (int ix = 0; ix < image_size[0]; ix++)
            if (aaafDenom[iz][iy][ix] > 0.0)
              aaafDest[iz][iy][ix] /= aaafDenom[iz][iy][ix];
      // Cleanup
      Dealloc3D(image_size, &afDenom, &aaafDenom);
    }
    else {
      // Optimization when no mask is supplied:
      // If there is no mask, but the user wants the result to be normalized,
      // then we convolve the filter with the rectangular box.  This is cheaper
      // because the convolution of a separable filter with a rectangular box 
      // shaped function is the product of the convolution with three 1-D 
      // functions which are 1 from 0..image_size[d], and 0 everywhere else.
      Scalar *aafDenom_precomputed[3];
      for (int d=0; d<3; d++) {
        Scalar *afAllOnes = new Scalar [image_size[d]];
        aafDenom_precomputed[d] = new Scalar [image_size[d]];
        for (int i=0; i < image_size[d]; i++)
          afAllOnes[i] = 1.0;
        aFilter[d].Apply(image_size[d], afAllOnes, aafDenom_precomputed[d]);
        delete [] afAllOnes;
      }
      for (int iz = 0; iz < image_size[2]; iz++) {
        for (int iy = 0; iy < image_size[1]; iy++) {
          for (int ix = 0; ix < image_size[0]; ix++) {
            Scalar denominator = (aafDenom_precomputed[0][ix] *
                                  aafDenom_precomputed[1][iy] *
                                  aafDenom_precomputed[2][iz]);
            aaafDest[iz][iy][ix] /= denominator;
          }
        }
      }
      // delete the array we created for storing the precomputed denominator:
      for (int d=0; d<3; d++)
        delete [] aafDenom_precomputed[d];
    } // else clause for "if (aaafMask)"
  } // if (normalize)


  // Estimate the peak of the filter.
  //    (Currently, below, I assume the peak is located at the center,
  //     of the filter window, but perhaps I should relax that assumption.)
  // Here's an example of a Gaussian filter:
  //           filter(x) =  A  * exp(-(1/2)*((x/σ_x)^2 + (y/σ_y)^2 + (z/σ_z)^2))
  //                     = A_x * exp(-(1/2)*(x/σ_x)^2) *
  //                       A_y * exp(-(1/2)*(y/σ_y)^2) *
  //                       A_z * exp(-(1/2)*(z/σ_z)^2)
  // The 3D Gaussian coefficient ("A") is the product of the
  // 1-D Gaussian coefficients in the X,Y,Z directions (A_x * A_y * A_z).
  // Those coefficients happen to equal the value of the corresponding 1-D
  // Gaussian evaluated at its peak, which is stored in the central entry

  Scalar A_coeff = (aFilter[0].afH[0] *
                    aFilter[1].afH[0] *
                    aFilter[2].afH[0]);

  return A_coeff;

} //ApplySeparable3D(aFilter)




/// @brief Apply a Gaussian filter (blur) to an image
/// @code h(x,y,z)=A*exp(-0.5*((x/σ_x)^2+(y/σ_y)^2+(z/σ_z)^2) @endcode
/// In this version, the user manually specifies the filter window width.
/// Voxels outside the boundary of the image or outside the mask are not
/// considered during the averaging (blurring) process, and the result after
/// blurring is weighted accordingly (if normalize=true, which it is default).
///    An explanation of normalization:
/// Gaussian filters are a form of weighted averaging of nearby voxels.
/// However sometimes these nearby voxels are unavailable because they either
/// lie outside the boundaries of the image, or they lie outside the mask
///   (ie. in regions where 0 <= aaafMask[iz][iy][ix] < 1.0).
/// In that case, we can "normalize" the resulting filtered image by dividing
/// this weighted average brightness of the voxels near a particular location,
/// ...by the sum of the weights which were used to calculate that average.
/// (at that location).
/// This function was not intended for public use.
///
/// @returns the "A" coefficient (determined by normalization)

template<class Scalar>
static Scalar
ApplyGauss3D(int const image_size[3], //!< image size in x,y,z directions
             Scalar const *const *const *aaafSource,   //!< source image (3D array)
             Scalar ***aaafDest,     //!< filtered (blurred) image stored here
             Scalar const *const *const *aaafMask,     //!< ignore voxels if aaafMask[i][j][k]==0
             Scalar const sigma[3],  //!< Gaussian sigma parameters σ_x,σ_y,σ_z
             int const truncate_halfwidth[3], //!< the filter window width
             bool normalize = true,           //!< normalize the average?
             ostream *pReportProgress = NULL  //!< print progress to the user?
             )
{
  assert(aaafSource);
  assert(aaafDest);
  //assert(aaafMask);
  //Allocate filters in all 3 directions.  (Later apply them sequentially.)
  Filter1D<Scalar, int> aFilter[3];
  for (int d=0; d < 3; d++)
    aFilter[d] = GenFilterGauss1D(sigma[d], truncate_halfwidth[d]);

  return ApplySeparable3D(image_size, 
                          aaafSource,
                          aaafDest,
                          aaafMask,
                          aFilter,
                          normalize,
                          pReportProgress);

} //ApplyGauss3D(sigma, truncate_halfwidth, ...)



/// @brief Apply a Gaussian filter (blur) to an image
/// @code h(x,y,z)=A*exp(-0.5*(x^2+y^2+z^2)/σ^2) @endcode
/// In this version, the user manually specifies the filter window width.
/// Voxels outside the boundary of the image or outside the mask are not
/// considered during the averaging (blurring) process, and the result after
/// blurring is weighted accordingly (if normalize=true, which it is default).
///    An explanation of normalization:
/// Gaussian filters are a form of weighted averaging of nearby voxels.
/// However sometimes these nearby voxels are unavailable because they either
/// lie outside the boundaries of the image, or they lie outside the mask
///   (ie. in regions where 0 <= aaafMask[iz][iy][ix] < 1.0).
/// In that case, we can "normalize" the resulting filtered image by dividing
/// this weighted average brightness of the voxels near a particular location,
/// ...by the sum of the weights which were used to calculate that average.
/// (at that location).
///
/// @returns the "A" coefficient (determined by normalization)
template<class Scalar>
static Scalar
ApplyGauss3D(int const image_size[3], //!< image size in x,y,z directions
             Scalar const *const *const *aaafSource,   //!< source image (3D array)
             Scalar ***aaafDest,     //!< filtered (blurred) image stored here
             Scalar const *const *const *aaafMask, //!< ignore voxels if aaafMask[i][j][k]==0
             Scalar sigma,                   //!< Gaussian sigma parameter σ
             int     truncate_halfwidth,      //!< the filter window width

             bool normalize = true,           //!< normalize the average?
             ostream *pReportProgress = NULL  //!< print progress to the user?
             )
{
  Scalar afSigma[3] = {sigma, sigma, sigma};
  int afTruncateHalfwidth[3] = {truncate_halfwidth,
                                truncate_halfwidth,
                                truncate_halfwidth};
  ApplyGauss3D(image_size,
               aaafSource,
               aaafDest,
               aaafMask,
               afSigma,
               afTruncateHalfwidth,
               normalize,
               pReportProgress);
}





/// @brief Apply a Gaussian filter (blur) to a 3D image
///
/// @code h(x,y,z)=A*exp(-0.5*((x/σ_x)^2+(y/σ_y)^2+(z/σ_z)^2) @endcode
/// The constant "A" is determined by normalization.
/// In this version, the user specifies the filter window-width in units of σ
/// Voxels outside the boundary of the image or outside the mask are not
/// considered during the averaging (blurring) process, and the resulting
/// after blurring is weighted accordingly.  (normalize=true by default)
/// This function was not intended for public use.
///
/// @returns the "A" coefficient (determined by normalization)

template<class Scalar>
static Scalar
ApplyGauss3D(int const image_size[3], //!< image size in x,y,z directions
             Scalar const *const *const *aaafSource,   //!< source image (3D array)
             Scalar ***aaafDest,     //!< filtered (blurred) image stored here
             Scalar const *const *const *aaafMask,     //!< ignore voxels if aaafMask[i][j][k]==0
             Scalar const sigma[3],  //!< Gaussian sigma parameters σ_x,σ_y,σ_z
             Scalar truncate_ratio = 2.5,  //!< how many sigma before truncating?
             bool normalize = true,           //!< normalize the average?
             ostream *pReportProgress = NULL  //!< print progress to the user?
             )
{
  Scalar A;
  int truncate_halfwidth[3];
  if (truncate_ratio > 0)
    for (int d=0; d < 3; d++)
      truncate_halfwidth[d] = floor(sigma[d] * truncate_ratio);

  A = ApplyGauss3D(image_size, 
                   aaafSource,
                   aaafDest,
                   aaafMask,
                   sigma,
                   truncate_halfwidth,
                   normalize,
                   pReportProgress);
  return A;
} // ApplyGauss3D(..., truncate_ratio, ...)



/// @brief Apply a Gaussian filter (blur) to a 3D image
///
/// @code h(x,y,z)=A*exp(-0.5*((x^2+y^2+z^2)/σ^2) @endcode
/// The constant "A" is determined by normalization.
/// In this version, the user specifies the filter window-width in units of σ
/// Voxels outside the boundary of the image or outside the mask are not
/// considered during the averaging (blurring) process, and the resulting
/// after blurring is weighted accordingly.  (normalize=true by default)
///
/// @returns the "A" coefficient (determined by normalization)

template<class Scalar>
Scalar
ApplyGauss3D(int const image_size[3], //!< image size in x,y,z directions
             Scalar const *const *const *aaafSource,   //!< source image (3D array)
             Scalar ***aaafDest,     //!< filtered (blurred) image stored here
             Scalar const *const *const *aaafMask,     //!< ignore voxels if aaafMask[i][j][k]==0
             Scalar sigma,                 //!< Gaussian sigma parameter σ
             Scalar truncate_ratio = 2.5,  //!< how many sigma before truncating?
             bool normalize = true,           //!< normalize the average?
             ostream *pReportProgress = NULL  //!< print progress to the user?
             )
{
  Scalar afSigma[3] = {sigma, sigma, sigma};
  ApplyGauss3D(image_size,
               aaafSource,
               aaafDest,
               aaafMask,
               afSigma,
               truncate_ratio,
               normalize,
               pReportProgress);
}






/// @brief Apply a Difference-of-Gaussians filter to a 3D image
///
/// @code
///  h(x,y,z) =   A*exp(-0.5 * ((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))
///             - B*exp(-0.5 * ((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))
/// @endcode
/// The constants "A" and "B" are determined by normalization.
/// In this version, the user manually specifies the filter window width.
/// Voxels outside the boundary of the image or outside the mask are not
/// considered during the averaging (blurring) process, and the resulting
/// after blurring is weighted accordingly.  (normalize=true by default)

template<class Scalar>
void
ApplyDog3D(int const image_size[3], //!< image size in x,y,z directions
           Scalar const *const *const *aaafSource,   //!< source image (3D array)
           Scalar ***aaafDest,     //!< filtered image stored here
           Scalar const *const *const *aaafMask,     //!< ignore voxels if aaafMask[i][j][k]==0
           Scalar const sigma_a[3], //!< (a_x, a_y, a_z) in the formula above
           Scalar const sigma_b[3], //!< (b_x, b_y, b_z) in the formula above
           int const truncate_halfwidth[3],//!<(half the filter window width along x,y,z)
           Scalar *pA = NULL, //!< Optional: report "A" (normalized coefficient) to the caller?
           Scalar *pB = NULL, //!< Optional: report "B" (normalized coefficient) to the caller?
           ostream *pReportProgress = NULL  //!< print progress to the user?
           )
{
  Scalar ***aaafTemp; //temporary array to store partially processed tomogram
  Scalar *afTemp;     //temporary array to store partially processed tomogram

  if (pReportProgress)
    *pReportProgress
      << " -- Attempting to allocate space for one more image.\n"
      << " -- (If this crashes your computer, find a computer with\n"
      << " --  more RAM and use \"ulimit\", OR use a smaller image.)\n";

  Alloc3D(image_size,
          &afTemp,
          &aaafTemp);

  Scalar A, B;        // let the user know what A B coefficients were used

  // Convolve the original source with the 1st Gaussian
  A = ApplyGauss3D(image_size,
                   aaafSource,
                   aaafDest,   // <-- save result here
                   aaafMask,
                   sigma_a,
                   truncate_halfwidth,
                   true,
                   pReportProgress);

  // Convolve the original source with the 2nd Gaussian
  B = ApplyGauss3D(image_size,
                   aaafSource,
                   aaafTemp,   // <-- save result here
                   aaafMask,
                   sigma_b,
                   truncate_halfwidth,
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

  // Report the A and B (normalization) coefficients to the caller?
  if (pA)
    *pA = A;
  if (pB)
    *pB = B;

} // ApplyDog3D()





/// @brief Apply a "scale free" Difference-of-Gaussians filter to a 3D image
///
/// @code
///   h(x,y,z) = scale * ( A*exp(-0.5*r_a^2) - B*exp(-0.5*r_b^2) )
///  where r_a = √((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))
///    and r_b = √((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))
///      a_x = σ_x*(1-0.5*δ), a_y = σ_y*(1-0.5*δ), a_z = σ_z*(1-0.5*δ)
///      b_x = σ_x*(1+0.5*δ), b_y = σ_y*(1+0.5*δ), b_z = σ_z*(1+0.5*δ)
///    scale = (1.0 / δ^2)
/// @endcode
/// The constants "A" and "B" are determined by normalization.
/// In this version, the user manually specifies the filter window width.
/// Voxels outside the boundary of the image or outside the mask are not
/// considered during the averaging (blurring) process, and the resulting
/// after blurring is weighted accordingly.  (normalize=true by default)

template<class Scalar>
void
ApplyDogScaleFree3D(int const image_size[3], //!< source image size
                    Scalar const *const *const *aaafSource,   //!< source image (3D array)
                    Scalar ***aaafDest,     //!< filtered image stored here
                    Scalar const *const *const *aaafMask,     //!< ignore voxels if aaafMask[i][j][k]==0
                    Scalar const sigma[3],  //!< Gaussian width in x,y,z drections
                    Scalar delta_sigma_over_sigma, //difference in Gauss widths
                    Scalar truncate_ratio,  //!< how many sigma before truncating?
                    Scalar *pA = NULL, //!< Optional: report "A" (normalized coefficient) to the caller?
                    Scalar *pB = NULL, //!< Optional: report "B" (normalized coefficient) to the caller?
                    ostream *pReportProgress = NULL  //!< print progress to the user?
                    )
{
  // "-log" approximates to the "Laplacian of a Gaussian" ("DOG") filter
  // with the difference of two Gaussians ("DOG") filter.
  // (The two Gaussians have widths which are slightly above and
  //  slightly below the width parameters specified by the user.)
  // For background details see:
  // https://en.wikipedia.org/wiki/Blob_detection
  // https://en.wikipedia.org/wiki/Difference_of_Gaussians
  // https://en.wikipedia.org/wiki/Mexican_hat_wavelet
  Scalar sigma_a[3];
  Scalar sigma_b[3];
  sigma_a[0] = sigma[0] * (1.0 - 0.5*delta_sigma_over_sigma);
  sigma_a[1] = sigma[1] * (1.0 - 0.5*delta_sigma_over_sigma);
  sigma_a[2] = sigma[2] * (1.0 - 0.5*delta_sigma_over_sigma);
  sigma_b[0] = sigma[0] * (1.0 + 0.5*delta_sigma_over_sigma);
  sigma_b[1] = sigma[1] * (1.0 + 0.5*delta_sigma_over_sigma);
  sigma_b[2] = sigma[2] * (1.0 + 0.5*delta_sigma_over_sigma);

  int truncate_halfwidth[3];
  for (int d=0; d < 3; d++)
    truncate_halfwidth[d] = floor(truncate_ratio *
                                  MAX(sigma_a[d],
                                      sigma_b[d]));

  ApplyDog3D(image_size,
             aaafSource,
             aaafDest,
             aaafMask,
             sigma_a,
             sigma_b,
             truncate_halfwidth,
             pA,
             pB,
             pReportProgress);

  // Do we care about the overall magnitude of the result?
  //In order to approximate the magnitude of the
  //  "Scale-Normalized Laplacian of a Gaussian"   (a.k.a. "Lnorm")
  //...it is necessary to multiply our result so far by "t / delta_t"
  //
  //https://en.wikipedia.org/wiki/Blob_detection#The_Laplacian_of_Gaussian
  //
  //...where "t" = (1/2) * log_width^2
  //         because filter_mrc includes factor of (1/2) in log_width
  //   and "delta_t" = t * delta_sigma_over_sigma^2
  //   (Note: The "(1/2)" is present in both "t" and "delta_t" so the ratio 
  //          between them is unaffected by its presence)

  Scalar t_over_delta_t = 1.0 / SQR(delta_sigma_over_sigma);

  for (int iz = 0; iz < image_size[2]; iz++)
    for (int iy = 0; iy < image_size[1]; iy++)
      for (int ix = 0; ix < image_size[0]; ix++)
        aaafDest[iz][iy][ix] *= t_over_delta_t;

  // Report the A and B (normalization) coefficients to the caller?
  // If so, we need to include the factor of t_over_delta_t that we used above.
  if (pA)
    *pA *= t_over_delta_t;
  if (pB)
    *pB *= t_over_delta_t;

  // Optional: The laplacian is the 2nd derivative w respect to voxel location
  //           Should we instead make it the 2nd derivative w respect to
  //           physical distance?  (This makes it easier to compare
  //           images taken at different resolutions.)
  //           If so, multiply A, B, and aaafDest result by 1/voxel_width^2
  //           (Don't worry about this for now.  Let the calller deal with this)
  //            

} // ApplyDogScaleFree3D()






/// @brief Apply a "scale free" Difference-of-Gaussians filter to a 3D image
///
/// @code
///   h(x,y,z) = scale * ( A*exp(-0.5*r_a^2) - B*exp(-0.5*r_b^2) )
///  where r_a = √((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))
///    and r_b = √((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))
///      a_x = σ_x*(1-0.5*δ), a_y = σ_y*(1-0.5*δ), a_z = σ_z*(1-0.5*δ)
///      b_x = σ_x*(1+0.5*δ), b_y = σ_y*(1+0.5*δ), b_z = σ_z*(1+0.5*δ)
///    scale = (1.0 / δ^2)
/// @endcode
/// The constants "A" and "B" are determined by normalization.
/// In this version, the user specifies the filter window width in units of σ.
/// Voxels outside the boundary of the image or outside the mask are not
/// considered during the averaging (blurring) process, and the resulting
/// after blurring is weighted accordingly.  (normalize=true by default)

template<class Scalar>
void
ApplyDogScaleFree3D(int const image_size[3], //!< source image size
                    Scalar const *const *const *aaafSource, //!< source image
                    Scalar ***aaafDest,     //!< save results here
                    Scalar const *const *const *aaafMask,  //!< ignore voxels where mask==0
                    Scalar sigma,  //!< Gaussian width in x,y,z drections
                    Scalar delta_sigma_over_sigma, //!< δ, difference in Gauss widths
                    Scalar truncate_ratio=2.5,  //!< how many sigma before truncating?
                    Scalar *pA = NULL, //!< Optional: report "A" (normalized coefficient) to the caller?
                    Scalar *pB = NULL, //!< Optional: report "B" (normalized coefficient) to the caller?
                    ostream *pReportProgress = NULL  //!< print progress to the user?
                    )
{
  Scalar sigma_xyz[3] = {sigma, sigma, sigma};
  ApplyDogScaleFree3D(image_size,
                      aaafSource,
                      aaafDest,
                      aaafMask,
                      sigma_xyz,
                      delta_sigma_over_sigma,
                      truncate_ratio,
                      pA,
                      pB,
                      pReportProgress);
}




/// @brief Find all scale-invariant blobs in the image as a function of sigma
///        regardless of overlap.  (Blobs with poor scores can be discarded)
///        Discard blobs with poor scores or large overlap with other blobs.
///
/// Algorithm described in:
///    Lindeberg,T., Int. J. Comput. Vision., 30(2):77-116, (1998)

template<class Scalar>
void
BlobDog(int const image_size[3], //!< source image size
        Scalar const *const *const *aaafSource,   //!< source image
        Scalar const *const *const *aaafMask,     //!< ignore voxels where mask==0
        const vector<Scalar>& blob_sigma, //!< blob widths to try, ordered
        vector<array<Scalar,3> >& minima_crds, //!< store minima x,y,z coords here
        vector<array<Scalar,3> >& maxima_crds, //!< store maxima x,y,z coords here
        vector<Scalar>& minima_sigma, //!< corresponding width for that minima
        vector<Scalar>& maxima_sigma, //!< corresponding width for that maxima
        vector<Scalar>& minima_scores, //!< what was the blob's score?
        vector<Scalar>& maxima_scores, //!< (score = intensity after filtering)
        Scalar delta_sigma_over_sigma=0.02,//!< δ param for approximating LOG with DOG
        Scalar truncate_ratio=2.8,      //!< how many sigma before truncating?
        Scalar minima_threshold=0.0,    //!< discard blobs with unremarkable scores
        Scalar maxima_threshold=0.0,    //!< discard blobs with unremarkable scores
        bool use_threshold_ratios=true, //!< threshold=ratio*best_score ?
        // optional arguments
        ostream *pReportProgress = NULL, //!< optional: report progress to the user?
        Scalar ****aaaafI = NULL, //!<optional: preallocated memory for filtered images (indexable)
        Scalar **aafI = NULL      //!<optional: preallocated memory for filtered images (contiguous)
        )

{

  // We need 3 images to store the result of filtering the image
  // using DOG filters with 3 different widths.  Store those images here:

  if (pReportProgress)
    *pReportProgress
      << " -- Attempting to allocate space for 3 more images.        --\n"
      << " -- (If this crashes your computer, find a computer with   --\n"
      << " --  more RAM and use \"ulimit\", OR use a smaller image.)   --\n";

  bool preallocated = ! (aaaafI == NULL);

  if (! preallocated) {
    assert(aaaafI == NULL);
    assert(aafI == NULL);
    aaaafI = new Scalar*** [3];
    aafI = new Scalar* [3];
    Alloc3D(image_size,
            &(aafI[0]),
            &(aaaafI[0]));
    Alloc3D(image_size,
            &(aafI[1]),
            &(aaaafI[1]));
    Alloc3D(image_size,
            &(aafI[2]),
            &(aaaafI[2]));
  }


  Scalar global_min_score = 1.0;  //impossible, minima < 0
  Scalar global_max_score = -1.0; //impossible, maxima > 0

  //REMOVE THIS CRUFT:
  //bool disable_thresholds = ((! use_threshold_ratios) &&
  //                           (minima_threshold >= maxima_threshold));

  if (pReportProgress)
    *pReportProgress << "\n----- Blob detection initiated using "
                     << blob_sigma.size() << " trial Gaussians -----\n\n";

  for (int ir = 0; ir < blob_sigma.size(); ir++) {

    if (pReportProgress)
      *pReportProgress
        << "--- Progress: "<<ir+1<<"/"<<blob_sigma.size() << "\n"
        << "--- Applying DOG filter using sigma[" << ir << "] = "
        << blob_sigma[ir] << " (in voxels) ---\n";

    // Run the image through a DOG filter with the currently selected width.
    // Compare the resulting filtered image with filtered images
    // smaller blob widths.
    // If a given voxel has a value which is larger than it's surrounding 26
    // voxels in this image, -as well as- the same 27 voxels which were filtered
    // using larger and smaller blob widths, respectively...
    //   ... -then- it is a local maximum in 4-dimensional X,Y,Z,R space.
    // Keep track of all of these X,Y,Z,R local maxima voxels and their values
    // (and the local minima as well).

    // We are going to apply a DOG filter to the original image and
    // store it in the aaaaf array.
    // Unfortunately, there probably won't be enough memory 
    // to keep track of all of the filtered versions of the image.
    // Instead we only need to keep track of the last 3 filtered images,
    // So the aaaafI[] array only has 3 entries.  We store the images in:
    //   aaaafI[ j2i[-1]]
    //   aaaafI[ j2i[0] ]
    //   aaaafI[ j2i[1] ]
    // ... and the values of j2i[-1], j2i[0], and j2i[1] cycle between 0,1,2

    int j2i_[3];       // an array of indices.  index into it using 0,1,2
    int *j2i = j2i_+1; // convenient to offset by 1 so we can index using -1,0,1
    j2i[-1] = (ir-2) % 3;
    j2i[0]  = (ir-1) % 3;
    j2i[1]  = ir % 3;

    //Apply the DOG filter using the most recent blob width:
    ApplyDogScaleFree3D(image_size,
                        aaafSource,
                        aaaafI[j2i[1]], //<-store the most recent filtered image
                        aaafMask,
                        blob_sigma[ir],
                        delta_sigma_over_sigma,
                        truncate_ratio);

    // We must have filtered at least 3 images using different blob widths
    if (ir < 2)
      continue;

    // The blob widths should be ordered consistently (increasing or decreasing)
    assert(((blob_sigma[ir-2] < blob_sigma[ir-1]) &&   //increasing order
            (blob_sigma[ir-1] < blob_sigma[ir])) ||    
           ((blob_sigma[ir-2] > blob_sigma[ir-1]) &&   //decreasing order
            (blob_sigma[ir-1] > blob_sigma[ir])));


    if (pReportProgress)
      *pReportProgress
        << "--- Searching for local minima & maxima with width[" << ir-1 << "] = "
        << blob_sigma[ir-1]*2.0*sqrt(3) << " ---\n";

    // As we are looking for local minima and maxima we need to
    // keep track of the best scores so far
    // and discard blobs with weak scores as early as possible.
    // If use_threshold_ratios=true (default), then any
    // local maxima which are not > maxima_threshold * highest_score_so_far
    // will be discarded as are minima which are not
    // < maxima_threshold * lowest_score_so_far, which should be < 0.
    // In practice this keep the list of minima and maxima of manageable size.
    // (Otherwise a significant fraction of the voxels in the image could end
    // up being a minima or maxima at some scale, and this could easily use up
    // an enormous amout of memory and be difficult to deal with later.
    // It's better to prevent this list from growing too large at any point.

    #pragma omp parallel
    {

      // (The following variables are private for each thread/processor.
      // Later, the global list of minima and maxima will be updated with the
      // information collected from each processor at the end of this iteration)

      vector<array<Scalar,3> > min_crds_proc; //store minima x,y,z coords here
      vector<array<Scalar,3> > max_crds_proc; //store maxima x,y,z coords here
      vector<Scalar> min_sigma_proc; //corresponding width for that minima
      vector<Scalar> max_sigma_proc; //corresponding width for that maxima
      vector<Scalar> min_scores_proc; //what was the blob's score?
      vector<Scalar> max_scores_proc; //(score = intensity after filtering)
      Scalar global_min_score_proc = global_min_score;
      Scalar global_max_score_proc = global_max_score;

      #pragma omp for collapse(2)
      for (int iz = 0; iz < image_size[2]; iz++) {
        for (int iy = 0; iy < image_size[1]; iy++) {
          for (int ix = 0; ix < image_size[0]; ix++) {
            // Search the 81-1 = 80 surrounding voxels to see if this voxel is
            // either a minima or a maxima in 4-dimensional x,y,z,r space
            bool is_minima = true;
            bool is_maxima = true;
            for (int jr = -1; jr <= 1; jr++) {
              for (int jz = -1; jz <= 1; jz++) {
                for (int jy = -1; jy <= 1; jy++) {
                  for (int jx = -1; jx <= 1; jx++) {
                    if (jx==0 && jy==0 && jz==0 && jr==0)
                      continue; // Skip the central voxel. Check neighbors only
                    int Ix = ix + jx;
                    int Iy = iy + jy;
                    int Iz = iz + jz;
                    // All neighbors must be available for comparison
                    if ((Ix < 0) || (Ix >= image_size[0]) ||
                        (Iy < 0) || (Iy >= image_size[1]) ||
                        (Iz < 0) || (Iz >= image_size[2]) ||
                        (aaafMask && (aaafMask[Iz][Iy][Ix] == 0))) {
                      is_minima = false;
                      is_maxima = false;
                      continue;
                    }
                    Scalar entry    = aaaafI[ j2i[0]  ][ iz ][ iy ][ ix ];
                    Scalar neighbor = aaaafI[ j2i[jr] ][ Iz ][ Iy ][ Ix ];
                    if (neighbor <= entry)
                      is_minima = false;
                    if (neighbor >= entry)
                      is_maxima = false;
                  }
                }
              }
            }

            Scalar score = aaaafI[ j2i[0] ][ iz ][ iy ][ ix ];

            if ((! aaafMask) || (aaafMask[iz][iy][ix] != 0))
            {
              Scalar minima_threshold_so_far = minima_threshold;
              if (use_threshold_ratios)
                minima_threshold_so_far=minima_threshold*global_min_score_proc;
              if (is_minima &&
                  (score < 0.0) &&
                  ((score < minima_threshold_so_far) // || disable_thresholds
                   ))   
              {
                array<Scalar, 3> ixiyiz;
                ixiyiz[0] = ix;
                ixiyiz[1] = iy;
                ixiyiz[2] = iz;
                min_crds_proc.push_back(ixiyiz);
                min_sigma_proc.push_back(blob_sigma[ir-1]);
                min_scores_proc.push_back(score);
                if (score < global_min_score_proc)
                  global_min_score_proc = score;
              }

              Scalar maxima_threshold_so_far = maxima_threshold;
              if (use_threshold_ratios)
                maxima_threshold_so_far=maxima_threshold*global_max_score_proc;
              if (is_maxima &&
                  (score > 0.0) &&
                  ((score > maxima_threshold_so_far) // || disable_thresholds
                   ))   
              {
                array<Scalar, 3> ixiyiz;
                ixiyiz[0] = ix;
                ixiyiz[1] = iy;
                ixiyiz[2] = iz;
                max_crds_proc.push_back(ixiyiz);
                max_sigma_proc.push_back(blob_sigma[ir-1]);
                max_scores_proc.push_back(score);
                if (score > global_max_score_proc)
                  global_max_score_proc = score;
              }
            }
            assert(! (is_minima && is_maxima));
          } //for (int ix=0; ix<image_size[0]; ix++) {
        } //for (int iy=0; iy<image_size[1]; iy++) {
      } //for (int iz=0; iz<image_size[2]; iz++) {

      #pragma omp critical
      {

        // Append the newly minima and maxima discovered by this processor
        // to the global list of minima and maxima:

        minima_crds.insert(minima_crds.end(),
                           min_crds_proc.begin(),
                           min_crds_proc.end());
        min_crds_proc.clear();
        minima_sigma.insert(minima_sigma.end(),
                            min_sigma_proc.begin(),
                            min_sigma_proc.end());
        min_sigma_proc.clear();
        minima_scores.insert(minima_scores.end(),
                             min_scores_proc.begin(),
                             min_scores_proc.end());
        min_scores_proc.clear();
        maxima_crds.insert(maxima_crds.end(),
                           max_crds_proc.begin(),
                           max_crds_proc.end());
        max_crds_proc.clear();
        maxima_sigma.insert(maxima_sigma.end(),
                            max_sigma_proc.begin(),
                            max_sigma_proc.end());
        max_sigma_proc.clear();
        maxima_scores.insert(maxima_scores.end(),
                             max_scores_proc.begin(),
                             max_scores_proc.end());
        max_scores_proc.clear();

        // Update the global minima and maxima as well:
        if (global_min_score > global_min_score_proc)
          global_min_score = global_min_score_proc;
        if (global_max_score < global_max_score_proc)
          global_max_score = global_max_score_proc;
      } //#pragma omp critical

    } //#pragma omp parallel

    if (pReportProgress)
      *pReportProgress
        << "--- (Found " << minima_crds.size ()
        << " and " << maxima_crds.size()
        << " local minima and maxima, respectively so far) ---\n" << endl;

    assert((minima_crds.size() == minima_sigma.size()) &&
           (minima_crds.size() == minima_scores.size()));
    assert((maxima_crds.size() == maxima_sigma.size()) &&
           (maxima_crds.size() == maxima_scores.size()));

  } //for (ir = 0; ir < blob_sigma.size(); ir++)


  if (use_threshold_ratios) {
    minima_threshold *= global_min_score;
    maxima_threshold *= global_max_score;
  }
  // Now that we know what the true global minima and maxima are,
  // go back and discard maxima whose scores are not higher than
  // maxima_threshold * global_max_score.
  // (Do the same for local minima as well.)
  int i;
  i = 0;
  while (i < minima_scores.size()) {
    assert(minima_scores[i] < 0.0);
    if (minima_scores[i] >= minima_threshold) {
      // delete this blob
      minima_crds.erase(minima_crds.begin() + i,
                        minima_crds.begin() + i + 1);
      minima_sigma.erase(minima_sigma.begin() + i,
                         minima_sigma.begin() + i + 1);
      minima_scores.erase(minima_scores.begin() + i,
                          minima_scores.begin() + i + 1);
    }
    else
      i++;
  }
  i = 0;
  while (i < maxima_scores.size()) {
    assert(maxima_scores[i] > 0.0);
    if (maxima_scores[i] <= maxima_threshold) {
      // delete this blob
      maxima_crds.erase(maxima_crds.begin() + i,
                        maxima_crds.begin() + i + 1);
      maxima_sigma.erase(maxima_sigma.begin() + i,
                         maxima_sigma.begin() + i + 1);
      maxima_scores.erase(maxima_scores.begin() + i,
                          maxima_scores.begin() + i + 1);
    }
    else
      i++;
  }


  if (! preallocated) {
    // Deallocate the temporary arrays we created earlier
    Dealloc3D(image_size,
              &aafI[0],
              &aaaafI[0]);
    Dealloc3D(image_size,
              &aafI[1],
              &aaaafI[1]);
    Dealloc3D(image_size,
              &aafI[2],
              &aaaafI[2]);
    delete [] aaaafI;
    delete [] aafI;
  }

} //BlobDog()



/// @brief Find all scale-invariant blobs in the image as a function of diameter
///        regardless of overlap.  (Blobs with poor scores can be discarded.)
///
/// Algorithm described in:
///    Lindeberg,T., Int. J. Comput. Vision., 30(2):77-116, (1998)
/// This version refers to blobs by their diameter (instead of "sigma").
/// This version can discard blobs which overlap with existing blobs.
/// (This is sometimes called "non-max suppression".)

template<class Scalar>
void
BlobDogD(int const image_size[3], //!<source image size
         Scalar const *const *const *aaafSource,   //!< source image
         Scalar const *const *const *aaafMask,     //!< ignore voxels where mask==0
         const vector<Scalar>& blob_diameters, //!<blob widths to try, ordered
         vector<array<Scalar,3> >& minima_crds, //!<store minima x,y,z coords here
         vector<array<Scalar,3> >& maxima_crds, //!<store maxima x,y,z coords here
         vector<Scalar>& minima_diameters, //!<corresponding width for that minima
         vector<Scalar>& maxima_diameters, //!<corresponding width for that maxima
         vector<Scalar>& minima_scores, //!<what was the blob's score?
         vector<Scalar>& maxima_scores, //!<(score = intensity after filtering)
         // optional arguments
         Scalar delta_sigma_over_sigma=0.02,//!<param for approximating LOG with DOG
         Scalar truncate_ratio=2.5,    //!<how many sigma before truncating?
         Scalar minima_threshold=0.5,  //!<discard blobs with unremarkable scores
         Scalar maxima_threshold=0.5,  //!<discard blobs with unremarkable scores
         bool    use_threshold_ratios=true, //!<threshold=ratio*best_score?
         ostream *pReportProgress = NULL, //!<report progress to the user?
         Scalar ****aaaafI = NULL, //!<preallocated memory for filtered images
         Scalar **aafI = NULL     //!<preallocated memory for filtered images (conserve memory)
         )
{

  vector<Scalar> minima_sigma;
  vector<Scalar> maxima_sigma;
  vector<Scalar> blob_sigma(blob_diameters.size());
  for (int i=0; i < blob_diameters.size(); i++)
    blob_sigma[i] = blob_diameters[i]/(2.0*sqrt(3));

  BlobDog(image_size,
          aaafSource,
          aaafMask,
          blob_sigma,
          minima_crds,
          maxima_crds,
          minima_sigma,
          maxima_sigma,
          minima_scores,
          maxima_scores,
          delta_sigma_over_sigma,
          truncate_ratio,
          minima_threshold,
          maxima_threshold,
          use_threshold_ratios,
          pReportProgress,
          aaaafI,
          aafI);

  minima_diameters.resize(minima_sigma.size());
  for (int i=0; i < minima_sigma.size(); i++)
    minima_diameters[i] = minima_sigma[i] * 2.0 * sqrt(3);

  maxima_diameters.resize(maxima_sigma.size());
  for (int i=0; i < maxima_sigma.size(); i++)
    maxima_diameters[i] = maxima_sigma[i] * 2.0 * sqrt(3);

} // BlobDogD()





/// @brief  Create a 3D image containing multiple spheres (or spherical shells)
///         various sizes, thicknesses, and locations, specified by the caller.
///         The resulting spheres can (optionally) be superimposed with the
///         existing image (if aaafImage contains data).  This is done by 
///         rescaling the background image voxel brightnesses.

template<class Scalar>
void
VisualizeBlobs(int const image_size[3], //!< image size
               Scalar ***aaafImage,  //!< array where we should write new image
               Scalar const *const *const *aaafMask,   //!< ignore voxels where mask==0
               vector<array<Scalar,3> > &centers, //!< center of each blob
               vector<Scalar> &diameters,         //!< diameter of eachs pherical shell (in voxels)
               vector<Scalar> &shell_thicknesses, //!< how thick is each spherical shell (in voxels)
               vector<Scalar> &voxel_intensities_foreground, //!< voxels in spherical shell get this value
               Scalar voxel_intensity_background = 0.0, //!< assign background voxels to this value
               Scalar voxel_intensity_background_rescale = 0.25, //!< superimpose with old image? Which weight?
               bool voxel_intensity_foreground_normalize = false, //!< divide brightnesses by number of voxels in spherical shell? (rarely useful)
               ostream *pReportProgress = NULL //!<optional: report progress to the user?
               )
{

  Scalar tomo_ave  =  AverageArr(image_size,
                                  aaafImage,
                                  aaafMask);
  Scalar tomo_stddev  =  StdDevArr(image_size,
                                    aaafImage,
                                    aaafMask);

  double score_sum_sq = 0.0;
  for (int i = 0; i < voxel_intensities_foreground.size(); i++)
    score_sum_sq += SQR(voxel_intensities_foreground[i]);

  double score_rms = sqrt(score_sum_sq / voxel_intensities_foreground.size());

  // The spherical shells will have brightnesses chosen according to their
  // scores. Rescale the tomogram intensities so that they are approximately
  // in the range of scores.  This way we can see both at the same time
  // (when viewing the tomogram using IMOD, for example)
  for (int iz = 0; iz < image_size[2]; iz++) {
    for (int iy = 0; iy < image_size[1]; iy++) {
      for (int ix = 0; ix < image_size[0]; ix++) {
        aaafImage[iz][iy][ix] =
          ((aaafImage[iz][iy][ix] - tomo_ave) / tomo_stddev)
          *
          score_rms
          *
          voxel_intensity_background_rescale;
      }
    }
  }


  for (int iz=0; iz<image_size[2]; iz++)
    for (int iy=0; iy<image_size[1]; iy++)
      for (int ix=0; ix<image_size[0]; ix++)
        aaafImage[iz][iy][ix] += voxel_intensity_background;


  for (int i = 0; i < centers.size(); i++) {
    if (pReportProgress)
      *pReportProgress << "processing coordinates " << i+1 << " / "
           << centers.size() << ": x,y,z(in_voxels)="
           << centers[i][0] << "," << centers[i][1] << "," << centers[i][2];

    if ((centers[i][0] < 0) || (centers[i][1] < 0) || (centers[i][2] < 0) ||
        (centers[i][0] >= image_size[0]) ||
        (centers[i][1] >= image_size[1]) ||
        (centers[i][2] >= image_size[2]))
        throw InputErr("Error: Coordinates in the text file lie outside the boundaries of the image.\n"
                       "       Did you set the voxel width correctly?\n"
                       "       (Did you use the \"-w\" argument?)\n");
      
    if (aaafMask &&
        (aaafMask[static_cast<int>(centers[i][2])]
                 [static_cast<int>(centers[i][1])]
                 [static_cast<int>(centers[i][0])] == 0.0))
      continue;

    int Rs = ceil(diameters[i]/2-0.5);
    if (Rs < 0) Rs = 0;
    Scalar Rssqr_max = SQR(diameters[i]/2);
    Scalar Rssqr_min = 0.0;
    if ((shell_thicknesses[i] > 0.0) && (diameters[i]/2 - shell_thicknesses[i] > 0.0))
      Rssqr_min = SQR(diameters[i]/2 - shell_thicknesses[i]);

    if (pReportProgress)
      *pReportProgress << ", diameter=" << diameters[i] << ", th=" << shell_thicknesses[i]
                       <<"\n";

    // Normalize the brightness of each sphere?
    // (ie by dividing the intensity by the number of voxels in the sphere)
    Scalar imultiplier = 1.0;
    long nvoxelspersphere = 1;
    if (voxel_intensity_foreground_normalize) {
      nvoxelspersphere = 0;
      for (int jz = -Rs; jz <= Rs; jz++) {
        for (int jy = -Rs; jy <= Rs; jy++) {
          for (int jx = -Rs; jx <= Rs; jx++) {
            Scalar rsqr = jx*jx + jy*jy + jz*jz;
            if ((Rssqr_min <= rsqr) && (rsqr <= Rssqr_max))
              nvoxelspersphere++;
          }
        }
      }
    }
    imultiplier = 1.0 / nvoxelspersphere;
    for (int jz = -Rs; jz <= Rs; jz++) {
      for (int jy = -Rs; jy <= Rs; jy++) {
        for (int jx = -Rs; jx <= Rs; jx++) {
          int rsqr = jx*jx + jy*jy + jz*jz;
          if (! ((Rssqr_min <= rsqr) && (rsqr <= Rssqr_max)))
            continue;
          else if ((centers[i][0] + jx < 0) ||
                   (centers[i][0] + jx >= image_size[0]) ||
                   (centers[i][1] + jy < 0) ||
                   (centers[i][1] + jy >= image_size[1]) ||
                   (centers[i][2] + jz < 0) ||
                   (centers[i][2] + jz >= image_size[2]))
            continue;
          else if (aaafMask
                   &&
                   (aaafMask[static_cast<int>(centers[i][2])+jz]
                            [static_cast<int>(centers[i][1])+jy]
                            [static_cast<int>(centers[i][0])+jx]
                    == 0.0))
            continue;
          else
            aaafImage[static_cast<int>(centers[i][2])+jz]
                     [static_cast<int>(centers[i][1])+jy]
                     [static_cast<int>(centers[i][0])+jx] =
                          voxel_intensities_foreground[i] * imultiplier;
        }
      }
    }
  } //for (int i = 0; i < centers.size(); i++) {

} //VisualizeBlobs()





/// @brief apply permutation in-place
///
/// credit: Raymond Chen
template<class T>
void
apply_permutation(vector<T>& v,
                  vector<size_t>& indices) {
  vector<T> v_copy(v);
  for (size_t i = 0; i < indices.size(); i++) {
    size_t j = indices[i];
    v[i] = v_copy[j];
  }
}




/// @brief sort blobs by their scores
template<class Scalar>
void
SortBlobs(vector<array<Scalar,3> >& blob_crds, //!< x,y,z of each blob's center
          vector<Scalar>& blob_diameters,  //!< the width of each blob
          vector<Scalar>& blob_scores,  //!< the score for each blob
          bool descending_order = true, //!<sort scores in ascending or descending order?
          ostream *pReportProgress = NULL //!< optional: report progress to the user?
          )
{ 
  size_t n_blobs = blob_crds.size();
  assert(n_blobs == blob_diameters.size());
  assert(n_blobs == blob_scores.size());
  vector<tuple<Scalar, size_t> > score_index(n_blobs);
  for (size_t i = 0; i < n_blobs; i++)
    score_index[i] = make_tuple(blob_scores[i], i);

  if (n_blobs > 0) {
    if (pReportProgress)
      *pReportProgress << "-- Sorting blobs according to their scores...";
    if (descending_order)
      sort(score_index.rbegin(),
           score_index.rend());
    else
      sort(score_index.begin(),
           score_index.end());
    vector<size_t> permutation(n_blobs);
    for (size_t i = 0; i < score_index.size(); i++)
      permutation[i] = get<1>(score_index[i]);
    score_index.clear();
    apply_permutation(blob_crds,  permutation);
    apply_permutation(blob_diameters, permutation);
    apply_permutation(blob_scores, permutation);
    if (pReportProgress)
      *pReportProgress << "done --" << endl;
  }
} //SortBlobs()





/// @brief  Calculate the volume of overlap between two spheres of radius 
///         Ri and Rj separated by a distance of rij.

template<class Scalar>
Scalar CalcVolOverlap(Scalar rij,//!<the distance between the spheres' centers
                       Scalar Ri, //!< the radius of sphere i
                       Scalar Rj  //!< the radius of sphere j
                       )
{
  // WLOG, assume Ri <= Rj.  Insure that below
  if (Ri > Rj) {
    Scalar tmp = Ri;
    Ri = Rj;
    Rj = tmp;
  }
  if (rij <= Ri) {
    return (4*M_PI/3) * Ri*Ri*Ri;
  }

  // "xi" and "xj" are the distances from the sphere centers
  // to the plane where the two spheres intersect.
  Scalar xi = 0.5 * (1.0/rij) * (rij*rij + Ri*Ri - Rj*Rj);
  Scalar xj = 0.5 * (1.0/rij) * (rij*rij + Rj*Rj - Ri*Ri);
  Scalar volume_overlap =
    (M_PI/3)*( Ri*Ri*Ri * (1 + SQR(xi/Ri)) * (1 - (xi/Ri)) +
               Rj*Rj*Rj * (1 + SQR(xj/Rj)) * (1 - (xj/Rj)) );
  return volume_overlap;
}



/// @brief these options tell "DiscardOverlappingBlobs" how to give
/// priority to different blobs based on their score (ie, minima or maxima?)
typedef enum eSortBlobCriteria {
  DO_NOT_SORT,
  PRIORITIZE_HIGH_SCORES,
  PRIORITIZE_LOW_SCORES
} SortBlobCriteria;



/// @brief nonmax suppression for blobs
///
///   Blobs can either be discarded because the distance between their centers
///   exceeds the sum of their (*min_radial_separation_ratio)
///   OR because the overlapping volume exceeds the volume of the larger sphere
///   (after multiplication by max_volume_overlap_large)
///   OR because the overlapping volume exceeds the volume of the smaller sphere
///   (after multiplication by max_volume_overlap_small)
///
///   In this function, the width of each blob is located in an array storing
///   their diameters (instead of their corresponding "sigma" values).
template<class Scalar>
void
DiscardOverlappingBlobs(vector<array<Scalar,3> >& blob_crds,
                        vector<Scalar>& blob_diameters, 
                        vector<Scalar>& blob_scores,
                        SortBlobCriteria sort_blob_criteria, //!< priority to high or low scoring blobs?
                        Scalar min_radial_separation_ratio, //!< discard blobs too close
                        Scalar max_volume_overlap_large, //!< discard blobs which overlap too much with the large blob
                        Scalar max_volume_overlap_small, //!< discard blobs which overlap too much with the small blob
                        ostream *pReportProgress = NULL, //!< report progress back to the user?
                        int scale=6 //!<occupancy_table_size shrunk by this much
                                    //!<relative to source (necessary to reduce memory usage)
                        )
{
  // Strategy:
  // 1) Sort the blobs in order of their scores: from good scores, to bad scores
  // 2) Beginning with the best scoring blobs and proceeding downwards,
  //    check to see if the current blob overlaps with any of the blobs which
  //    came before it (which score better).
  // 3) If there is an overlap, delete it from the list.

  if (sort_blob_criteria == PRIORITIZE_HIGH_SCORES)
    SortBlobs(blob_crds,
              blob_diameters, 
              blob_scores,
              true,
              pReportProgress);
  else if (sort_blob_criteria == PRIORITIZE_LOW_SCORES)
    SortBlobs(blob_crds,
              blob_diameters, 
              blob_scores,
              false,
              pReportProgress);


  if (pReportProgress)
    *pReportProgress
      << " -- Attempting to allocate space for one more image.       --\n"
      << " -- (If this crashes your computer, find a computer with   --\n"
      << " --  more RAM and use \"ulimit\", OR use a smaller image.)   --\n";

  // Occupancy table
  //     (originally named "bool ***aaabOcc")
  // This table stores a list of blobs which occupy a given voxel.

  int occupancy_table_size[3];
  int bounds_min[3] = {0, 0, 0};
  int bounds_max[3] = {-1, -1, -1};
  
  for (int i=0; i < blob_crds.size(); i++) {
    for (int d=0; d < 3; d++) {
      Scalar reff = ceil(blob_diameters[i]/2); //blob radii in units of voxels
      if ((blob_crds[i][d] - reff < bounds_min[d]) ||
          (bounds_min[d] > bounds_max[d]))
        bounds_min[d] = blob_crds[i][d] - reff;
      if ((blob_crds[i][d] + reff > bounds_max[d]) ||
          (bounds_min[d] > bounds_max[d]))
        bounds_max[d] = blob_crds[i][d] + reff;
    }
  }

  // Because it is such an enormous table, the size of the table does not equal
  // the size of the original image. Instead we reduce it by a factor of "scale"
  for (int d=0; d < 3; d++)
    occupancy_table_size[d] = (1 + bounds_max[d] - bounds_min[d]) / scale;


  //vector<vector<vector<vector size_t> > > >  
  //  vvvvOcc(occupancy_table_size[3],
  //          vector<
  //          vector<
  //          <vector size_t>(0)
  //          >(occupancy_table_size[1])
  //          >(occupancy_table_size[2])
  //          );
  //
  // The nested for-loop below is less confusing, so I use this instead:

  vector<vector<vector<vector<size_t> > > > vvvvOcc(occupancy_table_size[2]);
  for(int iz = 0; iz < occupancy_table_size[2]; iz++) {
    vvvvOcc[iz] = vector<vector<vector<size_t> > >(occupancy_table_size[1]);
    for(int iy = 0; iy < occupancy_table_size[1]; iy++) {
      vvvvOcc[iz][iy] = vector<vector<size_t> >(occupancy_table_size[0]);
      //for(int ix = 0; ix < occupancy_table_size[0]; ix++) {
      //  vvvvOcc[iz][iy][ix] = vector<size_t>(0);
      //}
    }
  }

  if (pReportProgress)
    *pReportProgress
      << "     -- done                                                --\n\n";


  if (pReportProgress)
    *pReportProgress
      << " Detecting collisions between blobs... \n";



  // Loop through all of the blobs and fill the occupancy table.  If a given
  // blob overlaps with an existing blob (presumabely with a better score),
  // then we check to see how much they overlap before deciding whether to
  // discard the new blob.

  size_t i = 0;
  while (i < blob_crds.size()) {
    size_t n_blobs = blob_crds.size();
    assert(n_blobs == blob_diameters.size());
    assert(n_blobs == blob_scores.size());
    bool discard = false;
    Scalar reff_ = blob_diameters[i]/2; //blob radii in units of voxels
    Scalar Reff_ = reff_ / scale; //blob radii expressed in "low rez" units
    int Reff = ceil(Reff_);         //round up
    int Reffsq = ceil(Reff_*Reff_);
    Scalar ix = blob_crds[i][0];      //coordinates of the center of the blob
    Scalar iy = blob_crds[i][1];
    Scalar iz = blob_crds[i][2];
    int Ix = floor((ix - bounds_min[0]) / scale); //blob center coords (lowrez version)
    int Iy = floor((iy - bounds_min[1]) / scale);
    int Iz = floor((iz - bounds_min[2]) / scale);
    for(int Jz = -Reff; Jz <= Reff && (! discard); Jz++) {
      for(int Jy = -Reff; Jy <= Reff && (! discard); Jy++) {
        for(int Jx = -Reff; Jx <= Reff && (! discard); Jx++) {
          int rsq = Jx*Jx + Jy*Jy + Jz*Jz;
          if (! ((0 <= Ix+Jx) && (Ix+Jx < occupancy_table_size[0]) &&
                 (0 <= Iy+Jy) && (Iy+Jy < occupancy_table_size[1]) &&
                 (0 <= Iz+Jz) && (Iz+Jz < occupancy_table_size[2])))
            continue;
          if (rsq > Reffsq)
            continue;

          // loop over the other blobs which occupy this location (if any)
          vector<size_t>& overlapping_blobs = vvvvOcc[Iz+Jz][Iy+Jy][Ix+Jx];
          for (size_t _k=0; _k < overlapping_blobs.size(); _k++) {
            size_t k = overlapping_blobs[_k];
            Scalar kx = blob_crds[k][0];
            Scalar ky = blob_crds[k][1];
            Scalar kz = blob_crds[k][2];
            Scalar rik = sqrt((ix-kx)*(ix-kx)+(iy-ky)*(iy-ky)+(iz-kz)*(iz-kz));
            Scalar vol_overlap = CalcVolOverlap(rik,
                                                 blob_diameters[i]/2,
                                                 blob_diameters[k]/2);
            Scalar ri = blob_diameters[i]/2;
            Scalar rk = blob_diameters[k]/2;
            if (rik < (ri + rk) * min_radial_separation_ratio)
              discard = true;
            Scalar vi = (4*M_PI/3)*(ri*ri*ri);
            Scalar vk = (4*M_PI/3)*(rk*rk*rk);
            Scalar v_large = vi;
            Scalar v_small = vk;
            if (vk > vi) {
              v_large = vk;
              v_small = vi;
            }
            if ((vol_overlap / v_large > max_volume_overlap_large) ||
                (vol_overlap / v_small > max_volume_overlap_small)) 
              discard = true;
          }
        }
      }
    }
    if (discard) {
      blob_crds.erase(blob_crds.begin() + i);
      blob_diameters.erase(blob_diameters.begin() + i);
      blob_scores.erase(blob_scores.begin() + i);
    }
    else {
      for(int Jz = -Reff; Jz <= Reff && (! discard); Jz++) {
        for(int Jy = -Reff; Jy <= Reff && (! discard); Jy++) {
          for(int Jx = -Reff; Jx <= Reff && (! discard); Jx++) {
            int rsq = Jx*Jx + Jy*Jy + Jz*Jz;
            if (! ((0 <= Ix+Jx) && (Ix+Jx < occupancy_table_size[0]) &&
                   (0 <= Iy+Jy) && (Iy+Jy < occupancy_table_size[1]) &&
                   (0 <= Iz+Jz) && (Iz+Jz < occupancy_table_size[2])))
              continue;
            if (rsq > Reffsq)
              continue;
            vvvvOcc[Iz+Jz][Iy+Jy][Ix+Jx].push_back(i);
          }
        }
      }
      i++;
    }
  } //while (i < blob_crds.size())

  if (pReportProgress)
    *pReportProgress << " done.                                 \n";
} //DiscardOverlappingBlobs()











template<class Scalar>

class CompactMultiChannelImage3D
{

private:

  Scalar **aafI;
  Scalar *afI;
  size_t n_good_voxels;
  int n_channels_per_voxel;
  int image_size[3];

public:

  Scalar ****aaaafI; // Stores the image data


  inline int
  nchannels() {
    return n_channels_per_voxel;
  }

  CompactMultiChannelImage3D(int set_n_channels_per_voxel) {
    n_channels_per_voxel = set_n_channels_per_voxel;
    n_good_voxels = 0;
    aaaafI = NULL;
  }

  CompactMultiChannelImage3D(int set_n_channels_per_voxel,
                             int const set_image_size[3],
                             Scalar const *const *const *aaafMask = NULL,
                             ostream *pReportProgress = NULL  //!< print progress to the user?
                             )
  {
    n_channels_per_voxel = set_n_channels_per_voxel;
    n_good_voxels = 0;
    aaaafI = NULL;
    Resize(image_size, aaafMask, pReportProgress);
  }

  inline void
  Resize(int const set_image_size[3],
         Scalar const *const *const *aaafMask = NULL,
         ostream *pReportProgress = NULL  //!< print progress to the user?
         )
  {
    if (aaaafI)
      Dealloc();
    Alloc(set_image_size, aaafMask, pReportProgress);
  }

  ~CompactMultiChannelImage3D() {
    Dealloc();
  }

private:

  inline void
  Alloc(int const set_image_size[3],
        Scalar const *const *const *aaafMask = NULL,
        ostream *pReportProgress = NULL  //!< print progress to the user?
        )
  {
    image_size[0] = set_image_size[0];
    image_size[1] = set_image_size[1];
    image_size[2] = set_image_size[2];

    if (pReportProgress)
      *pReportProgress
        << " -- Attempting to allocate space for a "
        << n_channels_per_voxel << "-channel image\n"
        << " -- (If this crashes your computer, find a computer with\n"
        << " --  more RAM and use \"ulimit\", OR use a smaller image.)\n";
    Alloc3D(image_size,
            &aafI,
            &aaaafI);
    for (int iz = 0; iz < image_size[2]; iz++)
      for (int iy = 0; iy < image_size[1]; iy++)
        for (int ix = 0; ix < image_size[0]; ix++)
          aaaafI[iz][iy][ix] = NULL;

    for (int iz = 0; iz < image_size[2]; iz++) {
      for (int iy = 0; iy < image_size[1]; iy++) {
        for (int ix = 0; ix < image_size[0]; ix++) {
          if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
            continue;
          n_good_voxels++;
        }
      }
    }
    afI = new Scalar[n_good_voxels * n_channels_per_voxel];
    int n = 0;
    for (int iz = 0; iz < image_size[2]; iz++) {
      for (int iy = 0; iy < image_size[1]; iy++) {
        for (int ix = 0; ix < image_size[0]; ix++) {
          if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
            continue;
          aaaafI[iz][iy][ix] =
            &(afI[n * n_channels_per_voxel]);
          n++;
        }
      }
    }

    if (pReportProgress)
      *pReportProgress
        << "        done\n" << endl;
  } //Alloc()


  inline void
  Dealloc()
  {
    delete [] afI;

    // Now delete the array of pointers (aaaafI, which pointed to afI)
    Dealloc3D(image_size,
              &aafI,
              &aaaafI);
    afI = NULL;
    aafI = NULL;
    aaaafI = NULL;
  }

}; //class CompactMultiChannelImage3D





using namespace selfadjoint_eigen3;

/// @brief  Calculate matrix of 2nd derivatives (the hessian)
///         as well as the the vector of 1st derivatives (the gradient)
///         of the source image (aaafSource), at every location where aaafMask
///         is non-zero (or everywhere if aaafMask is NULL)
///         Apply a Gaussian blur to the image (of width sigma) beforehand,
///         (truncating the blur filter at a distance of truncate_ratio*sigma
///          voxels from the center of the Gaussian).
///         Hessians and Gradients are saved in 3-dimensional arrays of
///         either "TensorContainers", or "VectorContainers".
///         Both of these objects must support array subscripting.
///         Each "TensorContainer" object is expected to behave like
///         a one-dimensional array of 6 scalars (hence the "aaaafHessian"
///         argument behaves like a 4-dimensional array)
///         Each "VectorContainer" object is expected to behave like
///         a one-dimensional array of 3 scalars.

template<class Scalar, class VectorContainer, class TensorContainer>
void
CalcHessian3D(int const image_size[3], //!< source image size
              Scalar const *const *const *aaafSource, //!< source image
              VectorContainer ***aaaafGradient,  //!< save results here (if not NULL)
              TensorContainer ***aaaafHessian, //!< save results here (if not NULL)
              Scalar const *const *const *aaafMask,  //!< ignore voxels where mask==0
              Scalar sigma,  //!< Gaussian width in x,y,z drections
              Scalar truncate_ratio=2.5,  //!< how many sigma before truncating?
              ostream *pReportProgress = NULL  //!< print progress to the user?
              )
{
  assert(aaafSource);
  assert(aaaafHessian);
  
  int truncate_halfwidth = floor(sigma * truncate_ratio);

  // Here we use the fast, sloppy way to compute gradients and Hessians:
  // First smooth the image,
  // Then infer derivatives from finite differences.
  //
  // This only works for large sigma (sigma at least 1.0).
  //
  // A slower, more rigorous approach would be to convolve the image
  // with the derivatives of a Gaussian.  (If we need to analyze images
  // with very thin, closely spaced objects, we can try that way later.)

  // First, apply the Gaussian filter to the image
  if (pReportProgress)
    *pReportProgress
      << " -- Attempting to allocate space for one more image.       --\n"
      << " -- (If this crashes your computer, find a computer with   --\n"
      << " --  more RAM and use \"ulimit\", OR use a smaller image.)   --\n";

  Scalar ***aaafSmoothed;
  Scalar *afSmoothed;
  Alloc3D(image_size,
          &afSmoothed,
          &aaafSmoothed);

  if (pReportProgress)
    *pReportProgress << "\n";

  ApplyGauss3D(image_size,
               aaafSource,
               aaafSmoothed,
               aaafMask,
               sigma,
               truncate_halfwidth,
               false,
               pReportProgress);

  assert(image_size[0] >= 3);
  assert(image_size[1] >= 3);
  assert(image_size[2] >= 3);

  // Now compute gradients and hessians

  for (int iz = 1; iz < image_size[2]-1; iz++) {
    #pragma omp parallel for collapse(2)
    for (int iy = 1; iy < image_size[1]-1; iy++) {
      for (int ix = 1; ix < image_size[0]-1; ix++) {
        if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
          continue;

        if (aaaafGradient) {
          Scalar gradient[3];
          gradient[0]=0.5*(aaafSmoothed[iz][iy][ix+1] - 
                           aaafSmoothed[iz][iy][ix-1]);
          gradient[1]=0.5*(aaafSmoothed[iz][iy+1][ix] - 
                           aaafSmoothed[iz][iy-1][ix]);
          gradient[2]=0.5*(aaafSmoothed[iz+1][iy][ix] - 
                           aaafSmoothed[iz-1][iy][ix]);

          // Optional: Insure that the resulting gradient is dimensionless:
          // (Lindeberg 1993 "On Scale Selection for Differential Operators")

          gradient[0] *= sigma;
          gradient[1] *= sigma;
          gradient[2] *= sigma;
          
          aaaafGradient[iz][iy][ix][0] = gradient[0];
          aaaafGradient[iz][iy][ix][1] = gradient[1];
          aaaafGradient[iz][iy][ix][2] = gradient[2];

          #ifndef NDEBUG
          if ((ix==image_size[0]/2) &&
              (iy==image_size[1]/2) &&
              (iz==image_size[2]/2))
          {
            cerr << "[iz][iy][ix]=["<<iz<<"]["<<iy<<"]["<<ix<<"], "
                 << "gradient[0] = " << aaaafGradient[iz][iy][ix][0] << endl;
          }
          #endif  //#ifndef NDEBUG

        }

        if (aaaafHessian) {
          Scalar hessian[3][3];
          hessian[0][0] = (aaafSmoothed[iz][iy][ix+1] + 
                           aaafSmoothed[iz][iy][ix-1] - 
                           2*aaafSmoothed[iz][iy][ix]);
          hessian[1][1] = (aaafSmoothed[iz][iy+1][ix] + 
                           aaafSmoothed[iz][iy-1][ix] - 
                           2*aaafSmoothed[iz][iy][ix]);
          hessian[2][2] = (aaafSmoothed[iz+1][iy][ix] + 
                           aaafSmoothed[iz-1][iy][ix] - 
                           2*aaafSmoothed[iz][iy][ix]);

          hessian[0][1] = 0.25 * (aaafSmoothed[iz][iy+1][ix+1] + 
                                  aaafSmoothed[iz][iy-1][ix-1] - 
                                  aaafSmoothed[iz][iy-1][ix+1] - 
                                  aaafSmoothed[iz][iy+1][ix-1]);
          hessian[1][0] = hessian[0][1];

          hessian[1][2] = 0.25 * (aaafSmoothed[iz+1][iy+1][ix] + 
                                  aaafSmoothed[iz-1][iy-1][ix] - 
                                  aaafSmoothed[iz-1][iy+1][ix] - 
                                  aaafSmoothed[iz+1][iy-1][ix]);
          hessian[2][1] = hessian[1][2];

          hessian[2][0] = 0.25 * (aaafSmoothed[iz+1][iy][ix+1] + 
                                  aaafSmoothed[iz-1][iy][ix-1] - 
                                  aaafSmoothed[iz+1][iy][ix-1] - 
                                  aaafSmoothed[iz-1][iy][ix+1]);
          hessian[0][2] = hessian[2][0];



          #ifndef NDEBUG
          // DEBUG: REMOVE THE NEXT IF STATMENT AFTER DEBUGGING IS FINISHED
          if ((ix==image_size[0]/2) &&
              (iy==image_size[1]/2) &&
              (iz==image_size[2]/2))
          {
            cerr << "[iz][iy][ix]=["<<iz<<"]["<<iy<<"]["<<ix<<"], "
                 << "hessian[0][0] = " << hessian[0][0] << endl;
          }
          #endif  //#ifndef NDEBUG


          // Optional: Insure that the result is dimensionless:
          // (Lindeberg 1993 "On Scale Selection for Differential Operators")
          for (int di=0; di < 3; di++)
            for (int dj=0; dj < 3; dj++)
              hessian[di][dj] *= sigma*sigma;

          // To reduce memory consumption,
          // save the resulting 3x3 matrix in a smaller 1-D array whose index
          // is given by MapIndices_3x3_to_linear[][]
          for (int di = 0; di < 3; di++)
            for (int dj = di; dj < 3; dj++)
              aaaafHessian[iz][iy][ix][ MapIndices_3x3_to_linear[di][dj] ]
                = hessian[di][dj];

        } //if (aaaafHessian)
      } //for (int ix = 1; ix < image_size[0]-1; ix++)
    } //for (int iy = 1; iy < image_size[1]-1; iy++)
  } //for (int iz = 1; iz < image_size[2]-1; iz++)

  Dealloc3D(image_size,
            &afSmoothed,
            &aaafSmoothed);

  if (pReportProgress)
    *pReportProgress << "done ----" << endl;

} //CalcHessian3D()






/// CalcMomentTensor3D()
/// This is almost certainly algebraically equivalent to CalcHessian3D()
/// However this version of the function might be more robust for small ridges.
/// (This is because I apply the derivative to the Gaussian filter
///  before applying the filter, ...instead of applying the Gaussian filter
///  and then taking differences afterwards.  If the width of the object being
///  detected is not much more than 3-voxel wide, the 3-voxel wide differences
///  used in the other implementation are a large source of error.)
/// Unfortunately this version is slower and needs much more memory however.
/// Eventually, I might elliminate one of these implementations.

template<class Scalar, class FirstMomentContainer, class SecondMomentContainer>
void
CalcMomentTensor3D(int const image_size[3], //!< source image size
                   Scalar const *const *const *aaafSource, //!< source image
                   FirstMomentContainer ***aaaaf1stMoment,  //!< save results here (if not NULL)
                   SecondMomentContainer ***aaaaf2ndMoment, //!< save results here (if not NULL)
                   Scalar const *const *const *aaafMask,  //!< ignore voxels where mask==0
                   Scalar sigma,  //!< Gaussian width in x,y,z drections
                   Scalar truncate_ratio=2.5,  //!< how many sigma before truncating?
                   ostream *pReportProgress = NULL  //!< print progress to the user?
                   )
{
  assert(image_size);
  assert(aaafSource);
  
  int truncate_halfwidth = floor(sigma * truncate_ratio);

  //calculate the filter used for the 0'th moment (ordinary Gaussian filter)
  Filter1D<Scalar, int> filter0 = GenFilterGauss1D(sigma,
                                                    truncate_halfwidth);

  //calculate the filter used for the 1st moment (Guassian(x) * x)
  Filter1D<Scalar, int> filter1 = GenFilterGauss1D(sigma,
                                                    truncate_halfwidth);
  for (int i = -truncate_halfwidth; i <= truncate_halfwidth; i++)
    filter1.afW[i] *= i;


  //calculate the filter used for the 2nd moment (Guassian(x) * x^2)
  Filter1D<Scalar, int> filter2 = GenFilterGauss1D(sigma,
                                                    truncate_halfwidth);
  for (int i = -truncate_halfwidth; i <= truncate_halfwidth; i++)
    filter2.afW[i] *= i*i;


  Scalar ***aaafNorm;
  Scalar *afNorm;
  Alloc3D(image_size,
          &afNorm,
          &aaafNorm);

  if (aaafMask) {
    // Calculate the normalization we need by blurring the mask by the Gaussian
    ApplyGauss3D(image_size,
                 aaafMask,   // <-- use the mask as the source image
                 aaafNorm,   // <-- save result here
                 aaafMask,
                 sigma, // width of Gaussian
                 truncate_halfwidth,
                 true,
                 pReportProgress);
  }


  Filter1D<Scalar, int> aFilter[3];


  if (aaaaf1stMoment) {

    if (pReportProgress)
      *pReportProgress
        << " -- Attempting to allocate space for 3 more images.        --\n"
        << " -- (If this crashes your computer, find a computer with   --\n"
        << " --  more RAM and use \"ulimit\", OR use a smaller image.) --\n";

    Scalar ***aaafIx;
    Scalar *afIx;
    Alloc3D(image_size,
            &afIx,
            &aaafIx);

    Scalar ***aaafIy;
    Scalar *afIy;
    Alloc3D(image_size,
            &afIy,
            &aaafIy);

    Scalar ***aaafIz;
    Scalar *afIz;
    Alloc3D(image_size,
            &afIz,
            &aaafIz);

    // calculate the x moment (=x^1 * y^0 * z^0)
    aFilter[0] = filter1;  // x^1
    aFilter[1] = filter0;  // y^0
    aFilter[2] = filter0;  // z^0
    return ApplySeparable3D(image_size, 
                            aaafSource,
                            aaafIx,
                            aaafMask,
                            aFilter,
                            (aaafMask == NULL), //don't normalize if theres a mask
                            pReportProgress);

    // calculate the y moment (=x^0 * y^1 * z^0)
    aFilter[0] = filter0;  // x^0
    aFilter[1] = filter1;  // y^1
    aFilter[2] = filter0;  // z^0
    return ApplySeparable3D(image_size, 
                            aaafSource,
                            aaafIy,
                            aaafMask,
                            aFilter,
                            (aaafMask == NULL), //don't normalize if theres a mask
                            pReportProgress);

    // calculate the x moment (=x^0 * y^0 * z^1)
    aFilter[0] = filter0;  // x^0
    aFilter[1] = filter0;  // y^0
    aFilter[2] = filter1;  // z^1
    return ApplySeparable3D(image_size, 
                            aaafSource,
                            aaafIz,
                            aaafMask,
                            aFilter,
                            (aaafMask == NULL), //don't normalize if theres a mask
                            pReportProgress);

    if (aaafMask) {
      // If a mask was specified, then divide by the contribution from the mask
      // (An infinite sized mask should result in a contribution/weight of 1)
      for(int iz=0; iz<image_size[2]; iz++) {
        for(int iy=0; iy<image_size[1]; iy++) {
          for(int ix=0; ix<image_size[0]; ix++) {
            if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
              continue;
            assert(aaafNorm[iz][iy][ix] > 0.0);
            aaafIx[iz][iy][ix] /= aaafNorm[iz][iy][ix];
            aaafIy[iz][iy][ix] /= aaafNorm[iz][iy][ix];
            aaafIz[iz][iy][ix] /= aaafNorm[iz][iy][ix];
          }
        }
      }
    }

    for (int iz = 1; iz < image_size[2]-1; iz++) {
      for (int iy = 1; iy < image_size[1]-1; iy++) {
        for (int ix = 1; ix < image_size[0]-1; ix++) {
          if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
            continue;
          Scalar first_deriv[3];
          first_deriv[0] = aaafIx[iz][iy][ix];
          first_deriv[1] = aaafIy[iz][iy][ix];
          first_deriv[2] = aaafIz[iz][iy][ix];

          // Optional: Insure that the resulting first_deriv is dimensionless:
          // (Lindeberg 1993 "On Scale Selection for Differential Operators")
          first_deriv[0] *= sigma;
          first_deriv[1] *= sigma;
          first_deriv[2] *= sigma;
          
          // Store the result in aaaaf1stMoment[]
          aaaaf1stMoment[iz][iy][ix][0] = first_deriv[0];
          aaaaf1stMoment[iz][iy][ix][1] = first_deriv[1];
          aaaaf1stMoment[iz][iy][ix][2] = first_deriv[2];
        }
      }
    }
    Dealloc3D(image_size,
              &afIx,
              &aaafIx);
    Dealloc3D(image_size,
              &afIy,
              &aaafIy);
    Dealloc3D(image_size,
              &afIz,
              &aaafIz);
  } //if (pFirstMoment)



  if (aaaaf2ndMoment) {
    if (pReportProgress)
      *pReportProgress << "\n"
        " ------ Calculating the average of nearby voxels: ------\n";
    // P = original image (after subtracting average nearby intensities):

    Scalar ***aaafP;
    Scalar *afP;
    Alloc3D(image_size,
            &afP,
            &aaafP);

    ApplyGauss3D(image_size,
                 aaafSource,
                 aaafP,   // <-- save result here
                 aaafMask,
                 sigma, // width of Gaussian
                 truncate_halfwidth,
                 true,
                 pReportProgress);

    // Subtract the average value from the image intensity, and store in P:
    for(int iz=0; iz<image_size[2]; iz++)
      for(int iy=0; iy<image_size[1]; iy++)
        for(int ix=0; ix<image_size[0]; ix++)
          aaafP[iz][iy][ix] = aaafSource[iz][iy][ix] - aaafP[iz][iy][ix];

    if (pReportProgress)
      *pReportProgress
        << " -- Attempting to allocate space for 6 more images.        --\n"
        << " -- (If this crashes your computer, find a computer with   --\n"
        << " --  more RAM and use \"ulimit\", OR use a smaller image.) --\n";

    Scalar ***aaafIxx;
    Scalar *afIxx;
    Alloc3D(image_size,
            &afIxx,
            &aaafIxx);

    Scalar ***aaafIyy;
    Scalar *afIyy;
    Alloc3D(image_size,
            &afIyy,
            &aaafIyy);

    Scalar ***aaafIzz;
    Scalar *afIzz;
    Alloc3D(image_size,
            &afIzz,
            &aaafIzz);

    Scalar ***aaafIxy;
    Scalar *afIxy;
    Alloc3D(image_size,
            &afIxy,
            &aaafIxy);

    Scalar ***aaafIyz;
    Scalar *afIyz;
    Alloc3D(image_size,
            &afIyz,
            &aaafIyz);

    Scalar ***aaafIxz;
    Scalar *afIxz;
    Alloc3D(image_size,
            &afIxz,
            &aaafIxz);

    // calculate the x*x moment of the innertia (=x^2 * y^0 * z^0)
    aFilter[0] = filter2;  // x^2
    aFilter[1] = filter0;  // y^0
    aFilter[2] = filter0;  // z^0
    return ApplySeparable3D(image_size, 
                            aaafP,
                            aaafIxx,
                            aaafMask,
                            aFilter,
                            (aaafMask == NULL), //don't normalize if theres a mask
                            pReportProgress);

    // calculate the y*y moment of the innertia (=x^0 * y^2 * z^0)
    aFilter[0] = filter0;  // x^0
    aFilter[1] = filter2;  // y^2
    aFilter[2] = filter0;  // z^0
    return ApplySeparable3D(image_size, 
                            aaafP,
                            aaafIyy,
                            aaafMask,
                            aFilter,
                            (aaafMask == NULL), //don't normalize if theres a mask
                            pReportProgress);

    // calculate the z*z moment of the innertia (=x^0 * y^0 * z^2)
    aFilter[0] = filter0;  // x^0
    aFilter[1] = filter0;  // y^0
    aFilter[2] = filter2;  // z^2
    return ApplySeparable3D(image_size, 
                            aaafP,
                            aaafIzz,
                            aaafMask,
                            aFilter,
                            (aaafMask==NULL), //don't normalize if theres a mask
                            pReportProgress);

    // calculate the x*y moment of the innertia (=x^1 * y^1 * z^0)
    aFilter[0] = filter1;  // x^1
    aFilter[1] = filter1;  // y^1
    aFilter[2] = filter0;  // z^0
    return ApplySeparable3D(image_size, 
                            aaafP,
                            aaafIxy,
                            aaafMask,
                            aFilter,
                            (aaafMask==NULL), //don't normalize if theres a mask
                            pReportProgress);

    // calculate the y*z moment of the innertia (=x^0 * y^1 * z^1)
    aFilter[0] = filter0;  // x^0
    aFilter[1] = filter1;  // y^1
    aFilter[2] = filter1;  // z^1
    return ApplySeparable3D(image_size, 
                            aaafP,
                            aaafIyz,
                            aaafMask,
                            aFilter,
                            (aaafMask==NULL), //don't normalize if theres a mask
                            pReportProgress);

    // calculate the x*z moment of the innertia (=x^1 * y^0 * z^1)
    aFilter[0] = filter1;  // x^1
    aFilter[1] = filter0;  // y^0
    aFilter[2] = filter1;  // z^1
    return ApplySeparable3D(image_size, 
                            aaafP,
                            aaafIxz,
                            aaafMask,
                            aFilter,
                            (aaafMask==NULL), //don't normalize if theres a mask
                            pReportProgress);

    if (aaafMask) {
      // If a mask was specified, then divide by the contribution from the mask
      // (An infinite sized mask should result in a contribution/weight of 1)
      for(int iz=0; iz<image_size[2]; iz++) {
        for(int iy=0; iy<image_size[1]; iy++) {
          for(int ix=0; ix<image_size[0]; ix++) {
            if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
              continue;
            assert(aaafNorm[iz][iy][ix] > 0.0);
            aaafIxx[iz][iy][ix] /= aaafNorm[iz][iy][ix];
            aaafIyy[iz][iy][ix] /= aaafNorm[iz][iy][ix];
            aaafIzz[iz][iy][ix] /= aaafNorm[iz][iy][ix];
            aaafIxy[iz][iy][ix] /= aaafNorm[iz][iy][ix];
            aaafIyz[iz][iy][ix] /= aaafNorm[iz][iy][ix];
            aaafIxz[iz][iy][ix] /= aaafNorm[iz][iy][ix];
          }
        }
      }
    }

    for(int iz=0; iz<image_size[2]; iz++) {
      for(int iy=0; iy<image_size[1]; iy++) {
        for(int ix=0; ix<image_size[0]; ix++) {
          if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
            continue;
          // Store the result in aaaaf2ndMoment[].  As usual, to save memory
          // we use the "MapIndices_3x3_to_linear[][]" function to store the
          // entries of the symmetric 3x3 matrix in a 1D array with only 6 entries.
          // (Symmetric 3x3 matrices can have at most 6 unique entries.)
          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[0][0] ] =
            aaafIxx[ix][iy][iz];

          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[0][1] ] =
            aaafIxy[ix][iy][iz];
    
          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[0][2] ] =
            aaafIxz[ix][iy][iz];

          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[1][0] ] =
            aaafIxy[ix][iy][iz];  // (redundant but harmless)

          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[1][1] ] =
            aaafIyy[ix][iy][iz];

          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[1][2] ] =
            aaafIyz[ix][iy][iz];

          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[2][0] ] =
            aaafIxz[ix][iy][iz];  // (redundant but harmless)

          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[2][1] ] =
            aaafIyz[ix][iy][iz];  // (redundant but harmless)

          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[2][2] ] =
            aaafIzz[ix][iy][iz];
        }
      }
    }

    if (pReportProgress)
      *pReportProgress << "done ----" << endl;

    Dealoc3D(image_size,
             &afP,
             &aaafP);

    Dealloc3D(image_size,
              &afIxx,
              &aaafIxx);
    Dealloc3D(image_size,
              &afIyy,
              &aaafIyy);
    Dealloc3D(image_size,
              &afIzz,
              &aaafIzz);
    Dealloc3D(image_size,
              &afIxy,
              &aaafIxy);
    Dealloc3D(image_size,
              &afIxz,
              &aaafIxz);
    Dealloc3D(image_size,
              &afIyz,
              &aaafIyz);
  } //if (pSecondMoment)

  Dealloc3D(image_size,
            &afNorm,
            &aaafNorm);

} // CalcMomentTensor3D()







/// @brief  Convert a volumetric 3D 6-channel image, where each voxel in the
///         image has the (non-redundant) components of a symmetrix 3x3 matrix.
///         The output of this function is another 3D 6-channel image, however
///         each voxel in this image contains the 3-eigenvalues as well as the
///         eigevectors (stored as 3 Shoemake coordinates).
///         If a non-NULL "aaafMask" argument was specified, voxels in the
///         image are ignored when aaafMask[iz][iy][ix] == 0.
/// @note   The "TensorContainer" object type is expected to behave like
///         a one-dimensional array of 6 scalars.

template<class Scalar, class TensorContainer>
void
DiagonalizeHessianImage3D(int const image_size[3], //!< source image size
                          TensorContainer const *const *const *aaaafSource, //!< input tensor
                          TensorContainer ***aaaafDest, //!< output tensors stored here (can be the same as aaaafSource)
                          Scalar const *const *const *aaafMask,  //!< ignore voxels where mask==0
                          EigenOrderType eival_order = selfadjoint_eigen3::INCREASING_EIVALS,
                          ostream *pReportProgress = NULL  //!< print progress to the user?
                          )
{
  assert(aaaafSource);
  if (pReportProgress && aaaafSource)
    *pReportProgress << "\n"
      "---- Diagonalizing the Hessian everywhere (within the mask)... "
                     << flush;

  #pragma omp parallel for collapse(2)
  for (int iz = 1; iz < image_size[2]-1; iz++) {
    for (int iy = 1; iy < image_size[1]-1; iy++) {
      for (int ix = 1; ix < image_size[0]-1; ix++) {
        if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
          continue;


        #ifndef NDEBUG
        // REMOVE THE NEXT IF STATEMENT AFTER YOU ARE THROUGH DEBUGGING:
        if ((ix==image_size[0]/2) &&
            (iy==image_size[1]/2) &&
            (iz==image_size[2]/2))
        {
          Scalar hessian[3][3];
          for (int di=0; di<3; di++)
            for (int dj=0; dj<3; dj++)
              hessian[di][dj] = aaaafSource[iz][iy][ix][ MapIndices_3x3_to_linear[di][dj] ];
          Scalar quat[4];
          cerr << "[iz][iy][ix]=["<<iz<<"]["<<iy<<"]["<<ix<<"]\n"
               << "hessian = \n"
               << "    "<<hessian[0][0]<<","<<hessian[0][1]<<","<<hessian[0][2]<<"\n"
               << "    "<<hessian[1][0]<<","<<hessian[1][1]<<","<<hessian[1][2]<<"\n"
               << "    "<<hessian[2][0]<<","<<hessian[2][1]<<","<<hessian[2][2]<<"\n";
          float eivals[3];
          float eivects[3][3];
          DiagonalizeSym3(hessian,
                          eivals,
                          eivects,
                          eival_order);
          if (Determinant3(eivects) < 0.0) {
            for (int d=0; d<3; d++)
              eivects[0][d] *= -1.0;
          }
          cerr << "eivects = \n"
               << "    "<<eivects[0][0]<<","<<eivects[0][1]<<","<<eivects[0][2]<<"\n"
               << "    "<<eivects[1][0]<<","<<eivects[1][1]<<","<<eivects[1][2]<<"\n"
               << "    "<<eivects[2][0]<<","<<eivects[2][1]<<","<<eivects[2][2]<<"\n";

          // Note: Each eigenvector is a currently row-vector in eivects[3][3];
          // It's optional, but I prefer to transpose this, because I think of
          // each eigenvector as a column vector.  Either way should work.
          Transpose3(eivects);
          Matrix2Quaternion(eivects, quat); //convert to 3x3 matrix
          cerr << "quat = " << quat[0]<<","<<quat[1]<<","<<quat[2]<<","<<quat[2]<<"\n";
        }
        #endif  //#ifndef NDEBUG




        DiagonalizeSymCompact3(aaaafSource[iz][iy][ix],
                               aaaafDest[iz][iy][ix],
                               eival_order);


        #ifndef NDEBUG
        // REMOVE THE NEXT IF STATEMENT AFTER YOU ARE THROUGH DEBUGGING:
        if ((ix==image_size[0]/2) &&
            (iy==image_size[1]/2) &&
            (iz==image_size[2]/2))
        {
          float shoemake[3]; // <- the eigenvectors stored in "Shoemake" format
          shoemake[0]       = aaaafDest[iz][iy][ix][3];
          shoemake[1]       = aaaafDest[iz][iy][ix][4];
          shoemake[2]       = aaaafDest[iz][iy][ix][5];
          float quat[4];
          Shoemake2Quaternion(shoemake, quat); //convert to 3x3 matrix
          float eivects[3][3];
          Quaternion2Matrix(quat, eivects); //convert to 3x3 matrix
          cerr << "[iz][iy][ix]=["<<iz<<"]["<<iy<<"]["<<ix<<"]\n"
               << "quat2 = " << quat[0]<<","<<quat[1]<<","<<quat[2]<<","<<quat[2]<<"\n"
               << "eivects = \n"
               << "    "<<eivects[0][0]<<","<<eivects[0][1]<<","<<eivects[0][2]<<"\n"
               << "    "<<eivects[1][0]<<","<<eivects[1][1]<<","<<eivects[1][2]<<"\n"
               << "    "<<eivects[2][0]<<","<<eivects[2][1]<<","<<eivects[2][2]<<"\n"
               << endl;
        }
        #endif  //#ifndef NDEBUG

      } //for (int ix = 1; ix < image_size[0]-1; ix++) {
    } //for (int iy = 1; iy < image_size[1]-1; iy++) {
  } //for (int iz = 1; iz < image_size[2]-1; iz++) {
} //DiagonalizeHessianImage3D()






/// @brief  Read a volumetric 3D 6-channel image, where each voxel in the
///         image has the (non-redundant) components of a symmetrix 3x3 matrix
///         which has been diagonalized using DiagonalizeHessianImage3D().
///         This function loops over all the voxels and performs the
///         inverse operation.
/// @note   The "TensorContainer" object type is expected to behave like
///         a one-dimensional array of 6 scalars.

template<class Scalar, class TensorContainer>
void
UndiagonalizeHessianImage3D(int const image_size[3],  //!< source image size
                            TensorContainer const *const *const *aaaafSource, //!< input tensor
                            TensorContainer ***aaaafDest, //!< output tensors stored here (can be the same as aaaafSource)
                            Scalar const *const *const *aaafMask,  //!< ignore voxels where mask==0
                            ostream *pReportProgress = NULL  //!< print progress to the user?
                            )
{
  assert(aaaafSource);

  if (pReportProgress)
    *pReportProgress << "\n"
      "---- Undiagonalizing the Hessian everywhere (within the mask)... "
                     << flush;

  #pragma omp parallel for collapse(2)
  for (int iz = 1; iz < image_size[2]-1; iz++) {
    for (int iy = 1; iy < image_size[1]-1; iy++) {
      for (int ix = 1; ix < image_size[0]-1; ix++) {
        if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
          continue;

        UndiagonalizeSymCompact3(aaaafSource, aaaafDest);

      } //for (int ix = 1; ix < image_size[0]-1; ix++) {
    } //for (int iy = 1; iy < image_size[1]-1; iy++) {
  } //for (int iz = 1; iz < image_size[2]-1; iz++) {
} //UndiagonalizeHessianImage3D()





/// @brief: Calculate how "plane"-like a feature along a ridge is
///         from the hessian (the matrix of 2nd derivatives).
///         This function assumes the hessian matrix has been diagonalized
///         and that the first 3 entries of the "diagonalizedHessian"
///         argument are the eigenvalues of the original hessian matrix.

template<class TensorContainer, class VectorContainer>
double
ScoreHessianPlanar(TensorContainer diagonalizedHessian,
                   VectorContainer gradient)
{
  //typedef decltype(diagonalizedHessian[0]) Scalar;
  double lambda1 = diagonalizedHessian[0];
  double lambda2 = diagonalizedHessian[1];
  double lambda3 = diagonalizedHessian[2];

  // not needed:
  //Scalar grad_sqd = 0.0;
  //if (gradient)
  //  grad_sqd = DotProduct3(gradient, gradient);

  // REMOVE THIS CRUFT
  // The "score_ratio" variable is the score function used in Eq(5) of
  // Martinez-Sanchez++Fernandez_JStructBiol2011.  IT PERFORMS POORLY.
  //float score_ratio;
  //score_ratio = ((abs(lambda1) - sqrt(abs(lambda2*lambda3)))
  //               / grad_sqd);
  //score_ratio *= score_ratio;

  // REMOVE THIS CRUFT:
  // The following "linear" metric produces interesting results, but
  // the resulting membrane structures that are detected are not well
  // separated from the huge amount of background noise.  DON'T USE:
  //float Linear_norm = lambda1 - lambda2;
  //score = Linear_norm / grad_sqd;


  // I decided to try the "Ngamma_norm" metric proposed on p.26 of
  // Lindeberg Int.J.ComputVis.1998,
  // "Edge and ridge detection with automatic scale selection"
  double Nnorm = lambda1*lambda1 - lambda2*lambda2;
  Nnorm *= Nnorm;

  return Nnorm;
} // ScoreHessianPlanar()





/// @brief: Calculate how "plane"-like a feature along a ridge is
///         from the tensor created by the process of tensor voting.
///         This function assumes that "diagonalizedMatrix3x3" has been 
///         diagonalized, and that its first 3 entries 
///         are the eigenvalues of the original hessian matrix.
template<class TensorContainer>
double
ScoreTensorPlanar(const TensorContainer diagonalizedMatrix3)
{
  //return ScoreHessianPlanar(diagonalizedMatrix3x3,
  //                          NULL,
  //                          multiplier);
  double lambda1 = diagonalizedMatrix3[0];
  double lambda2 = diagonalizedMatrix3[1];
  return lambda1 - lambda2;  // the "stickness" (See TensorVoting paper)
}




/// @brief  A class for performing simple tensor-voting image processing 
///         tasksin 3D.  Currently only "stick" voting is supported.
///         (Other kinds of tensor voting, such as "plate" and "ball" are not.)
///         This can perform tensor voting for both types of
///         "stick" fields in 3D:
///         (1) Stick-fields corresponding to planar-surface-like features,
///         (2) Stick-fields corresponding to curve-like features.

template<class Scalar, class Integer, class VectorContainer, class TensorContainer>

class TV3D {

private:

  Scalar sigma;
  Integer exponent;
  Integer halfwidth[3];
  Integer array_size[3];
  Filter3D<Scalar, Integer> radial_decay_lookup;
  array<Scalar, 3> ***aaaafDisplacement;
  array<Scalar, 3> *aafDisplacement;

public:

  TV3D():radial_decay_lookup() {
    Init();
  }

  TV3D(Scalar set_sigma,
       Integer set_exponent,
       Scalar filter_cutoff_ratio=2.5):radial_decay_lookup() {
    Init();
    SetExponent(set_exponent);
    SetSigma(set_sigma, filter_cutoff_ratio);
  }

  inline TV3D(const Filter3D<Scalar, Integer>& source) {
    Resize(source.halfwidth); // allocates and initializes afH and aaafH
    std::copy(source.aafDisplacement,
              source.aafDisplacement + (array_size[0] * array_size[1] * array_size[2]),
              aafDisplacement);
  }

  ~TV3D() {
    DeallocDisplacement();
  }

  void SetExponent(Scalar set_exponent) {
    exponent = set_exponent;
  }

  void SetSigma(Scalar set_sigma, Scalar filter_cutoff_ratio=2.5)
  {
    sigma = set_sigma;
    for (int d=0; d<3; d++)
      halfwidth[d] = floor(sigma * filter_cutoff_ratio);
    Resize(halfwidth);
  }


  /// @brief  Perform dense stick-voting, using every voxel in the image
  ///         (with non-zero correspondin entries in the aaafMaskSource array)
  ///         as a source, and collecting votes at every voxel
  ///         (with non-zero correspondin entries in the aaafMaskDest array)
  ///         in the aaaafDest array.
  ///         Other kinds of tensor voting ("plate" and "ball")
  ///         are not supported by this function.
  ///         This function can perform tensor voting for both types of
  ///         "stick" fields in 3D:
  ///         (1) Stick-fields corresponding to planar-surface-like features,
  ///         (2) Stick-fields corresponding to curve-like features.
  ///         This function expects a 3D array of normalized "vectors"
  ///         (one vector for each voxel in the original image, 3 numbers each).
  ///         These "vectors" can have user-defined type (implementation),
  ///         however they must support 1-dimensional subscripting (i=0,1,2).
  ///         This function also expects an array of numbers ("saliencies")
  ///         which store the strength of each vector.
  ///         (Note: Do not multiply each vector by its saliency.
  ///                The length of each vector is expected to be 1.0)
  /// After this function is invoked, aaaafDest will store an array of
  /// diagonalized tensors, one tensor for each voxel in the original image
  /// (unless aaafMaskDest!=NULL and the corresponding entry there is 0).

  void
  TVDenseStick(Integer const image_size[3],  //!< source image size
               Scalar const *const *const *aaafSaliency,  //!< saliency (score) of each voxel (usually based on Hessian eigenvalues)
               VectorContainer const *const *const *aaaafV,  //!< vector associated with each voxel
               TensorContainer ***aaaafDest,  //!< votes will be collected here
               Scalar const *const *const *aaafMaskSource=NULL,  //!< ignore voxels in source where mask==0
               Scalar const *const *const *aaafMaskDest=NULL,  //!< don't cast votes wherever mask==0
               bool detect_curves_not_surfaces=false,
               //Scalar saliency_threshold = 0.0,
               bool normalize=true,
               ostream *pReportProgress=NULL  //!< print progress to the user?
               )
  {
    if (pReportProgress)
      *pReportProgress
        << " -- Attempting to allocate space for one more image.\n"
        << " -- (If this crashes your computer, find a computer with\n"
        << " --  more RAM and use \"ulimit\", OR use a smaller image.)\n";

    float *afDenominator = NULL;
    float ***aaafDenominator = NULL;

    if (normalize) {
      Alloc3D(image_size,
              &afDenominator,
              &aaafDenominator);
    }

    TVDenseStick(image_size,
                 aaafSaliency,
                 aaaafV,
                 aaaafDest,
                 aaafMaskSource,
                 aaafMaskDest,
                 detect_curves_not_surfaces,
                 //saliency_threshold,
                 aaafDenominator,
                 pReportProgress);


    // If any of the tensor voting sums were incomplete
    // (due to image boundaries, or mask boundaries).
    // then normalize the resulting magnitudes of the filter
    // at pixels located near the image boundaries (or mask boundaries).
    // OPTIONAL (probably not useful to most users)
    if (normalize) {
      if (pReportProgress)
        *pReportProgress << "  Normalizing the result of tensor voting" << endl;

      // REMOVE THIS CRUFT:
      //if (aaafMaskSource) {

      assert(aaafDenominator);
      for (Integer iz=0; iz<image_size[2]; iz++) {
        #pragma omp parallel for collapse(2)
        for (Integer iy=0; iy<image_size[1]; iy++) {
          for (Integer ix=0; ix<image_size[0]; ix++) {
            if ((! aaafMaskDest) || (aaafMaskDest[iz][iy][ix] == 0))
              continue;
            assert(aaafDenominator[iz][iy][ix] > 0.0);
            for (int di = 0; di < 3; di++) {
              for (int dj = di; dj < 3; dj++) {
                aaaafDest[iz][iy][ix][ MapIndices_3x3_to_linear[di][dj] ]
                  /= aaafDenominator[iz][iy][ix];
              }
            }
          }
        }
      } //for (Integer iz=0; iz<image_size[2]; iz++)
      Dealloc3D(image_size,
                &afDenominator,
                &aaafDenominator);

      // REMOVE THIS CRUFT:
      //} // if (aaafMask)
      //else {
      //  // When no mask is supplied, 
      //  // If there is no mask, but the user wants the result to be normalized,
      //  // then we convolve the filter with the rectangular box.  This is cheaper
      //  // because the convolution of a separable filter with a rectangular box 
      //  // shaped function is the product of the convolution with three 1-D 
      //  // functions which are 1 from 0..image_size[d], and 0 everywhere else.
      //  Filter1D<Scalar, int> filter1d = 
      //    GenFilterGauss1D(sigma, truncate_halfwidth);
      //  Scalar *aafDenom_precomputed[3];
      //  for (int d=0; d<3; d++) {
      //    Scalar *afAllOnes = new Scalar [image_size[d]];
      //    aafDenom_precomputed[d] = new Scalar [image_size[d]];
      //    for (Integer i=0; i < image_size[d]; i++)
      //      afAllOnes[i] = 1.0;
      //    filter1d.Apply(image_size[d], afAllOnes, aafDenom_precomputed[d]);
      //    delete [] afAllOnes;
      //  }
      //  for (Integer iz = 0; iz < image_size[2]; iz++) {
      //    #pragma omp parallel for collapse(2)
      //    for (Integer iy = 0; iy < image_size[1]; iy++) {
      //      for (Integer ix = 0; ix < image_size[0]; ix++) {
      //        if ((! aaafMaskDest) || (aaafMaskDest[iz][iy][ix] == 0))
      //          continue;
      //        Scalar denominator = (aafDenom_precomputed[0][ix] *
      //                              aafDenom_precomputed[1][iy] *
      //                              aafDenom_precomputed[2][iz]);
      //        for (int di=0; di<3; di++) {
      //          for (int dj=0; dj<3; dj++) {
      //            aaaafDest[iz][iy][ix][ MapIndices_3x3_to_linear[di][dj] ]
      //              /= denominator;
      //          }
      //        }
      //      }
      //    }
      //  }
      //  // delete the array we created for storing the precomputed denominator:
      //  for (int d=0; d<3; d++)
      //    delete [] aafDenom_precomputed[d];
      //} // else clause for "if (aaafMaskSource)"

    } //if (normalize)

  } // TVDenseStick()


  void
  TVDenseStick(Integer const image_size[3],  //!< source image size
               Scalar const *const *const *aaafSaliency,  //!< saliency (score) of each voxel (usually based on Hessian eigenvalues)
               VectorContainer const *const *const *aaaafV,  //!< vector associated with each voxel
               TensorContainer ***aaaafDest,  //!< votes will be collected here
               Scalar const *const *const *aaafMaskSource=NULL,  //!< ignore voxels in source where mask==0
               Scalar const *const *const *aaafMaskDest=NULL,  //!< don't cast votes wherever mask==0
               bool detect_curves_not_surfaces=false,
               //Scalar saliency_threshold = 0.0,
               Scalar ***aaafDenominator=NULL,
               ostream *pReportProgress=NULL  //!< print progress to the user?
               )
  {
    assert(aaafSaliency);
    assert(aaaafV);
    assert(aaaafDest);


    //optional: count the number of voxels which can
    //          cast votes (useful for benchmarking)
    if (pReportProgress) {
      size_t n_salient = 0;
      size_t n_all = 0;
      for (Integer iz=0; iz<image_size[2]; iz++) {
        for (Integer iy=0; iy<image_size[1]; iy++) {
          for (Integer ix=0; ix<image_size[0]; ix++) {
            n_all++;
            //if (aaafSaliency[iz][iy][ix] > saliency_threshold)
            if (aaafSaliency[iz][iy][ix] != 0.0)
              n_salient++;
          }
        }
      }
      *pReportProgress << "  (fraction of salient voxels = "
                       << (static_cast<double>(n_salient) / n_all)
                       << ")\n"
                       << "  (use aggressive thresholding to improve speed)"
                       << endl;
    }


    // First, initialize the arrays which will store the results with zeros.
    for (Integer iz=0; iz<image_size[2]; iz++)
      for (Integer iy=0; iy<image_size[1]; iy++)
        for (Integer ix=0; ix<image_size[0]; ix++)
          for (int di=0; di<3; di++)
            for (int dj=0; dj<3; dj++)
              aaaafDest[iz][iy][ix][ MapIndices_3x3_to_linear[di][dj] ] = 0.0;

    if (aaafDenominator) {
      // aaafDenominator[][][] keeps track of how much of the sum of
      // tensor-voting contributions was available at each voxel location.
      // (If some voxels were unavailable or outside the boundaries of 
      //  the image, they cannot cast votes at this voxel location.)
      //  This array keeps track of that.)  We should initialize this array too.
      for (Integer iz=0; iz<image_size[2]; iz++)
        for (Integer iy=0; iy<image_size[1]; iy++)
          for (Integer ix=0; ix<image_size[0]; ix++)
            aaafDenominator[iz][iy][ix] = 0.0;
    }

    // REMOVE THIS CRUFT
    //assert(pv);
    //assert(pV->nchannels() == 3);
    //pV->Resize(image_size, aaafMaskSource, pReportProgress);

    if (pReportProgress)
      *pReportProgress << "---- Begin Tensor Voting (dense, stick) ----\n"
                       << "  progress: processing plane#" << endl;

    // The mask should be 1 everywhere we want to consider, and 0 elsewhere.
    // Multiplying the density in the tomogram by the mask removes some of 
    // the voxels from consideration later on when we do the filtering.
    // (Later, we will adjust the weight of the average we compute when we
    //  apply the filter in order to account for the voxels we deleted now.)

    for (Integer iz=0; iz<image_size[2]; iz++) {

      if (pReportProgress)
        *pReportProgress << "  " << iz+1 << " / " << image_size[2] << "\n";

      #pragma omp parallel for collapse(2)
      for (Integer iy=0; iy<image_size[1]; iy++) {

        for (Integer ix=0; ix<image_size[0]; ix++) {

          // ------------ Uncomment whichever is faster: -----------
          // EITHER uncomment TVCastStickVotes() or TVReceiveStickVotes()
          // (VERSION1 OR VERSION2), BUT NOT BOTH

          // VERSION1:
          //           Have the voxel at ix,iy,iz cast votes at nearby voxels:

          TVCastStickVotes(ix, iy, iz,
                           image_size,
                           aaafSaliency,
                           aaaafV,
                           aaaafDest,
                           aaafMaskSource,
                           aaafMaskDest,
                           detect_curves_not_surfaces,
                           //saliency_threshold,
                           aaafDenominator);

          // VERSION 2: Have the voxel at ix,iy,iz receive votes
          //            from nearby voxels
          //
          //TVReceiveStickVotes(ix, iy, iz,
          //                    image_size,
          //                    aaafSaliency,
          //                    aaaafV,
          //                    aaaafDest,
          //                    aaafMaskSource,
          //                    aaafMaskDest,
          //                    detect_curves_not_surfaces,
          //                    //saliency_threshold,
          //                    (aaafDenominator
          //                     ? &(aaafDenominator[iz][iy][ix])
          //                     : NULL));

        }
      }
    }

    if (pReportProgress)
      *pReportProgress << "---- Diagonalizing Tensor Voting results ----" << endl;

    // Diagonalize the resulting tensor.
    // The resulting eigenvalues and eigenvectors can be analyzed by the caller.
    DiagonalizeHessianImage3D(image_size,
                              aaaafDest, //<--undiagonalized tensors (input)
                              aaaafDest, //<--diagonalized tensors (output)
                              aaafMaskDest,
                              selfadjoint_eigen3::DECREASING_EIVALS,
                              pReportProgress);


    // DEBUG: assert that all eigenvalues are nonnegative
    for (Integer iz=0; iz<image_size[2]; iz++) {
      for (Integer iy=0; iy<image_size[1]; iy++) {
        for (Integer ix=0; ix<image_size[0]; ix++) {
          if ((! aaafMaskDest) || (aaafMaskDest[iz][iy][ix] == 0.0))
            continue;
          // The eigenvalues are in the first 3 entries of aaaafDest[iz][iy][ix]
          float *eivals = aaaafDest[iz][iy][ix];
          assert(eivals);
          for (int d=0; d<3; d++)
            assert(eivals[d] >= 0.0);
        }
      }
    }

  } //TVDenseStick()


private:
  void
  TVCastStickVotes(Integer ix,  //!< coordinates of the voter
                   Integer iy,  //!< coordinates of the voter
                   Integer iz,  //!< coordinates of the voter
                   Integer const image_size[3],
                   Scalar const *const *const *aaafSaliency, //!< saliency (score) of each voxel (usually calculated from Hessian eigenvalues)
                   VectorContainer const *const *const *aaaafV,  //!< vector associated with each voxel
                   TensorContainer ***aaaafDest,  //!< votes will be collected here
                   Scalar const *const *const *aaafMaskSource,  //!< ignore voxels in source where mask==0
                   Scalar const *const *const *aaafMaskDest,  //!< ignore voxels in dest where mask==0
                   bool detect_curves_not_surfaces = false,
                   //Scalar saliency_threshold = 0.0,
                   Scalar ***aaafDenominator = NULL) const
  {
    assert(aaafSaliency);
    assert(aaaafV);
    assert(aaaafDest);

    Scalar saliency = aaafSaliency[iz][iy][ix];
    if (saliency == 0.0)
    //if (saliency <= saliency_threshold)
      return;

    Scalar mask_val = 1.0;
    if (aaafMaskSource) {
      mask_val = aaafMaskSource[iz][iy][ix];
      if (mask_val == 0.0)
        return;
    }

    Scalar n[3]; // the direction of the stick tensor (see below)
    for (int d=0; d<3; d++)
      n[d] = aaaafV[iz][iy][ix][d];

    //Note: The "filter_val" also is needed to calculate
    //      the denominator used in normalization.
    //      It is unusual to use a mask unless you intend
    //      to normalize the result later, but I don't enforce this

    for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++) {
      Integer iz_jz = iz+jz;
      if ((iz_jz < 0) || (image_size[2] <= iz_jz))
        continue;

      for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {
        Integer iy_jy = iy+jy;
        if ((iy_jy < 0) || (image_size[1] <= iy_jy))
          continue;

        for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {
          Integer ix_jx = ix+jx;
          if ((ix_jx < 0) || (image_size[0] <= ix_jx))
            continue;

          if ((aaafMaskDest) && (aaafMaskDest[iz_jz][iy_jy][ix_jx] == 0.0))
            continue;

          // The function describing how the vote-strength falls off with
          // distance has been precomputed and is stored in
          // radial_decay_lookup.aaafH[jz][jy][jx];
          // In most tensor-voting implementations, this is a Gaussian.
          Scalar filter_val = radial_decay_lookup.aaafH[jz][jy][jx];
          if (aaafMaskSource)
            filter_val *= mask_val;

          Scalar decay_radial = filter_val;
          if (decay_radial == 0.0)
            continue;

          Scalar r[3];
          for (int d=0; d<3; d++)
            r[d] = aaaafDisplacement[jz][jy][jx][d];

          //
          //   .
          //   :\
          //   : \
          // ->:  \<- theta
          //   :   \
          //   :    \ 
          //   :     \
          //   :      \
          //   :       \  
          //   ^        \
          //   |         \ 
          // n |          \
          //   |        .-' r
          //   |,,..--''
          //  0
          //
          // "n" = the direction of the stick tensor.  If the object being 
          //       detected is a surface, n is perpendicular to that surface.
          //       If the object is a curve, then n points tangent to the curve.
          // "r" = the position of the vote receiver relative to the voter
          //
          // "theta" = the angle of r relative to the plane perpendicular to n.
          //           (NOT the angle of r relative to n.  This is a confusing
          //            convention, but this is how it is normally defined.)

          Scalar sintheta = DotProduct3(r, n); //(sin() not cos(), see diagram)
          Scalar sinx2 = sintheta * 2.0;
          Scalar sin2 = sintheta * sintheta;
          Scalar cos2 = 1.0 - sin2;
          Scalar angle_dependence2 = cos2;
          if (detect_curves_not_surfaces) {
            // If we are detecting 1-dimensional curves (polymers, etc...)
            // instead of 2-dimensional surfaces (membranes, ...)
            // then the stick direction is assumed to be along the 
            // of the curve, (as opposed to perpendicular to the 2D surface).
            // As such
            angle_dependence2 = sin2;
          }

          Scalar decay_angular;

          switch(exponent) {
          case 2:
            decay_angular = angle_dependence2;
            break;
          case 4:
            decay_angular = angle_dependence2 * angle_dependence2;
            break;
          default:
            //Scalar angle_dependence = sqrt(angle_dependence2);
            decay_angular = pow(angle_dependence2, 0.5*exponent);
            break;
          }

          Scalar n_rotated[3];
          for (int d=0; d<3; d++) {
            if (detect_curves_not_surfaces)
              n_rotated[d] = n[d] - sinx2*r[d];
            else
              n_rotated[d] = sinx2*r[d] - n[d];
          }

          Scalar tensor_vote[3][3];
          for (int di = 0; di < 3; di++) {
            for (int dj = di; dj < 3; dj++) {
              tensor_vote[di][dj] = (saliency *
                                     decay_radial *
                                     decay_angular *
                                     n_rotated[di] * n_rotated[dj]);

              // OLD CODE: I used to implement aaaafDest as a 5-D array.
              //
              //aaaafDest[iz_jz][iy_jy][ix_jx][di][dj] += tensor_vote[di][dj];
              //
              // NEW CODE:
              // The ix_jx,iy_jy,iz_jz'th entry in aaaafDest should be
              // a 3x3 matrix. Since this matrix is symmentric, it contains
              // only 6 non-redundant entries.  Consequently there is no
              // need to store 9 numbers if only 6 are needed.
              // So I implemented a version of the 3x3 matrix which has
              // only 6 entries, arranged in a 1-D array of size 6.
              // To access these entries, use "MapIndices_3x3_to_linear[][]".

              aaaafDest[iz_jz][iy_jy][ix_jx][MapIndices_3x3_to_linear[di][dj]]
                += tensor_vote[di][dj];
            }
          }

          if (aaafDenominator)
            aaafDenominator[iz_jz][iy_jy][ix_jx] += filter_val;

        } // for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++)
      } // for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++)
    } // for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++)


  } // TVCastStickVotes()




  void
  TVReceiveStickVotes(Integer ix,  //!< coordinates of the receiver voxel
                      Integer iy,  //!< coordinates of the receiver voxel
                      Integer iz,  //!< coordinates of the receiver voxel
                      Integer const image_size[3],
                      Scalar const *const *const *aaafSaliency, //!< saliency (score) of each voxel (usually calculated from Hessian eigenvalues)
                      VectorContainer const *const *const *aaaafV,  //!< vector associated with each voxel
                      TensorContainer ***aaaafDest,  //!< votes will be collected here
                      Scalar const *const *const *aaafMaskSource,  //!< ignore voxels in source where mask==0
                      Scalar const *const *const *aaafMaskDest,  //!< ignore voxels in dest where mask==0
                      bool detect_curves_not_surfaces = false,
                      //Scalar saliency_threshold = 0.0,
                      Scalar *pDenominator = NULL) const
  {
    assert(aaafSaliency);
    assert(aaaafV);
    assert(aaaafDest);

    if (aaafMaskDest && (aaafMaskDest[iz][iy][ix] == 0.0))
      return;

    Scalar denominator = 0.0;

    for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++) {
      Integer iz_jz = iz-jz;
      if ((iz_jz < 0) || (image_size[2] <= iz_jz))
        continue;

      for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {
        Integer iy_jy = iy-jy;
        if ((iy_jy < 0) || (image_size[1] <= iy_jy))
          continue;

        for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {
          Integer ix_jx = ix-jx;
          if ((ix_jx < 0) || (image_size[0] <= ix_jx))
            continue;

          // The function describing how the vote-strength falls off with
          // distance has been precomputed and is stored in
          // radial_decay_lookup.aaafH[jz][jy][jx];
          // In most tensor-voting implementations, this is a Gaussian.
          Scalar filter_val = radial_decay_lookup.aaafH[jz][jy][jx];

          if (aaafMaskSource) {
            Scalar mask_val = aaafMaskSource[iz_jz][iy_jy][ix_jx];
            if (mask_val == 0.0)
              continue;
            filter_val *= mask_val;
          }
          //Note: The "filter_val" also is needed to calculate
          //      the denominator used in normalization.
          //      It is unusual to use a mask unless you intend
          //      to normalize the result later, but I don't enforce this

          Scalar saliency = aaafSaliency[iz_jz][iy_jy][ix_jx];
          //if (saliency <= saliency_threshold)
          if (saliency == 0.0)
            continue;

          Scalar decay_radial = filter_val;
          if (decay_radial == 0.0)
            continue;

          Scalar r[3];
          Scalar n[3];
          for (int d=0; d<3; d++) {
            r[d] = aaaafDisplacement[jz][jy][jx][d];
            n[d] = aaaafV[iz_jz][iy_jy][ix_jx][d];
          }

          //
          //   .
          //   :\
          //   : \
          // ->:  \<- theta
          //   :   \
          //   :    \ 
          //   :     \
          //   :      \
          //   :       \  
          //   ^        \
          //   |         \ 
          // n |          \
          //   |        .-' r
          //   |,,..--''
          //  0
          //
          // "n" = the direction of the stick tensor.  If the object being 
          //       detected is a surface, n is perpendicular to that surface.
          //       If the object is a curve, then n points tangent to the curve.
          // "r" = the position of the vote receiver relative to the voter
          //
          // "theta" = the angle of r relative to the plane perpendicular to n.
          //           (NOT the angle of r relative to n.  This is a confusing
          //            convention, but this is how it is normally defined.)

          Scalar sintheta = DotProduct3(r, n); //(sin() not cos(), see diagram)
          Scalar sinx2 = sintheta * 2.0;
          Scalar sin2 = sintheta * sintheta;
          Scalar cos2 = 1.0 - sin2;
          Scalar angle_dependence2 = cos2;
          if (detect_curves_not_surfaces) {
            // If we are detecting 1-dimensional curves (polymers, etc...)
            // instead of 2-dimensional surfaces (membranes, ...)
            // then the stick direction is assumed to be along the 
            // of the curve, (as opposed to perpendicular to the 2D surface).
            // As such
            angle_dependence2 = sin2;
          }

          Scalar decay_angular;

          switch(exponent) {
          case 2:
            decay_angular = angle_dependence2;
            break;
          case 4:
            decay_angular = angle_dependence2 * angle_dependence2;
            break;
          default:
            //Scalar angle_dependence = sqrt(angle_dependence2);
            decay_angular = pow(angle_dependence2, 0.5*exponent);
            break;
          }

          Scalar n_rotated[3];
          for (int d=0; d<3; d++) {
            if (detect_curves_not_surfaces)
              n_rotated[d] = n[d] - sinx2*r[d];
            else
              n_rotated[d] = sinx2*r[d] - n[d];
          }

          Scalar tensor_vote[3][3];
          for (int di = 0; di < 3; di++) {
            for (int dj = di; dj < 3; dj++) {
              tensor_vote[di][dj] = (saliency *
                                     decay_radial *
                                     decay_angular *
                                     n_rotated[di] * n_rotated[dj]);

              // OLD CODE: I used to implement aaaafDest as a 5-D array.
              //
              //aaaafDest[iz][iy][ix][di][dj] += tensor_vote[di][dj];
              //
              // NEW CODE:
              // The ix,iy,iz'th entry in aaaafDest should be a 3x3 matrix.
              // Since this matrix is symmentric, it contains only 6
              // non-redundant entries.  Consequently there is no need to
              // store 9 numbers if only 6 are needed.
              // So I implemented a version of the 3x3 matrix which has
              // only 6 entries, arranged in a 1-D array of size 6.
              // To access these entries, use "MapIndices_3x3_to_linear[][]".

              aaaafDest[iz][iy][ix][ MapIndices_3x3_to_linear[di][dj] ]
                += tensor_vote[di][dj];
            }
          }

          if (pDenominator)
            denominator += filter_val;
        } // for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++)
      } // for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++)
    } // for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++)

    if (pDenominator)
      *pDenominator = denominator;

  } // TVReceiveStickVotes()



  inline void swap(TV3D<Scalar,Integer,VectorContainer,TensorContainer> &other)
  {
    std::swap(sigma, other.sigma);
    std::swap(exponent, other.exponent);
    std::swap(radial_decay_lookup, other.radial_decay_lookup);
    std::swap(aafDisplacement, other.aafDisplacement);
    std::swap(aaaafDisplacement, other.aaaafDisplacement);
    std::swap(halfwidth, other.halfwidth);
    std::swap(array_size, other.array_size);
  }


  inline TV3D<Scalar,Integer,VectorContainer,TensorContainer>&
    operator = (TV3D<Scalar,Integer,VectorContainer,TensorContainer> source)
  {
    this->swap(source);
    return *this;
  }

private:

  void Init() {
    sigma = 0.0;
    halfwidth[0] = -1;
    halfwidth[1] = -1;
    halfwidth[2] = -1;
    array_size[0] = -1;
    array_size[1] = -1;
    array_size[2] = -1;
    aafDisplacement = NULL;
    aaaafDisplacement = NULL;
  }


  void Resize(Integer set_halfwidth[3]) {
    for (int d=0; d<3; d++) {
      array_size[d] = 2*halfwidth[d] + 1;
    }

    Scalar sigmas[3] = {sigma, sigma, sigma};
    radial_decay_lookup =
      GenFilterGenGauss3D(sigmas,
                          static_cast<float>(2.0),
                          halfwidth);

    AllocDisplacement();
    PrecalcDisplacement(aaaafDisplacement);
  }


  void DeallocDisplacement() {
    if (aaaafDisplacement) {
      for (int d=0; d<3; d++) {
        array_size[d] = 2*halfwidth[d] + 1;
      }
      //shift pointers back to normal
      for (Integer iz = -halfwidth[2]; iz <= halfwidth[2]; iz++) {
        for (Integer iy = -halfwidth[1]; iy <= halfwidth[1]; iy++) {
          aaaafDisplacement[iz][iy] -= halfwidth[0];
        }
        aaaafDisplacement[iz] -= halfwidth[1];
      }
      aaaafDisplacement -= halfwidth[2];
      Dealloc3D(array_size,
                &aafDisplacement,
                &aaaafDisplacement);
    }
  }

  void AllocDisplacement() {
    if (aaaafDisplacement)
      Dealloc3D(array_size,
                &aafDisplacement,
                &aaaafDisplacement);
    Alloc3D(array_size,
            &aafDisplacement,
            &aaaafDisplacement);
    //shift pointers to enable indexing from i = -halfwidth .. +halfwidth
    for (Integer iz = 0; iz < array_size[2]; iz++) {
      for (Integer iy = 0; iy < array_size[1]; iy++) {
        aaaafDisplacement[iz][iy] += halfwidth[0];
      }
      aaaafDisplacement[iz] += halfwidth[1];
    }
    aaaafDisplacement += halfwidth[2];
  }

  void PrecalcDisplacement(array<Scalar, 3> ***aaaafDisplacement) {
    // pre-compute the normalized radius unit vector
    for (Integer iz = -halfwidth[2]; iz <= halfwidth[2]; iz++) {
      for (Integer iy = -halfwidth[1]; iy <= halfwidth[1]; iy++) {
        for (Integer ix = -halfwidth[0]; ix <= halfwidth[0]; ix++) {
          Scalar length = sqrt(ix*ix + iy*iy + iz*iz);
          if (length == 0)
            length = 1.0;
          aaaafDisplacement[iz][iy][ix][0] = ix / length;
          aaaafDisplacement[iz][iy][ix][1] = iy / length;
          aaaafDisplacement[iz][iy][ix][2] = iz / length;
        }
      }
    }
  }

}; // class TV3D


#endif //#ifndef _FILTER3D_H
