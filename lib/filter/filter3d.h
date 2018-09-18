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
using namespace std;
#include <filter1d.h>  // defines "Filter1D" (used in "ApplyGauss()")
#include <filter3d_utils.h> // defines "AverageArr()", and similar functions...







/// @brief A simple class for general linear convolutional filters in 3D
///
/// @note  In practice, this class is not used often because separable filters
///        based on Gaussians are much much faster.
///        -A 2018-9-11

template<class RealNum, class Integer>

class Filter3D {
public:
  RealNum *afH;         //!< contiguous block of memory storing the filter array
  RealNum ***aaafH;     //!< the same array which can be indexed using [i][j][k] notation
  Integer halfwidth[3]; //!< num pixels from filter center to edge in x,y,z directions
  Integer array_size[3]; //!<size of the array in x,y,z directions (2*halfwidth+1)

 public:

  /// @brief  Apply the filter to a 3D image (aaafSource[][][]).
  /// Save the results in the "aafDest" array.  (A "mask" is optional.)
  /// All arrays are 3D and assumed to be the same size.
  ///     
  /// @code
  /// If aafMask == NULL, then filter computes g(i):
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
  /// where: f(i) is the original image at position i (ix,iy,iz)
  ///          i  is shorthand for ix,iy,iz
  ///        g(i) is the image after the filter has been applied
  ///        h(j) is the filter
  ///        mask(i) selects the pixels we care about (usually either 0 or 1)
  /// @endcode
  ///
  /// @param size_source contains size of the source image (in the x,y,z directions)
  /// @param aaafSource[][][] is the source array (source image) <==> f(i)
  /// @param aaafDest[][][] will store the image after filtering <==> g(i)
  /// @param aaafMask[][][]==0 whenever we want to ignore entries in afSource[]. Optional.
  /// @param normalize  This boolean parameter = true if you want to divide g(i) by the sum of the weights considered. Optional.
  /// (useful if the sum of your filter elements, h(j), is 1, and if the sum was
  ///  not complete because some entries lie outside the mask or the boundary.)

  void Apply(Integer const size_source[3],
             RealNum const *const *const *aaafSource,
             RealNum ***aaafDest,
             RealNum const *const *const *aaafMask = NULL,
             bool normalize = false,
             ostream *pReportProgress = NULL) const
  {
    RealNum *afDenominator = NULL;
    RealNum ***aaafDenominator = NULL;
    if (normalize)
      Alloc3D(size_source, &afDenominator, &aaafDenominator);

    Apply(size_source,
          aaafSource,
          aaafDest,
          aaafMask,
          aaafDenominator,
          pReportProgress);

    if (normalize) {
      for (int iz=0; iz < size_source[2]; iz++)
        for (int iy=0; iy < size_source[1]; iy++)
          for (int ix=0; ix < size_source[0]; ix++)
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
  ///          i  is shorthand for ix,iy,iz
  ///        g(i) is the data after the filter has been applied
  ///        h(j) is the filter
  ///        mask(i) is usually either 0 or 1
  ///        d(i) is the "denominator" = sum of the filter weights after masking
  /// @endcode
  ///
  /// @param size_source contains size of the source image (in the x,y,z directions)
  /// @param aaafSource[][][] is the source array (source image) <==> g(i)
  /// @param aaafDest[][][] will store the image after filtering <==> g(i)
  /// @param aaafMask[][][]==0 whenever we want to ignore entries in afSource[][]. Optional.
  /// @param aaafDenominator[][][] will store d(i) if you supply a non-NULL pointer

  void Apply(Integer const size_source[3],
             RealNum const *const *const *aaafSource,
             RealNum ***aaafDest,
             RealNum const *const *const *aaafMask = NULL,
             RealNum ***aaafDenominator = NULL,
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



  inline Filter3D(const Filter3D<RealNum, Integer>& source) {
    Init();
    Resize(source.halfwidth); // allocates and initializes afH and aaafH
    //for(Int iz=-halfwidth[2]; iz<=halfwidth[2]; iz++)
    //  for(Int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
    //    for(Int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
    //      aaafH[iz][iy][ix] = source.aaafH[iz][iy][ix];
    // -- Use memcpy() instead: --
    //memcpy(afH,
    //       source.afH,
    //       array_size[0] * array_size[1] * array_size[2],
    //       *sizeof(RealNum));
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


  inline void swap(Filter3D<RealNum, Integer> &other) {
    std::swap(afH, other.afH);
    std::swap(aaafH, other.aaafH);
    std::swap(halfwidth, other.halfwidth);
    std::swap(array_size, other.array_size);
  }


  inline Filter3D<RealNum, Integer>&
    operator = (Filter3D<RealNum, Integer> source) {
    this->swap(source);
    return *this;
  }


  /// @ brief   Make sure the sum of the filter weights (in aaafH) is 1
  void Normalize() {
    RealNum total = 0.0;
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
  RealNum Average(RealNum const *const *const *aaafW=NULL) const {
    return AverageArr(array_size, aaafH, aaafW);
  }

  /// @brief Calculate the (weighted) average squared values in the filter array, aaafH
  /// @param aaafW optional weights used when calculating the averages
  /// @return the (weighted) average squared values in the filter array aaafH
  RealNum AverageSqr(RealNum const *const *const *aaafW=NULL) const {
    return _AveSqrArr(array_size, aaafH, aaafW);
  }

  /// @brief Calculate the (weighted) standard deviation of the filter array values
  /// @param aaafW optional weights used when calculating the standard deviation
  /// @return the (weighted) standard deviation of the filter array aaafH
  RealNum StdDev(RealNum const *const *const *aaafW=NULL) const {
    return StdDevArr(array_size, aaafH, aaafW);
  }

  /// @brief Calculate the (weighted) sum of the filter array values, aaafH
  /// @param aaafW optional weights used when calculating the sum
  /// @return the (weighted) sum of the filter array values
  RealNum Sum(RealNum const *const *const *aaafW=NULL) const {
    return _SumArr(array_size, aaafH, aaafW);
  }

  /// @brief Calculate the (weighted) sum of the squared filter array values
  /// @param aaafW optional weights used when calculating the sum
  /// @return the (weighted) sum of the squared filter array values
  RealNum SumSqr(RealNum const *const *const *aaafW=NULL) const {
    return _SumSqrArr(array_size, aaafH, aaafW);
  }

  /// @brief Add a number to all of the filter array values, aaafH
  /// @param offset  the number to add
  void AddScalar(RealNum offset) {
    _AddScalarArr(offset, array_size, aaafH);
  }

  /// @brief multiply all of the filter array values (aaafH) by a number
  /// @param offset  the number to multiply
  void MultiplyScalar(RealNum scale) {
    _MultiplyScalarArr(scale, array_size, aaafH);
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
  RealNum ApplyToVoxel(Integer ix,
                       Integer iy,
                       Integer iz,
                       Integer const size_source[3],
                       RealNum const *const *const *aaafSource,
                       RealNum const *const *const *aaafMask = NULL,
                       RealNum *pDenominator = NULL) const

  {
    RealNum g = 0.0;
    RealNum denominator = 0.0;

    for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++) {
      int iz_jz = iz-jz;
      if ((iz_jz < 0) || (size_source[2] <= iz_jz))
        continue;

      for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {
        int iy_jy = iy-jy;
        if ((iy_jy < 0) || (size_source[1] <= iy_jy))
          continue;

        for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {
          int ix_jx = ix-jx;
          if ((ix_jx < 0) || (size_source[0] <= ix_jx))
            continue;

          RealNum filter_val = aaafH[jz][jy][jx];

          if (aaafMask)
              filter_val *= aaafMask[iz_jz][iy_jy][ix_jx];
              //Note: The "filter_val" also is needed to calculate
              //      the denominator used in normalization.
              //      It is unusual to use a mask unless you intend
              //      to normalize the result later, but I don't enforce this

          RealNum delta_g = 
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
    for(Integer d=0; d < 3; d++) {
      halfwidth[d] = set_halfwidth[d];
      array_size[d] = 1 + 2*halfwidth[d];
    }
    Alloc3D(array_size, &afH, &aaafH);
    for (int iz = 0; iz < array_size[2]; iz++)
      for (int iy = 0; iy < array_size[1]; iy++)
        for (int ix = 0; ix < array_size[0]; ix++)
          aaafH[iz][iy][ix] = -1.0e38; //(if uninitiliazed memory read, we will know)

    //shift pointers to enable indexing from i = -halfwidth .. +halfwidth
    aaafH += halfwidth[2];
    for (int iz = 0; iz < array_size[2]; iz++) {
      aaafH[iz] += halfwidth[1];
      for (int iy = 0; iy < array_size[1]; iy++) {
        aaafH[iz][iy] += halfwidth[0];
      }
    }
  }


  /// @brief allocate space used by the filter array
  void Dealloc() {
    //Integer array_size[3];
    for(Integer d=0; d < 3; d++) {
      array_size[d] = 1 + 2*halfwidth[d];
      halfwidth[d] = -1;
    }
    //shift pointers back to normal
    aaafH -= halfwidth[2];
    for (int iz = 0; iz < array_size[2]; iz++) {
      aaafH[iz] -= halfwidth[1];
      for (int iy = 0; iy < array_size[1]; iy++) {
        aaafH[iz][iy] -= halfwidth[0];
      }
    }
    //then deallocate
    Dealloc3D(array_size, &afH, &aaafH);
    for(Integer d=0; d < 3; d++)
      array_size[d] = -1;
  }


  void Resize(Integer const set_halfwidth[3]) {
    Dealloc();
    Alloc(set_halfwidth);
  }


  void Init() {
    halfwidth[0] = -1;
    halfwidth[1] = -1;
    halfwidth[2] = -1;
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
  ///         This function was not intended for public use.
  /// @return the weighted sum of the error (raised to the power of
  ///         err_exponent, which is 2 by default) summed over all entries
  ///         in the template.  Weights are typically chosen so that they
  ///         decay to 0 near the boundaries of the aaafH array.
  ///         This function is not yet intended for public use.
  RealNum
  _TemplateError(Integer ix, //!< voxel's position in the x,y,z directions
                 Integer iy, //!< voxel's position in the x,y,z directions
                 Integer iz, //!< voxel's position in the x,y,z directions
                 Integer const size_source[3], //!< size of the source array
                 RealNum const *const *const *aaafSource,   //!< the array 
                 RealNum scale, //!< how much to scale template intensities
                 RealNum const *const *const *aaafW, //!< weights used in all summations
                 RealNum err_exponent=2, //!< exponent used for calculating error
                 RealNum const *const *const *aaafMask = NULL, //!< optional array storing voxels we should ignore (if 0)
                 ostream *pReportProgress = NULL //!< optional ostream for printing out progress of the calculation
                 ) const 
  {
      RealNum g = 0.0;
      RealNum denominator = 0.0;

      for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++) {
      int iz_jz = iz-jz;
      if ((iz_jz < 0) || (size_source[2] <= iz_jz))
        continue;

      for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {
        int iy_jy = iy-jy;
        if ((iy_jy < 0) || (size_source[1] <= iy_jy))
          continue;

        for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {
          int ix_jx = ix-jx;
          if ((ix_jx < 0) || (size_source[0] <= ix_jx))
            continue;

          RealNum delta_g = 
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
                     RealNum const *const *const *aaafSource, //!< source array
                     RealNum ***aaafDest, //!< store the template error results here
                     RealNum ***aaafC, //!< how much to scale template intensities
                     RealNum const *const *const *aaafW, //!< weights used in all summations
                     RealNum err_exponent=2, //!< exponent used for calculating error
                     RealNum const *const *const *aaafMask = NULL, //!< optional: indicate which entries should be ignored
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

}; // class Filter






/// @brief Create a 3D filter and fill it with a "generalized Gaussian" function
///   @verbatim h(x,y,z) = A*exp(-r^m)  @endverbatim
/// where  @verbatim r  = sqrt((x/σ_x)^2 + (y/σ_y)^2 + (z/σ_z)^2) @endverbatim
///   and  "A" is determined by normalization of the discrete sum
/// @note "A" is equal to the value stored in the middle of the array,
///       The caller can determine what "A" is by looking at this value.

template<class RealNum>
Filter3D<RealNum, int>
GenFilterGenGauss3D(RealNum width[3],    //!< "σ_x", "σ_y", "σ_z" parameters
                    RealNum m_exp,       //!< "m" exponent parameter
                    int truncate_halfwidth[3], //!< size of filter window
                    RealNum *pA=NULL,    //!< optional:report A coeff to user
                    ostream *pReportEquation=NULL//!< optional:report equation used to the user
                    )
{
  RealNum window_threshold = 1.0;
  for (int d=0; d<3; d++) {
    RealNum h = ((width[d]>0) ? exp(-pow(truncate_halfwidth[d]/width[d], m_exp)) : 1.0);
    if (h < window_threshold)
      window_threshold = h;
  }
  Filter3D<RealNum, int> filter(truncate_halfwidth);
  RealNum total = 0;
  for (int iz=-filter.halfwidth[2]; iz<=filter.halfwidth[2]; iz++) {
    for (int iy=-filter.halfwidth[1]; iy<=filter.halfwidth[1]; iy++) {
      for (int ix=-filter.halfwidth[0]; ix<=filter.halfwidth[0]; ix++) {
        RealNum r = sqrt(SQR(ix/width[0])+SQR(iy/width[1])+SQR(iz/width[2]));
        RealNum h = ((r>0) ? exp(-pow(r, m_exp)) : 1.0);
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
  RealNum A = filter.aaafH[0][0][0];


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

template<class RealNum>
Filter3D<RealNum, int>
GenFilterGenGauss3D(RealNum width[3],    //!< "σ_x", "σ_y", "σ_z" parameters
                    RealNum m_exp,       //!< "m" parameter in formula
                    RealNum filter_cutoff_ratio=2.5, //!< how many sigma (σ) before truncating?
                    RealNum *pA=NULL,    //!< optional:report A coeff to user
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
/// This function was not intended for public use

template<class RealNum>
static
Filter3D<RealNum, int> 
_GenFilterDogg3D(RealNum width_a[3],  //!< "a" parameter in formula
                 RealNum width_b[3],  //!< "b" parameter in formula
                 RealNum m_exp,  //!< "m" parameter in formula
                 RealNum n_exp,  //!< "n" parameter in formula
                 Filter3D<RealNum, int>& filter_A, //!< filters for the two
                 Filter3D<RealNum, int>& filter_B, //!< gaussians
                 RealNum *pA=NULL, //!< optional:report A,B coeffs to user
                 RealNum *pB=NULL, //!< optional:report A,B coeffs to user
                 ostream *pReportEquation = NULL //!< optional: report equation to the user
                 )
{
  RealNum A, B;
  //A, B = height of the central peak
  A = filter_A.aaafH[0][0][0];
  B = filter_B.aaafH[0][0][0];


  // The "difference of gaussians" filter is the difference between
  // these two (generalized) gaussian filters.
  int halfwidth[3];
  halfwidth[0] = MAX(filter_A.halfwidth[0], filter_B.halfwidth[0]);
  halfwidth[1] = MAX(filter_A.halfwidth[1], filter_B.halfwidth[1]);
  halfwidth[2] = MAX(filter_A.halfwidth[2], filter_B.halfwidth[2]);
  Filter3D<RealNum, int> filter(halfwidth);

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

        //FOR DEBUGGING REMOVE EVENTUALLY
        if (pReportEquation)
          *pReportEquation << "GenDogg: aaafH["<<iz<<"]["<<iy<<"]["<<ix<<"] = "
                           << filter.aaafH[iz][iy][ix] << endl;

        //*pReportEquation << aaafH[iz][iy][ix];
        //                         
        //if (ix == 0) pReportEquation << "\n"; else pReportEquation << " ";

      } // for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
    } // for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
  } // for (int iz=-halfwidth[2]; iz<=halfwidth[2]; iz++) {


  //FOR DEBUGGING REMOVE EVENTUALLY
  if (pReportEquation)
    *pReportEquation << "\n";

  // COMMENTING OUT the factor of 1/(A-B):
  //A = A/(A-B);
  //B = B/(A-B);

  if (pA && pB) {
    *pA = A; // Rescale A and B numbers returned to the caller
    *pB = B; // (because we divided the array entries by (A-B) earlier)
  }

  if (pReportEquation) {
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
template<class RealNum>
Filter3D<RealNum, int> 
GenFilterDogg3D(RealNum width_a[3],   //!< "a" parameter in formula
                RealNum width_b[3],   //!< "b" parameter in formula
                RealNum m_exp,        //!< "m" parameter in formula
                RealNum n_exp,        //!< "n" parameter in formula
                int halfwidth[3],     //!< the width of the filter
                RealNum *pA=NULL,     //!< optional:report A,B coeffs to user
                RealNum *pB=NULL,     //!< optional:report A,B coeffs to user
                ostream *pReportEquation = NULL //!< optional: print params used?
                )
{
  Filter3D<RealNum, int> filter_A =
    GenFilterGenGauss3D(width_a,      //"a_x", "a_y" gaussian width parameters
                        m_exp,        //"m" exponent parameter
                        halfwidth);

  Filter3D<RealNum, int> filter_B =
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




// @brief _ApplyGauss3D applies a Gaussian filter on a 3D array.
//        It assumes separate Filter1D objects have already been created 
//        which will blur the image in each direction (x,y,z).
//        This function was not intended for public use
//        Please use "ApplyGauss3d()" instead.  (See below.)
template<class RealNum>
static
RealNum
_ApplyGauss3D(int const image_size[3], 
              RealNum const *const *const *aaafSource,
              RealNum ***aaafDest,
              RealNum const *const *const *aaafMask,
              Filter1D<RealNum, int> aFilter[3], //preallocated 1D filters
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

  // Create an array to store the denominator.
  // (The "denominator" stores the normalization used when aaafMask != NULL)
  RealNum ***aaafDenom;
  RealNum *afDenom;
  if (normalize) {
    Alloc3D(image_size, &afDenom, &aaafDenom);
    for (int iz = 0; iz < image_size[2]; iz++)
      for (int iy = 0; iy < image_size[1]; iy++)
        for (int ix = 0; ix < image_size[0]; ix++)
          aaafDenom[iz][iy][ix] = 1.0; //(default value)
  }

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
    RealNum *afDest_tmp   = new RealNum [image_size[d]];
    RealNum *afSource_tmp = new RealNum [image_size[d]];
    RealNum *afMask_tmp   = NULL;
    if (aaafMask)
      afMask_tmp = new RealNum [image_size[d]];
    RealNum *afDenom_tmp = NULL;
    if (normalize)
      afDenom_tmp = new RealNum [image_size[d]];

    #pragma omp for collapse(2)
    for (int iy = 0; iy < image_size[1]; iy++) {
      for (int ix = 0; ix < image_size[0]; ix++) {

        // copy the data we need to the temporary arrays
        for (int iz = 0; iz < image_size[2]; iz++) {
          afSource_tmp[iz] = aaafDest[iz][iy][ix];  //copy from prev aaafDest
          if (aaafMask)
            afMask_tmp[iz] = aaafMask[iz][iy][ix];
          if (normalize)
            afDenom_tmp[iz] = aaafDenom[iz][iy][ix];
        }

        // apply the filter to the 1-D temporary array (afSource_tmp)
        aFilter[d].Apply(image_size[d],
                         afSource_tmp,
                         afDest_tmp,  //<-store filtered result here
                         afMask_tmp,
                         afDenom_tmp);//<-store sum of weights considered here

        // copy the results from the temporary filters back into the 3D arrays
        for (int iz = 0; iz < image_size[2]; iz++) {
          aaafDest[iz][iy][ix] = afDest_tmp[iz];
          if (normalize)
            aaafDenom[iz][iy][ix] = afDenom_tmp[iz]; //copy back into aaafDenom
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
    RealNum *afDest_tmp   = new RealNum [image_size[d]];
    RealNum *afSource_tmp = new RealNum [image_size[d]];
    //RealNum *afMask_tmp   = NULL;
    //if (aaafMask)
    //  afMask_tmp = new RealNum [image_size[d]];
    RealNum *afDenom_src_tmp = NULL;
    RealNum *afDenom_tmp = NULL;
    if (normalize) {
      afDenom_src_tmp = new RealNum [image_size[d]];
      afDenom_tmp     = new RealNum [image_size[d]];
    }

    #pragma omp for collapse(2)
    for (int iz = 0; iz < image_size[2]; iz++) {
      for (int ix = 0; ix < image_size[0]; ix++) {

        // copy the data we need to the temporary arrays
        for (int iy = 0; iy < image_size[1]; iy++) {
          afSource_tmp[iy] = aaafDest[iz][iy][ix];  //copy from prev aaafDest
          //if (aaafMask)
          //  afMask_tmp[iy] = aaafMask[iz][iy][ix];
          if (normalize)
            afDenom_src_tmp[iy] = aaafDenom[iz][iy][ix];
        }

        // apply the filter to the 1-D temporary array (afSource_tmp)
        aFilter[d].Apply(image_size[d],
                         afSource_tmp,
                         afDest_tmp); //<-store filtered result here
        if (normalize)
          // apply the filter to the 1-D array for the denominator as well
          aFilter[d].Apply(image_size[d],
                           afDenom_src_tmp, //<-weights so far (summed along z)
                           afDenom_tmp);//<-store sum of weights considered here

        // copy the results from the temporary filters back into the 3D arrays
        for (int iy = 0; iy < image_size[1]; iy++) {
          aaafDest[iz][iy][ix] = afDest_tmp[iy];
          if (normalize)
            aaafDenom[iz][iy][ix] = afDenom_tmp[iy]; //copy back into aaafDenom
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
    RealNum *afDest_tmp   = new RealNum [image_size[d]];
    RealNum *afSource_tmp = new RealNum [image_size[d]];
    //RealNum *afMask_tmp   = NULL;
    //if (aaafMask)
    //  afMask_tmp = new RealNum [image_size[d]];
    RealNum *afDenom_src_tmp = NULL;
    RealNum *afDenom_tmp = NULL;
    if (normalize) {
      afDenom_src_tmp = new RealNum [image_size[d]];
      afDenom_tmp     = new RealNum [image_size[d]];
    }

    #pragma omp for collapse(2)
    for (int iz = 0; iz < image_size[2]; iz++) {
      for (int iy = 0; iy < image_size[1]; iy++) {

        // copy the data we need to the temporary arrays
        for (int ix = 0; ix < image_size[0]; ix++) {
          afSource_tmp[ix] = aaafDest[iz][iy][ix];  //copy from prev aaafDest
          //if (aaafMask)
          //  afMask_tmp[ix] = aaafMask[iz][iy][ix];
          if (normalize)
            afDenom_src_tmp[ix] = aaafDenom[iz][iy][ix];
        }

        // apply the filter to the 1-D temporary array (afSource_tmp)
        aFilter[d].Apply(image_size[d],
                         afSource_tmp,
                         afDest_tmp); //<-store filtered result here
        if (normalize)
          // apply the filter to the 1-D array for the denominator as well
          aFilter[d].Apply(image_size[d],
                           afDenom_src_tmp, //<-weights so far (summed along z)
                           afDenom_tmp);//<-store sum of weights considered here

        // copy the results from the temporary filters back into the 3D arrays
        for (int ix = 0; ix < image_size[0]; ix++) {
          aaafDest[iz][iy][ix] = afDest_tmp[ix];
          if (normalize)
            aaafDenom[iz][iy][ix] = afDenom_tmp[ix]; //copy back into aaafDenom
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
  } //#pragma omp parallel private(afDest_tmp, ...)



  if (normalize) {
    assert(aaafDenom);
    for (int iz = 0; iz < image_size[2]; iz++)
      for (int iy = 0; iy < image_size[1]; iy++)
        for (int ix = 0; ix < image_size[0]; ix++)
          if (aaafDenom[iz][iy][ix] > 0.0)
            aaafDest[iz][iy][ix] /= aaafDenom[iz][iy][ix];
    // Cleanup
    Dealloc3D(image_size, &afDenom, &aaafDenom);
  }

  // Optional: filter(x) =  A  * exp(-(1/2)*((x/σ_x)^2 + (y/σ_y)^2 + (z/σ_z)^2))
  //                     = A_x * exp(-(1/2)*(x/σ_x)^2) *
  //                       A_y * exp(-(1/2)*(y/σ_y)^2) *
  //                       A_z * exp(-(1/2)*(z/σ_z)^2)
  // The 3D Gaussian coefficient ("A") is the product of the
  // 1-D Gaussian coefficients in the X,Y,Z directions (A_x * A_y * A_z).
  // Those coefficients happen to equal the value of the corresponding 1-D
  // Gaussian evaluated at its peak, which is stored in the central entry

  RealNum A_coeff = (aFilter[0].afH[0] *
                     aFilter[1].afH[0] *
                     aFilter[2].afH[0]);

  return A_coeff;

} //_ApplyGauss3D(aFilter)




/// @brief Apply a Gaussian filter (blur) to an image
///
/// @code h(x,y,z)=A*exp(-0.5*((x/σ_x)^2+(y/σ_y)^2+(z/σ_z)^2) @endcode
/// In this version, the user manually specifies the filter window width.
/// Voxels outside the boundary of the image or outside the mask are not
/// considered during the averaging (blurring) process, and the resulting
/// after blurring is weighted accordingly.  (normalize=true by default)
///
/// @returns the "A" coefficient (determined by normalization)

template<class RealNum>
RealNum
ApplyGauss3D(int const image_size[3], //!< image size in x,y,z directions
             RealNum const *const *const *aaafSource,   //!< source image (3D array)
             RealNum ***aaafDest,     //!< filtered (blurred) image stored here
             RealNum const *const *const *aaafMask,     //!< ignore voxels if aaafMask[i][j][k]==0
             RealNum const sigma[3],  //!< Gaussian sigma parameters σ_x,σ_y,σ_z
             int const truncate_halfwidth[3], //!< the filter window width
             bool normalize = true,           //!< normalize the average?
             ostream *pReportProgress = NULL  //!< print progress to the user?
             )
{
  assert(aaafSource);
  assert(aaafDest);
  //assert(aaafMask);
  //Allocate filters in all 3 directions.  (Later apply them sequentially.)
  Filter1D<RealNum, int> aFilter[3];
  for (int d=0; d < 3; d++)
    aFilter[d] = GenFilterGauss1D(sigma[d], truncate_halfwidth[d]);

  return _ApplyGauss3D(image_size, 
                       aaafSource,
                       aaafDest,
                       aaafMask,
                       aFilter,
                       normalize,
                       pReportProgress);

} //ApplyGauss3D(sigma, truncate_halfwidth, ...)




/// @brief Apply a Gaussian filter (blur) to a 3D image
///
/// @code h(x,y,z)=A*exp(-0.5*((x/σ_x)^2+(y/σ_y)^2+(z/σ_z)^2) @endcode
/// The constant "A" is determined by normalization.
/// In this version, the user specifies the filter window-width in units of σ
/// Voxels outside the boundary of the image or outside the mask are not
/// considered during the averaging (blurring) process, and the resulting
/// after blurring is weighted accordingly.  (normalize=true by default)
///
/// @returns the "A" coefficient (determined by normalization)

template<class RealNum>
RealNum
ApplyGauss3D(int const image_size[3], //!< image size in x,y,z directions
             RealNum const *const *const *aaafSource,   //!< source image (3D array)
             RealNum ***aaafDest,     //!< filtered (blurred) image stored here
             RealNum const *const *const *aaafMask,     //!< ignore voxels if aaafMask[i][j][k]==0
             RealNum const sigma[3],  //!< Gaussian sigma parameters σ_x,σ_y,σ_z
             RealNum truncate_ratio = 2.5,  //!< how many sigma before truncating?
             bool normalize = true,           //!< normalize the average?
             ostream *pReportProgress = NULL  //!< print progress to the user?
             )
{
  RealNum A;
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

template<class RealNum>
void
ApplyDog3D(int const image_size[3], //!< image size in x,y,z directions
           RealNum const *const *const *aaafSource,   //!< source image (3D array)
           RealNum ***aaafDest,     //!< filtered image stored here
           RealNum const *const *const *aaafMask,     //!< ignore voxels if aaafMask[i][j][k]==0
           RealNum const sigma_a[3], //!< (a_x, a_y, a_z) in the formula above
           RealNum const sigma_b[3], //!< (b_x, b_y, b_z) in the formula above
           int const truncate_halfwidth[3],//!<(half the filter window width along x,y,z)
           RealNum *pA = NULL, //!< Optional: report "A" (normalized coefficient) to the caller?
           RealNum *pB = NULL, //!< Optional: report "B" (normalized coefficient) to the caller?
           ostream *pReportProgress = NULL  //!< print progress to the user?
           )
{
  RealNum ***aaafTemp; //temporary array to store partially processed tomogram
  RealNum *afTemp;     //temporary array to store partially processed tomogram

  if (pReportProgress)
    *pReportProgress
      << " -- Attempting to allocate space for one more image. --\n"
      << " -- (If this crashes your computer, find a computer  --\n"
      << " --  with more RAM and use \"ulimit\")                 --\n";

  Alloc3D(image_size,
          &afTemp,
          &aaafTemp);

  RealNum A, B;        // let the user know what A B coefficients were used

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

template<class RealNum>
void
ApplyDogScaleFree3D(int const image_size[3], //!< source image size
                    RealNum const *const *const *aaafSource,   //!< source image (3D array)
                    RealNum ***aaafDest,     //!< filtered image stored here
                    RealNum const *const *const *aaafMask,     //!< ignore voxels if aaafMask[i][j][k]==0
                    RealNum const sigma[3],  //!< Gaussian width in x,y,z drections
                    RealNum delta_sigma_over_sigma, //difference in Gauss widths
                    RealNum truncate_ratio,  //!< how many sigma before truncating?
                    RealNum *pA = NULL, //!< Optional: report "A" (normalized coefficient) to the caller?
                    RealNum *pB = NULL, //!< Optional: report "B" (normalized coefficient) to the caller?
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
  RealNum sigma_a[3];
  RealNum sigma_b[3];
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

  RealNum t_over_delta_t = 1.0 / SQR(delta_sigma_over_sigma);

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

template<class RealNum>
void
ApplyDogScaleFree3D(int const image_size[3], //!< source image size
                    RealNum const *const *const *aaafSource, //!< source image
                    RealNum ***aaafDest,     //!< save results here
                    RealNum const *const *const *aaafMask,  //!< ignore voxels where mask==0
                    RealNum sigma,  //!< Gaussian width in x,y,z drections
                    RealNum delta_sigma_over_sigma, //!< δ, difference in Gauss widths
                    RealNum truncate_ratio=2.5,  //!< how many sigma before truncating?
                    RealNum *pA = NULL, //!< Optional: report "A" (normalized coefficient) to the caller?
                    RealNum *pB = NULL, //!< Optional: report "B" (normalized coefficient) to the caller?
                    ostream *pReportProgress = NULL  //!< print progress to the user?
                    )
{
  RealNum sigma_xyz[3] = {sigma, sigma, sigma};
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

template<class RealNum>
void
BlobDog(int const image_size[3], //!< source image size
        RealNum const *const *const *aaafSource,   //!< source image
        RealNum const *const *const *aaafMask,     //!< ignore voxels where mask==0
        const vector<RealNum>& blob_sigma, //!< blob widths to try, ordered
        vector<array<RealNum,3> >& minima_crds, //!< store minima x,y,z coords here
        vector<array<RealNum,3> >& maxima_crds, //!< store maxima x,y,z coords here
        vector<RealNum>& minima_sigma, //!< corresponding width for that minima
        vector<RealNum>& maxima_sigma, //!< corresponding width for that maxima
        vector<RealNum>& minima_scores, //!< what was the blob's score?
        vector<RealNum>& maxima_scores, //!< (score = intensity after filtering)
        RealNum delta_sigma_over_sigma=0.02,//!< δ param for approximating LOG with DOG
        RealNum truncate_ratio=2.8,      //!< how many sigma before truncating?
        RealNum minima_threshold=0.0,    //!< discard blobs with unremarkable scores
        RealNum maxima_threshold=0.0,    //!< discard blobs with unremarkable scores
        bool use_threshold_ratios=true, //!< threshold=ratio*best_score ?
        // optional arguments
        ostream *pReportProgress = NULL, //!< optional: report progress to the user?
        RealNum ****aaaafI = NULL, //!<optional: preallocated memory for filtered images (indexable)
        RealNum **aafI = NULL      //!<optional: preallocated memory for filtered images (contiguous)
        )

{

  // We need 3 images to store the result of filtering the image
  // using DOG filters with 3 different widths.  Store those images here:

  if (pReportProgress)
    *pReportProgress
      << " -- Attempting to allocate space for 3 more images.  --\n"
      << " -- (If this crashes your computer, find a computer  --\n"
      << " --  with more RAM and use \"ulimit\")               --\n";

  bool preallocated = ! (aaaafI == NULL);

  if (! preallocated) {
    assert(aaaafI == NULL);
    assert(aafI == NULL);
    aaaafI = new RealNum*** [3];
    aafI = new RealNum* [3];
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


  RealNum global_min_score = 1.0;  //impossible, minima < 0
  RealNum global_max_score = -1.0; //impossible, maxima > 0

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

      vector<array<RealNum,3> > min_crds_proc; //store minima x,y,z coords here
      vector<array<RealNum,3> > max_crds_proc; //store maxima x,y,z coords here
      vector<RealNum> min_sigma_proc; //corresponding width for that minima
      vector<RealNum> max_sigma_proc; //corresponding width for that maxima
      vector<RealNum> min_scores_proc; //what was the blob's score?
      vector<RealNum> max_scores_proc; //(score = intensity after filtering)
      RealNum global_min_score_proc = global_min_score;
      RealNum global_max_score_proc = global_max_score;

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
                    RealNum entry    = aaaafI[ j2i[0]  ][ iz ][ iy ][ ix ];
                    RealNum neighbor = aaaafI[ j2i[jr] ][ Iz ][ Iy ][ Ix ];
                    if (neighbor <= entry)
                      is_minima = false;
                    if (neighbor >= entry)
                      is_maxima = false;
                  }
                }
              }
            }

            RealNum score = aaaafI[ j2i[0] ][ iz ][ iy ][ ix ];

            if ((! aaafMask) || (aaafMask[iz][iy][ix] != 0))
            {
              RealNum minima_threshold_so_far = minima_threshold;
              if (use_threshold_ratios)
                minima_threshold_so_far=minima_threshold*global_min_score_proc;
              if (is_minima &&
                  (score < 0.0) &&
                  ((score < minima_threshold_so_far) // || disable_thresholds
                   ))   
              {
                array<RealNum, 3> ixiyiz;
                ixiyiz[0] = ix;
                ixiyiz[1] = iy;
                ixiyiz[2] = iz;
                min_crds_proc.push_back(ixiyiz);
                min_sigma_proc.push_back(blob_sigma[ir-1]);
                min_scores_proc.push_back(score);
                if (score < global_min_score_proc)
                  global_min_score_proc = score;
              }

              RealNum maxima_threshold_so_far = maxima_threshold;
              if (use_threshold_ratios)
                maxima_threshold_so_far=maxima_threshold*global_max_score_proc;
              if (is_maxima &&
                  (score > 0.0) &&
                  ((score > maxima_threshold_so_far) // || disable_thresholds
                   ))   
              {
                array<RealNum, 3> ixiyiz;
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

template<class RealNum>
void
BlobDogD(int const image_size[3], //!<source image size
         RealNum const *const *const *aaafSource,   //!< source image
         RealNum const *const *const *aaafMask,     //!< ignore voxels where mask==0
         const vector<RealNum>& blob_diameters, //!<blob widths to try, ordered
         vector<array<RealNum,3> >& minima_crds, //!<store minima x,y,z coords here
         vector<array<RealNum,3> >& maxima_crds, //!<store maxima x,y,z coords here
         vector<RealNum>& minima_diameters, //!<corresponding width for that minima
         vector<RealNum>& maxima_diameters, //!<corresponding width for that maxima
         vector<RealNum>& minima_scores, //!<what was the blob's score?
         vector<RealNum>& maxima_scores, //!<(score = intensity after filtering)
         // optional arguments
         RealNum delta_sigma_over_sigma=0.02,//!<param for approximating LOG with DOG
         RealNum truncate_ratio=2.5,    //!<how many sigma before truncating?
         RealNum minima_threshold=0.5,  //!<discard blobs with unremarkable scores
         RealNum maxima_threshold=0.5,  //!<discard blobs with unremarkable scores
         bool    use_threshold_ratios=true, //!<threshold=ratio*best_score?
         ostream *pReportProgress = NULL, //!<report progress to the user?
         RealNum ****aaaafI = NULL, //!<preallocated memory for filtered images
         RealNum **aafI = NULL     //!<preallocated memory for filtered images (conserve memory)
         )
{

  vector<RealNum> minima_sigma;
  vector<RealNum> maxima_sigma;
  vector<RealNum> blob_sigma(blob_diameters.size());
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

template<class RealNum>
void
VisualizeBlobs(int const image_size[3], //!< image size
               RealNum ***aaafImage,  //!< array where we should write new image
               RealNum const *const *const *aaafMask,   //!< ignore voxels where mask==0
               vector<array<RealNum,3> > &centers, //!< center of each blob
               vector<RealNum> &diameters,         //!< diameter of eachs pherical shell (in voxels)
               vector<RealNum> &shell_thicknesses, //!< how thick is each spherical shell (in voxels)
               vector<RealNum> &voxel_intensities_foreground, //!< voxels in spherical shell get this value
               RealNum voxel_intensity_background = 0.0, //!< assign background voxels to this value
               RealNum voxel_intensity_background_rescale = 0.25, //!< superimpose with old image? Which weight?
               bool voxel_intensity_foreground_normalize = false, //!< divide brightnesses by number of voxels in spherical shell? (rarely useful)
               ostream *pReportProgress = NULL //!<optional: report progress to the user?
               )
{

  RealNum tomo_ave  =  AverageArr(image_size,
                                  aaafImage,
                                  aaafMask);
  RealNum tomo_stddev  =  StdDevArr(image_size,
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
    RealNum Rssqr_max = SQR(diameters[i]/2);
    RealNum Rssqr_min = 0.0;
    if ((shell_thicknesses[i] > 0.0) && (diameters[i]/2 - shell_thicknesses[i] > 0.0))
      Rssqr_min = SQR(diameters[i]/2 - shell_thicknesses[i]);

    if (pReportProgress)
      *pReportProgress << ", diameter=" << diameters[i] << ", th=" << shell_thicknesses[i]
                       <<"\n";

    // Normalize the brightness of each sphere?
    // (ie by dividing the intensity by the number of voxels in the sphere)
    RealNum imultiplier = 1.0;
    long nvoxelspersphere = 1;
    if (voxel_intensity_foreground_normalize) {
      nvoxelspersphere = 0;
      for (int jz = -Rs; jz <= Rs; jz++) {
        for (int jy = -Rs; jy <= Rs; jy++) {
          for (int jx = -Rs; jx <= Rs; jx++) {
            RealNum rsqr = jx*jx + jy*jy + jz*jz;
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
template<class RealNum>
void
SortBlobs(vector<array<RealNum,3> >& blob_crds, //!< x,y,z of each blob's center
          vector<RealNum>& blob_diameters,  //!< the width of each blob
          vector<RealNum>& blob_scores,  //!< the score for each blob
          bool descending_order = true, //!<sort scores in ascending or descending order?
          ostream *pReportProgress = NULL //!< optional: report progress to the user?
          )
{ 
  size_t n_blobs = blob_crds.size();
  assert(n_blobs == blob_diameters.size());
  assert(n_blobs == blob_scores.size());
  vector<tuple<RealNum, size_t> > score_index(n_blobs);
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

template<class RealNum>
RealNum CalcVolOverlap(RealNum rij,//!<the distance between the spheres' centers
                       RealNum Ri, //!< the radius of sphere i
                       RealNum Rj  //!< the radius of sphere j
                       )
{
  // WLOG, assume Ri <= Rj.  Insure that below
  if (Ri > Rj) {
    RealNum tmp = Ri;
    Ri = Rj;
    Rj = tmp;
  }
  if (rij <= Ri) {
    return (4*M_PI/3) * Ri*Ri*Ri;
  }

  // "xi" and "xj" are the distances from the sphere centers
  // to the plane where the two spheres intersect.
  RealNum xi = 0.5 * (1.0/rij) * (rij*rij + Ri*Ri - Rj*Rj);
  RealNum xj = 0.5 * (1.0/rij) * (rij*rij + Rj*Rj - Ri*Ri);
  RealNum volume_overlap =
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
template<class RealNum>
void
DiscardOverlappingBlobs(vector<array<RealNum,3> >& blob_crds,
                        vector<RealNum>& blob_diameters, 
                        vector<RealNum>& blob_scores,
                        SortBlobCriteria sort_blob_criteria, //!< priority to high or low scoring blobs?
                        RealNum min_radial_separation_ratio, //!< discard blobs too close
                        RealNum max_volume_overlap_large, //!< discard blobs which overlap too much with the large blob
                        RealNum max_volume_overlap_small, //!< discard blobs which overlap too much with the small blob
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
      << "     -- Attempting to allocate space for one more image map --\n"
      << "     -- (If this crashes your computer, find a computer     --\n"
      << "     --  with more RAM and use \"ulimit\")                    --\n";

  // Occupancy table
  //     (originally named "bool ***aaabOcc")
  // This table stores a list of blobs which occupy a given voxel.

  int occupancy_table_size[3];
  int bounds_min[3] = {0, 0, 0};
  int bounds_max[3] = {-1, -1, -1};
  
  for (int i=0; i < blob_crds.size(); i++) {
    for (int d=0; d < 3; d++) {
      RealNum reff = ceil(blob_diameters[i]/2); //blob radii in units of voxels
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
    RealNum reff_ = blob_diameters[i]/2; //blob radii in units of voxels
    RealNum Reff_ = reff_ / scale; //blob radii expressed in "low rez" units
    int Reff = ceil(Reff_);         //round up
    int Reffsq = ceil(Reff_*Reff_);
    RealNum ix = blob_crds[i][0];      //coordinates of the center of the blob
    RealNum iy = blob_crds[i][1];
    RealNum iz = blob_crds[i][2];
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
            RealNum kx = blob_crds[k][0];
            RealNum ky = blob_crds[k][1];
            RealNum kz = blob_crds[k][2];
            RealNum rik = sqrt((ix-kx)*(ix-kx)+(iy-ky)*(iy-ky)+(iz-kz)*(iz-kz));
            RealNum vol_overlap = CalcVolOverlap(rik,
                                                 blob_diameters[i]/2,
                                                 blob_diameters[k]/2);
            RealNum ri = blob_diameters[i]/2;
            RealNum rk = blob_diameters[k]/2;
            if (rik < (ri + rk) * min_radial_separation_ratio)
              discard = true;
            RealNum vi = (4*M_PI/3)*(ri*ri*ri);
            RealNum vk = (4*M_PI/3)*(rk*rk*rk);
            RealNum v_large = vi;
            RealNum v_small = vk;
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






#endif //#ifndef _FILTER3D_H
