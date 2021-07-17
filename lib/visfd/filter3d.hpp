///   @file filter3d.hpp
///   @brief classes and functions that apply filters to 3D arrays
///   @author Andrew Jewett
///   @date 2018-2-26

#ifndef _FILTER3D_HPP
#define _FILTER3D_HPP

#include <cstring>
#include <ostream>
#include <cassert>
#include <cmath>
#include <limits>
#include <algorithm>
#ifndef CXX17_UNSUPPORTED
  #include <execution>
#endif
using namespace std;
#include <alloc3d.hpp>    // defines Alloc3D() and Dealloc3D()
#include <filter1d.hpp>   // defines "Filter1D" (used in ApplySeparable())
#include <err_visfd.hpp>  // defines the "VisfdErr" exception type




namespace visfd {


/// @class Filter3D
/// @brief A class for general linear (convolutional) filters in 3D
///
/// @note  In practice, this class is not used often because most operations
///        we need can be performed with separable filters which are faster.

template<typename Scalar, typename Integer>

class Filter3D {

public:
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
  /// where: f[i] is the original image at position i (ix,iy,iz)
  ///      i = shorthand for (ix,iy,iz) = a location in the filtered image
  ///      j = shorthand for (jx,jy,jz) is summed over the entries in the filter
  ///    i-j = shorthand for (ix-jx, iy-jy, iz-jz)
  ///        g[i] is the image after the filter has been applied
  ///        h[j] is the filter
  ///        mask[i] selects the pixels we care about (usually either 0 or 1)
  ///          (If not specified, assume it is 1 inside the image, 0 outside.)
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
             Scalar const *const *const *aaafMask = nullptr,
             bool normalize = false,
             ostream *pReportProgress = nullptr) const
  {
    Scalar *afDenominator = nullptr;
    Scalar ***aaafDenominator = nullptr;
    if (normalize)
      aaafDenominator = Alloc3D<Scalar>(size_source);

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
      Dealloc3D(aaafDenominator);
      aaafDenominator = nullptr;
    }
  }

  /// @brief  Apply the filter to a 3D image (aaafSource[][][]).
  ///         This version is identical to the other version of Apply()
  ///         except that this version returns both d(i) and g(i) whenever
  ///         you supply a non-null afDenominator[] argument (see below).
  ///         It also does NOT normalize the result (by dividing g(i) / d(i)).
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
  /// If aaafDenominator!=nullptr, then it also computes the "denominator":
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
  ///        d[i] is the "denominator" = sum of all the filter weights h[j]...
  ///             ...considered when calculating g[i] (weighted by the mask)
  ///        mask[i] selects the pixels we care about (usually either 0 or 1)
  ///          (If not specified, assume it is 1 inside the image, 0 outside.)
  /// @endcode
  ///
  /// @param size_source contains size of the source image (in the x,y,z directions)
  /// @param aaafSource[][][] is the source array (source image) <==> "f[i]"
  /// @param aaafDest[][][] will store the image after filtering <==> "g[i]"
  /// @param aaafMask[][][]==0 whenever we want to ignore entries in afSource[][]. Optional.
  /// @param aaafDenominator[][][] will store d[i] if you supply a non-null pointer

  void Apply(Integer const size_source[3],
             Scalar const *const *const *aaafSource,
             Scalar ***aaafDest,
             Scalar const *const *const *aaafMask = nullptr,
             Scalar ***aaafDenominator = nullptr,
             ostream *pReportProgress = nullptr) const
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
            if (aaafDenominator[iz][iy][ix])
              aaafDenominator[iz][iy][ix] = 0.0;
            continue;
          }

          aaafDest[iz][iy][ix] =
            ApplyToVoxel(ix, iy, iz,
                         size_source,
                         aaafSource,
                         aaafMask,
                         (aaafDenominator
                          ? &(aaafDenominator[iz][iy][ix])
                          : nullptr));
        }
      }
    }
  } // Apply()


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


  Filter3D(Integer const set_halfwidth[3]) {
    Init();
    Resize(set_halfwidth);
  }


  Filter3D() {
    Init();
  }


  // destructor, copy and move constructor, swap, and assignment operator

  ~Filter3D() {
    Dealloc();
  }


  Filter3D(const Filter3D<Scalar, Integer>& source) {
    Init();
    Resize(source.halfwidth); // allocates and initializes afH and aaafH
    for(Integer iz=-halfwidth[2]; iz<=halfwidth[2]; iz++)
      for(Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
        for(Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
          aaafH[iz][iy][ix] = source.aaafH[iz][iy][ix];
  }


  void swap(Filter3D<Scalar, Integer> &other) {
    std::swap(aaafH, other.aaafH);
    std::swap(halfwidth, other.halfwidth);
    std::swap(array_size, other.array_size);
  }

  // Move constructor (C++11)
  Filter3D(Filter3D<Scalar, Integer>&& other) {
    Init();
    this->swap(other);
  }

  // Using the "copy-swap" idiom for the assignment operator
  Filter3D<Scalar, Integer>&
    operator = (Filter3D<Scalar, Integer> source) {
    this->swap(source);
    return *this;
  }



  // Miscellaneous functions
  //    (COMMENT: Most of these functions are no longer needed.
  //              I may remove them in the future. -A 2021-7-11)
  

  /// @brief Calculate the (weighted) average value of the filter array aaafH
  /// @param aaafW optional weights used when calculating the averages
  /// @return the (weighted) average value of the filter array aaafH
  Scalar Average(Scalar const *const *const *aaafW=nullptr) const {
    Scalar numer = 0.0;
    Scalar denom = 0.0;
    for (int iz = -halfwidth[2]; iz <= halfwidth[2]; iz++) {
      for (int iy = -halfwidth[1]; iy <= halfwidth[1]; iy++) {
        for (int ix = -halfwidth[0]; ix <= halfwidth[0]; ix++) {
          Scalar w = 1.0;
          if (aaafW)
            w = aaafW[iz][iy][ix];
          numer += w * aaafH[iz][iy][ix];
          denom += w;
        }
      }
    }
    return numer / denom;
  }

  /// @brief Calculate the (weighted) average squared values in the filter array, aaafH
  /// @param aaafW optional weights used when calculating the averages
  /// @return the (weighted) average squared values in the filter array aaafH
  Scalar AverageSqr(Scalar const *const *const *aaafW=nullptr) const {
    Scalar numer = 0.0;
    Scalar denom = 0.0;
    for (int iz = -halfwidth[2]; iz <= halfwidth[2]; iz++) {
      for (int iy = -halfwidth[1]; iy <= halfwidth[1]; iy++) {
        for (int ix = -halfwidth[0]; ix <= halfwidth[0]; ix++) {
          Scalar w = 1.0;
          if (aaafW)
            w = aaafW[iz][iy][ix];
          numer += w * SQR(aaafH[iz][iy][ix]);
          denom += w;
        }
      }
    }
    return numer / denom;
  }

  /// @brief Calculate the (weighted) standard deviation of the filter array values
  /// @param aaafW optional weights used when calculating the standard deviation
  /// @return the (weighted) standard deviation of the filter array aaafH
  Scalar StdDev(Scalar const *const *const *aaafW=nullptr) const {
    Scalar ave = Average(aaafW);
    Scalar numer = 0.0;
    Scalar denom = 0.0;
    for (int iz = -halfwidth[2]; iz <= halfwidth[2]; iz++) {
      for (int iy = -halfwidth[1]; iy <= halfwidth[1]; iy++) {
        for (int ix = -halfwidth[0]; ix <= halfwidth[0]; ix++) {
          Scalar w = 1.0;
          if (aaafW)
            w = aaafW[iz][iy][ix];
          numer += w * SQR(aaafH[iz][iy][ix] - ave);
          denom += w;
        }
      }
    }
    return static_cast<Scalar>(sqrt(numer / denom));
  }

  /// @brief Calculate the (weighted) sum of the filter array values, aaafH
  /// @param aaafW optional weights used when calculating the sum
  /// @return the (weighted) sum of the filter array values
  Scalar Sum(Scalar const *const *const *aaafW=nullptr) const {
    Scalar sum = 0.0;
    for (int iz = -halfwidth[2]; iz <= halfwidth[2]; iz++) {
      for (int iy = -halfwidth[1]; iy <= halfwidth[1]; iy++) {
        for (int ix = -halfwidth[0]; ix <= halfwidth[0]; ix++) {
          Scalar w = 1.0;
          if (aaafW)
            w = aaafW[iz][iy][ix];
          sum += w * aaafH[iz][iy][ix];
        }
      }
    }
    return sum;
  }

  /// @brief Calculate the (weighted) sum of the squared filter array values
  /// @param aaafW optional weights used when calculating the sum
  /// @return the (weighted) sum of the squared filter array values
  Scalar SumSqr(Scalar const *const *const *aaafW=nullptr) const {
    Scalar sum = 0.0;
    for (int iz = -halfwidth[2]; iz <= halfwidth[2]; iz++) {
      for (int iy = -halfwidth[1]; iy <= halfwidth[1]; iy++) {
        for (int ix = -halfwidth[0]; ix <= halfwidth[0]; ix++) {
          Scalar w = 1.0;
          if (aaafW)
            w = aaafW[iz][iy][ix];
          sum += w * SQR(aaafH[iz][iy][ix]);
        }
      }
    }
    return sum;
  }

  /// @brief Add a number to all of the filter array values, aaafH
  /// @param offset  the number to add
  void AddScalar(Scalar offset) {
    for (int iz = -halfwidth[2]; iz <= halfwidth[2]; iz++)
      for (int iy = -halfwidth[1]; iy <= halfwidth[1]; iy++)
        for (int ix = -halfwidth[0]; ix <= halfwidth[0]; ix++)
          aaafH[iz][iy][ix] += offset;
  }

  /// @brief multiply all of the filter array values (aaafH) by a number
  /// @param offset  the number to multiply
  void MultiplyScalar(Scalar scale) {
    for (int iz = -halfwidth[2]; iz <= halfwidth[2]; iz++)
      for (int iy = -halfwidth[1]; iy <= halfwidth[1]; iy++)
        for (int ix = -halfwidth[0]; ix <= halfwidth[0]; ix++)
          aaafH[iz][iy][ix] *= scale;
  }



 private:

  /// @brief Apply a filter to the source image (aaafSource) at a particular
  ///        voxel location.
  /// @param ix  the voxel's position in that 3D image
  /// @param iy  the voxel's position in that 3D image
  /// @param iz  the voxel's position in that 3D image
  /// @param size_source contains size of the source image (in the x,y,z directions)
  /// @param aaafSource[][][] is the source image
  /// @param optional. aaafMask[i][j][k]==0 for voxels we want to exclude from consideration
  /// @param pDenominator=if you want to store the sum of the weights considered, pass a pointer to a number
  /// (useful if the sum was not complete due to some voxels being masked out,
  ///  or because the filter extends beyond the boundaries of the image)
  /// @note: This is a private function, not intended for public use.
  
  Scalar ApplyToVoxel(Integer ix,
                      Integer iy,
                      Integer iz,
                      Integer const size_source[3],
                      Scalar const *const *const *aaafSource,
                      Scalar const *const *const *aaafMask = nullptr,
                      Scalar *pDenominator = nullptr) const

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
            if (mask_val == 0.0)
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
    aaafH = Alloc3D<Scalar>(array_size);
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
      Init();
      return;
    }
    // Integer array_size[3];
    for(int d=0; d < 3; d++) {
      array_size[d] = 1 + 2*halfwidth[d];
    }
    // shift pointers back to normal
    for (Integer iz = -halfwidth[2]; iz <= halfwidth[2]; iz++) {
      for (Integer iy = -halfwidth[1]; iy <= halfwidth[1]; iy++) {
        aaafH[iz][iy] -= halfwidth[0];
      }
      aaafH[iz] -= halfwidth[1];
    }
    aaafH -= halfwidth[2];
    // then deallocate
    Dealloc3D(aaafH);
    aaafH = nullptr;
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
    aaafH = nullptr;
  }
  
}; // class Filter3D






/// @brief Create a 3D filter and fill it with a "generalized Gaussian" function
///   @verbatim h(x,y,z) = A*exp(-r^m)  @endverbatim
/// where  @verbatim r  = sqrt((x/σ_x)^2 + (y/σ_y)^2 + (z/σ_z)^2) @endverbatim
///   and  "A" is determined by normalization of the discrete sum
/// @note "A" is equal to the value stored in the middle of the array,
///       The caller can determine what "A" is by looking at this value.

template<typename Scalar>

Filter3D<Scalar, int>
GenFilterGenGauss3D(const Scalar width[3], //!< "σ_x", "σ_y", "σ_z" parameters
                    Scalar m_exp,          //!< "m" exponent parameter
                    const int truncate_halfwidth[3], //!< size of filter window
                    Scalar *pA=nullptr     //!< optional:report A coeff to user
                    )
{
  Scalar truncate_threshold = 1.0;
  for (int d=0; d<3; d++) {
    Scalar h = ((width[d]>0)
                 ? exp(-pow(truncate_halfwidth[d]/width[d], m_exp))
                 : 1.0);
    if (h < truncate_threshold)
      truncate_threshold = h;
  }
  Filter3D<Scalar, int> filter(truncate_halfwidth);
  Scalar total = 0;
  for (int iz=-filter.halfwidth[2]; iz<=filter.halfwidth[2]; iz++) {
    for (int iy=-filter.halfwidth[1]; iy<=filter.halfwidth[1]; iy++) {
      for (int ix=-filter.halfwidth[0]; ix<=filter.halfwidth[0]; ix++) {
        //Scalar x = ix/width[0];
        //Scalar y = iy/width[1];
        //Scalar z = iz/width[2];
        Scalar x = ((! ((width[0] == 0.0) && (ix==0))) ? ix/width[0] : 0.0);
        Scalar y = ((! ((width[1] == 0.0) && (iy==0))) ? iy/width[1] : 0.0);
        Scalar z = ((! ((width[2] == 0.0) && (iz==0))) ? iz/width[2] : 0.0);
        Scalar r = sqrt(x*x + y*y + z*z);
        Scalar h = ((r>0) ? exp(-pow(r, m_exp)) : 1.0);
        if (std::abs(h) < truncate_threshold)
          h = 0.0; //This eliminates corner entries which fall below threshold
                   //(and eliminates anisotropic artifacts due to these corners)
                   //There's no reason to keep any entries less than min value.
        filter.aaafH[iz][iy][ix] = h;
                    
        total += h;
      }
    }
  }

  // normalize:
  for (int iz=-filter.halfwidth[2]; iz<=filter.halfwidth[2]; iz++)
    for (int iy=-filter.halfwidth[1]; iy<=filter.halfwidth[1]; iy++)
      for (int ix=-filter.halfwidth[0]; ix<=filter.halfwidth[0]; ix++)
        filter.aaafH[iz][iy][ix] /= total;

  // The coefficient in front of the Gaussian ("A")
  // equals the height of its central peak,
  // which is located in the middle of the array
  Scalar A = filter.aaafH[0][0][0];

  if (pA) {
    *pA = A;
  }

  return filter;
} //GenFilterGenGauss3D(width, m_exp, truncate_halfwidth)







/// @brief Create a 3D filter and fill it with a "generalized Gaussian" function
///   @verbatim h(x,y,z) = A*exp(-r^m)  @endverbatim
/// where  @verbatim r  = sqrt((x/σ_x)^2 + (y/σ_y)^2 + (z/σ_z)^2) @endverbatim
///   and  "A" is determined by normalization of the discrete sum
/// @note "A" is equal to the value stored in the middle of the array (aaafH),
///       The caller can determine what "A" is by looking at aaafH there
/// @note The width of the filter window in each direction equals the "width"
///       in that direction multiplied by "filter_cutoff_ratio"

template<typename Scalar>

Filter3D<Scalar, int>
GenFilterGenGauss3D(const Scalar width[3], //!< "σ_x", "σ_y", "σ_z" parameters
                    Scalar m_exp,          //!< "m" parameter in formula
                    const Scalar filter_cutoff_ratio=2.5, //!< how many sigma (σ) before truncating?
                    Scalar *pA=nullptr     //!< optional:report A coeff to user
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
                             pA);
} //GenFilterGenGauss3D(width, m_exp, filter_cutoff_ratio)











///@brief ApplySeparable() applies a separable filter on a 3D array: aaafSource.
///        It assumes separate Filter1D objects have already been created 
///        which will blur the image successively in each direction (x,y,z).
///        The caller can specify an optional mask image (aaafMask), 
///        which allows us to exclude certain voxels from consideration.
///        This important function is invoked by ApplyGauss(), which
///        itself is invoked by almost every function in this library.
///        (However this function can be used for other non-Gaussian filters.)
///
/// @note: This is a low level function.  Most users should ignore it.
///
/// @return  The resulting blurred image is stored in the aaafDest[][][] array.
///          This function returns the effective height of the central peak of
///          the 3D filter.  (This is just the product of the central peaks of
///          each of the 1D filters.)  For normalized Gaussian shaped filters,
///          this peak height can be used to calculate the effective width
///          of the Gaussian to the caller (in case that information is useful).
///
/// @note  One unusual feature of this function is its ability to efficiently
///        normalize the resulting filtered image in the presence of a mask (as
///        well as near the image boundaries). Consequently, the image does not
///        necessarily fade to black near the boundaries of the image (or the
///        mask) but instead fades to the shade of the remaining nearby voxels.
///        This feature is enabled whenever the "normalize" argument is "true",
///        however this will slow the calculation. (If the aaafMask array is not
///        nullptr, then it will slow the calculation by a factor of up to 1.83)



// (Private comment:)
//     (Alternatively, this normalization could be handled by applying
//      the filter (blur) to both the mask, and the masked original image
//      and then dividing the two blurred images by each other.  But the
//      implementation here is a little bit faster than that: ~17% faster.)


template<typename Scalar>

Scalar
ApplySeparable(int const image_size[3],              //!<number of voxels in x,y,z directions
               Scalar const *const *const *aaafSource, //!<image to which we want to apply the filter
               Scalar ***aaafDest,                   //!<store the filtered image here
               Scalar const *const *const *aaafMask, //!<if not nullptr, ignore voxels if aaafMask[iz][iy][ix]!=0
               Filter1D<Scalar, int> aFilter[3],     //!<preallocated 1D filters
               bool normalize = true, //!< normalize the result near the boundaries?
               ostream *pReportProgress = nullptr) //!< print out progress to the user?
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
  Scalar ***aaafDenom = nullptr;
  Scalar *afDenom = nullptr;

  if (normalize) {
    if (aaafMask) {
      aaafDenom = Alloc3D<Scalar>(image_size);
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
    Scalar *afMask_tmp   = nullptr;
    if (aaafMask)
      afMask_tmp = new Scalar [image_size[d]];
    Scalar *afDenom_tmp = nullptr;
    if (normalize && aaafMask)
      afDenom_tmp = new Scalar [image_size[d]];

    #pragma omp for collapse(2)
    for (int iy = 0; iy < image_size[1]; iy++) {
      for (int ix = 0; ix < image_size[0]; ix++) {


        // --------------------------------------------
        // Heisenbug ?
        //     A "std::bad_alloc()" exception was thrown here using OpenMP 
        //     on 2018-10-18.  The problem was not apparent before this time,
        //     and does not occur when compiled without OpenMP support
        //     (More specifically, it does not occur when using the
        //      "setup_gcc_linux_dbg.sh" compile options which lack -fopenmp.)
        //
        //     This mysteriously seems to have "fixed" itself on 2018-10-24.
        //     I'm leaving a record of the problem here in case the problem
        //     resumes later on.
        //     (Incidentally, the location of the bug was determined 
        //      by filling the code with "fprintf()" statements.  
        //      Placing an "fprintf()" after this line had no effect.)
        // --------------------------------------------


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
            // (Note: if aaafMask==nullptr then we normalize using a faster method)
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
    //Scalar *afMask_tmp   = nullptr;
    //if (aaafMask)
    //  afMask_tmp = new Scalar [image_size[d]];
    Scalar *afDenom_src_tmp = nullptr;
    Scalar *afDenom_tmp = nullptr;
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
        // aaafDest[][][] and also afSource_tmp[].  Now we want to
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
            // (Note: if aaafMask==nullptr then we normalize using a faster method)
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
    //Scalar *afMask_tmp   = nullptr;
    //if (aaafMask)
    //  afMask_tmp = new Scalar [image_size[d]];
    Scalar *afDenom_src_tmp = nullptr;
    Scalar *afDenom_tmp = nullptr;
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
            // (Note: if aaafMask==nullptr then we normalize using a faster method)
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
      Dealloc3D(aaafDenom);
      aaafDenom = nullptr;
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
  // It might be useful or convenient to return this information to the caller.

  Scalar A_coeff = (aFilter[0].afH[0] *
                    aFilter[1].afH[0] *
                    aFilter[2].afH[0]);

  return A_coeff;

} //ApplySeparable(aFilter)




/// @brief Apply an ellipsoidal Gaussian filter (blur) to a 3D image
/// @code h(x,y,z)=A*exp(-0.5*((x/σ_x)^2+(y/σ_y)^2+(z/σ_z)^2) @endcode
/// In this version, the caller can specify all 3 components of σ (σ_x,σ_y,σ_z).
/// The filter must eventually be truncated.  In this version the caller
/// specifies the filter truncation window width (truncation_halfwidth)
/// in units of voxels.
///
/// If the user specifies a mask image (if aaafMask != nullptr), then voxels whose
/// corresponding entry in the aaafMask[][][] array equals 0 are ignored,
/// If the entry in the aaafMask[][][] is non-zero, it is used as a weight in
/// the averaging process.  (Usually the entries lie in the range from 0 to 1).
/// Voxels outside the boundary of the image (or outside the mask) are not
/// considered during the averaging (blurring) process, and the result after
/// blurring is weighted accordingly (if normalize=true, which it is default).
///
///    An explanation of normalization:
/// Gaussian filters are a form of weighted averaging of nearby voxels.
/// However sometimes these nearby voxels are unavailable because they either
/// lie outside the boundaries of the image, or they lie outside the mask
///   (ie. in regions where 0 <= aaafMask[iz][iy][ix] < 1.0).
/// In that case, we can "normalize" the resulting filtered image by dividing
/// this weighted average brightness of the voxels near a particular location,
/// ...by the sum of the weights which were used to calculate that average.
/// (at that location).
/// This will prevent the resulting filtered image from fading to black
/// near the image boundaries (or mask boundaries), and improve the chance
/// that subsequent features you plan to detect there are not penalized
/// simply due to lying close to the edge of the image.
///
/// @returns the "A" coefficient (determined by normalization)

template<typename Scalar>

Scalar
ApplyGauss(const int image_size[3], //!< image size in x,y,z directions
           Scalar const *const *const *aaafSource,   //!< source image (3D array)
           Scalar ***aaafDest,     //!< filtered (blurred) image stored here
           Scalar const *const *const *aaafMask,     //!< ignore voxels if aaafMask[i][j][k]==0
           Scalar const sigma[3],  //!< Gaussian sigma parameters σ_x,σ_y,σ_z
           const int truncate_halfwidth[3], //!< the filter window width
           bool normalize = true,           //!< normalize the average?
           ostream *pReportProgress = nullptr  //!< print progress to the user?
           )
{
  assert(aaafSource);
  assert(aaafDest);
  assert(sigma[0] >= 0.0);
  assert(sigma[1] >= 0.0);
  assert(sigma[2] >= 0.0);
  assert(truncate_halfwidth[0] > 0);
  assert(truncate_halfwidth[1] > 0);
  assert(truncate_halfwidth[2] > 0);

  //assert(aaafMask);
  //Allocate filters in all 3 directions.  (Later apply them sequentially.)
  Filter1D<Scalar, int> aFilter[3];
  for (int d=0; d < 3; d++)
    aFilter[d] = GenFilterGauss1D(sigma[d], truncate_halfwidth[d]);

  Scalar A;
  A = ApplySeparable(image_size, 
                     aaafSource,
                     aaafDest,
                     aaafMask,
                     aFilter,
                     normalize,
                     pReportProgress);
  return A;

} //ApplyGauss(sigma, truncate_halfwidth, ...)




/// @brief Apply a Gaussian filter (blur) to an image
/// @code h(x,y,z)=A*exp(-0.5*(x^2+y^2+z^2)/σ^2) @endcode
/// The filter must eventually be truncated.  
/// In this version, the Gaussian is isotropic
/// (in other words, it has the same width in the x,y,z directions),
/// and the caller specifies the filter truncation window width
/// (truncation_halfwidth) directly in units of voxels.
///
/// If the user specifies a mask image (if aaafMask != nullptr), then voxels whose
/// corresponding entry in the aaafMask[][][] array equals 0 are ignored,
/// If the entry in the aaafMask[][][] is non-zero, it is used as a weight in
/// the averaging process.  (Usually the entries lie in the range from 0 to 1).
/// Voxels outside the boundary of the image (or outside the mask) are not
/// considered during the averaging (blurring) process, and the result after
/// blurring is weighted accordingly (if normalize=true, which it is default).
///
///    An explanation of normalization:
/// Gaussian filters are a form of weighted averaging of nearby voxels.
/// However sometimes these nearby voxels are unavailable because they either
/// lie outside the boundaries of the image, or they lie outside the mask
///   (ie. in regions where 0 <= aaafMask[iz][iy][ix] < 1.0).
/// In that case, we can "normalize" the resulting filtered image by dividing
/// this weighted average brightness of the voxels near a particular location,
/// ...by the sum of the weights which were used to calculate that average.
/// (at that location).
/// This will prevent the resulting filtered image from fading to black
/// near the image boundaries (or mask boundaries), and improve the chance
/// that subsequent features you plan to detect there are not penalized
/// simply due to lying close to the edge of the image.
///
/// @returns the "A" coefficient (determined by normalization)

template<typename Scalar>

Scalar
ApplyGauss(const int image_size[3], //!< image size in x,y,z directions
           Scalar const *const *const *aaafSource,   //!< source image (3D array)
           Scalar ***aaafDest,     //!< filtered (blurred) image stored here
           Scalar const *const *const *aaafMask, //!< ignore voxels if aaafMask[i][j][k]==0
           Scalar sigma,                   //!< Gaussian sigma parameter σ
           int    truncate_halfwidth,      //!< the filter window width

           bool normalize = true,           //!< normalize the average?
           ostream *pReportProgress = nullptr  //!< print progress to the user?
           )
{
  Scalar afSigma[3] = {sigma, sigma, sigma};
  int afTruncateHalfwidth[3] = {truncate_halfwidth,
                                truncate_halfwidth,
                                truncate_halfwidth};
  return ApplyGauss(image_size,
                    aaafSource,
                    aaafDest,
                    aaafMask,
                    afSigma,
                    afTruncateHalfwidth,
                    normalize,
                    pReportProgress);
}





/// @brief Apply an ellipsoidal Gaussian filter (blur) to a 3D image
///
/// @code h(x,y,z)=A*exp(-0.5*((x/σ_x)^2+(y/σ_y)^2+(z/σ_z)^2) @endcode
/// The constant "A" is determined by normalization.
/// In this version, the caller can specify all 3 components of σ (σ_x,σ_y,σ_z).
/// The filter must eventually be truncated.  In this version the caller
/// specifies the filter truncation window width ("truncate_ratio") in units
/// of the Gaussian width parameter σ (as opposed to units of voxels).
///
/// If the user specifies a mask image (if aaafMask != nullptr), then voxels whose
/// corresponding entry in the aaafMask[][][] array equals 0 are ignored,
/// If the entry in the aaafMask[][][] is non-zero, it is used as a weight in
/// the averaging process.  (Usually the entries lie in the range from 0 to 1).
/// Voxels outside the boundary of the image (or outside the mask) are not
/// considered during the averaging (blurring) process, and the result after
/// blurring is weighted accordingly (if normalize=true, which it is default).
///
///    An explanation of normalization:
/// Gaussian filters are a form of weighted averaging of nearby voxels.
/// However sometimes these nearby voxels are unavailable because they either
/// lie outside the boundaries of the image, or they lie outside the mask
///   (ie. in regions where 0 <= aaafMask[iz][iy][ix] < 1.0).
/// In that case, we can "normalize" the resulting filtered image by dividing
/// this weighted average brightness of the voxels near a particular location,
/// ...by the sum of the weights which were used to calculate that average.
/// (at that location).
/// This will prevent the resulting filtered image from fading to black
/// near the image boundaries (or mask boundaries), and improve the chance
/// that subsequent features you plan to detect there are not penalized
/// simply due to lying close to the edge of the image.
///
/// @returns the "A" coefficient (determined by normalization)

template<typename Scalar>

Scalar
ApplyGauss(const int image_size[3], //!< image size in x,y,z directions
           Scalar const *const *const *aaafSource,  //!< source image (3D array)
           Scalar ***aaafDest,     //!< filtered (blurred) image stored here
           Scalar const *const *const *aaafMask,    //!< ignore voxels if aaafMask[i][j][k]==0
           const Scalar sigma[3],  //!< Gaussian sigma parameters σ_x,σ_y,σ_z
           Scalar truncate_ratio = 2.5,  //!< how many sigma before truncating?
           bool normalize = true,        //!< normalize the average?
           ostream *pReportProgress = nullptr  //!< print progress to the user?
           )
{
  Scalar A;
  int truncate_halfwidth[3];
  if (truncate_ratio > 0) {
    for (int d=0; d < 3; d++) {
      truncate_halfwidth[d] = floor(sigma[d] * truncate_ratio);
      if (truncate_halfwidth[d] < 1)
        truncate_halfwidth[d] = 1;
    }
  }

  A = ApplyGauss(image_size, 
                 aaafSource,
                 aaafDest,
                 aaafMask,
                 sigma,
                 truncate_halfwidth,
                 normalize,
                 pReportProgress);
  return A;
} // ApplyGauss(..., truncate_ratio, ...)




/// @brief Apply a spherically-symmetric Gaussian filter (blur) to a 3D image.
///
/// @note This is the most useful, popular variant of the ApplyGauss() function.
///
/// @code h(x,y,z)=A*exp(-0.5*((x^2+y^2+z^2)/σ^2) @endcode
/// The constant "A" is determined by normalization.
/// The filter must eventually be truncated.  In this version the caller
/// specifies the filter truncation window width ("truncate_ratio") in units
/// of the Gaussian width parameter σ (as opposed to units of voxels).
///
/// If the user specifies a mask image (if aaafMask != nullptr), then voxels whose
/// corresponding entry in the aaafMask[][][] array equals 0 are ignored,
/// If the entry in the aaafMask[][][] is non-zero, it is used as a weight in
/// the averaging process.  (Usually the entries lie in the range from 0 to 1).
/// Voxels outside the boundary of the image (or outside the mask) are not
/// considered during the averaging (blurring) process, and the result after
/// blurring is weighted accordingly (if normalize=true, which it is default).
///
///    An explanation of normalization:
/// Gaussian filters are a form of weighted averaging of nearby voxels.
/// However sometimes these nearby voxels are unavailable because they either
/// lie outside the boundaries of the image, or they lie outside the mask
///   (ie. in regions where 0 <= aaafMask[iz][iy][ix] < 1.0).
/// In that case, we can "normalize" the resulting filtered image by dividing
/// this weighted average brightness of the voxels near a particular location,
/// ...by the sum of the weights which were used to calculate that average.
/// (at that location).
/// This will prevent the resulting filtered image from fading to black
/// near the image boundaries (or mask boundaries), and improve the chance
/// that subsequent features you plan to detect there are not penalized
/// simply due to lying close to the edge of the image.
///
/// @returns the "A" coefficient (determined by normalization)

template<typename Scalar>

Scalar
ApplyGauss(const int image_size[3],   //!< image size in x,y,z directions
           Scalar const *const *const *aaafSource,  //!< source image (3D array)
           Scalar ***aaafDest,        //!< filtered (blurred) image stored here
           Scalar const *const *const *aaafMask,    //!< ignore voxels if aaafMask[i][j][k]==0
           Scalar sigma,                 //!< Gaussian sigma parameter σ
           Scalar truncate_ratio = 2.5,  //!< how many sigma before truncating?
           bool normalize = true,        //!< normalize the average?
           ostream *pReportProgress = nullptr  //!< print progress to the user?
           )
{
  Scalar afSigma[3] = {sigma, sigma, sigma};
  ApplyGauss(image_size,
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
/// In this version, the user specifies the filter truncation width in voxels.
/// Voxels outside the boundary of the image or outside the mask are not
/// considered during the averaging (blurring) process, and the resulting
/// after blurring is weighted accordingly.  (normalize=true by default)

template<typename Scalar>

void
ApplyDog(const int image_size[3], //!< image size in x,y,z directions
         Scalar const *const *const *aaafSource,   //!< source image (3D array)
         Scalar ***aaafDest,      //!< filtered image stored here
         Scalar const *const *const *aaafMask,     //!< ignore voxels if aaafMask[i][j][k]==0
         Scalar const sigma_a[3], //!< (a_x, a_y, a_z) in the formula above
         Scalar const sigma_b[3], //!< (b_x, b_y, b_z) in the formula above
         const int truncate_halfwidth[3],//!<(half the filter window width along x,y,z)
         Scalar *pA = nullptr, //!< Optional: report "A" (normalized coefficient) to the caller?
         Scalar *pB = nullptr, //!< Optional: report "B" (normalized coefficient) to the caller?
         ostream *pReportProgress = nullptr  //!< print progress to the user?
         )
{
  Scalar ***aaafTemp; //temporary array to store partially processed tomogram
  Scalar *afTemp;     //temporary array to store partially processed tomogram

  if (pReportProgress)
    *pReportProgress
      << " -- Attempting to allocate space for one more image.\n"
      << " -- (If this crashes your computer, find a computer with\n"
      << " --  more RAM and use \"ulimit\", OR use a smaller image.)\n";

  aaafTemp = Alloc3D<Scalar>(image_size);

  Scalar A, B;        // let the user know what A B coefficients were used

  // Convolve the original source with the 1st Gaussian
  A = ApplyGauss(image_size,
                 aaafSource,
                 aaafDest,   // <-- save result here
                 aaafMask,
                 sigma_a,
                 truncate_halfwidth,
                 true,
                 pReportProgress);

  // Convolve the original source with the 2nd Gaussian
  B = ApplyGauss(image_size,
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
  Dealloc3D(aaafTemp);
  aaafTemp = nullptr;

  // Report the A and B (normalization) coefficients to the caller?
  if (pA)
    *pA = A;
  if (pB)
    *pB = B;

} // ApplyDog()





/// @brief Apply a Laplacian-of-Gaussians filter to a 3D image.
///        The LoG filter is approximated with the following DoG filter:
/// @code
///   h(x,y,z) = scale * ( A*exp(-0.5*r_a^2) - B*exp(-0.5*r_b^2) )
///  where r_a = √((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))
///    and r_b = √((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))
///      a_x = σ_x*(1-0.5*δ), a_y = σ_y*(1-0.5*δ), a_z = σ_z*(1-0.5*δ)
///      b_x = σ_x*(1+0.5*δ), b_y = σ_y*(1+0.5*δ), b_z = σ_z*(1+0.5*δ)
///    scale = (1.0 / σ^2)
/// @endcode
/// To approximate the LoG, the δ parameter should not be larger than 0.05.
/// The constants "A" and "B" are determined by normalization.
/// In this version, the Gaussian-width (σ) can vary in the x,y,z directions.
/// (There's another version of this function where σ is a scalar.)
/// Voxels outside the boundary of the image or outside the mask are not
/// considered during the averaging (blurring) process, and the resulting
/// after blurring is weighted accordingly.  (normalize=true by default)
/// In this version, the user specifies the filter truncation width
/// (truncate_ratio) in units of σ_x, σ_y, σ_z (in the x,y,z directions).

template<typename Scalar>

void
ApplyLog(const int image_size[3], //!< source image size
         Scalar const *const *const *aaafSource,   //!< source image (3D array)
         Scalar ***aaafDest,     //!< filtered image stored here
         Scalar const *const *const *aaafMask,     //!< ignore voxels if aaafMask[i][j][k]==0
         const Scalar sigma[3],  //!< Gaussian width in x,y,z drections
         Scalar delta_sigma_over_sigma=0.02, //δ, difference in Gauss widths (approximation to LoG)
         Scalar truncate_ratio=2.5,  //!< how many sigma before truncating?
         Scalar *pA = nullptr, //!< Optional: report "A" (normalized coefficient) to the caller?
         Scalar *pB = nullptr, //!< Optional: report "B" (normalized coefficient) to the caller?
         ostream *pReportProgress = nullptr  //!< print progress to the user?
         )
{
  // "-log" approximates to the "Laplacian of a Gaussian" ("DoG") filter
  // with the difference of two Gaussians ("DoG") filter.
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
                                  std::max(sigma_a[d],
                                           sigma_b[d]));

  ApplyDog(image_size,
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
  //  Also see:
  //https://en.wikipedia.org/wiki/Blob_detection#The_difference_of_Gaussians_approach
  //
  //...where "t" = (1/2) * log_width^2
  //         because filter_mrc includes factor of (1/2) in log_width
  //   and "delta_t" = t * delta_sigma_over_sigma^2
  //               (where "delta_sigma_over_sigma" = δ in the equations above)
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

} // ApplyLog()






/// @brief Apply a Difference-of-Gaussians filter to a 3D image
///        The LoG filter is approximated with the following DoG filter:
/// @code
///   h(x,y,z) = scale * ( A*exp(-0.5*r_a^2) - B*exp(-0.5*r_b^2) )
///  where r_a = √((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))
///    and r_b = √((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))
///      a_x = σ_x*(1-0.5*δ), a_y = σ_y*(1-0.5*δ), a_z = σ_z*(1-0.5*δ)
///      b_x = σ_x*(1+0.5*δ), b_y = σ_y*(1+0.5*δ), b_z = σ_z*(1+0.5*δ)
///    scale = (1.0 / δ^2)
/// @endcode
/// To approximate the LoG, the δ parameter should not be larger than 0.05.
/// The constants "A" and "B" are determined by normalization.
/// In this version, the Gaussian width, σ, is the same in the xyz directions.
/// Voxels outside the boundary of the image or outside the mask are not
/// considered during the averaging (blurring) process, and the resulting
/// after blurring is weighted accordingly.  (normalize=true by default)

template<typename Scalar>

void
ApplyLog(const int image_size[3], //!< source image size
         Scalar const *const *const *aaafSource, //!< source image
         Scalar ***aaafDest,     //!< save results here
         Scalar const *const *const *aaafMask,  //!< ignore voxels where mask==0
         Scalar sigma,  //!< Gaussian width in x,y,z drections
         Scalar delta_sigma_over_sigma=0.02, //!< δ, difference in Gauss widths (approximation to LoG)
         Scalar truncate_ratio=2.5,  //!< how many sigma before truncating?
         Scalar *pA = nullptr, //!< Optional: report "A" (normalized coefficient) to the caller?
         Scalar *pB = nullptr, //!< Optional: report "B" (normalized coefficient) to the caller?
         ostream *pReportProgress = nullptr  //!< print progress to the user?
         )
{
  Scalar sigma_xyz[3] = {sigma, sigma, sigma};
  ApplyLog(image_size,
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




#ifndef CXX17_UNSUPPORTED

/// @brief Computes the median filter with an arbitrary "footprint" (see below),///        as defind here:
///        https://en.wikipedia.org/wiki/Median_filter
///        The "footprint" defines the nearby region over which the median
///        is calculated.  It is implemented as a vector of tuples.
///        Each tuple contains the relative coordinates of a voxel that
///        should be added to the set of voxels whose median brightness
///        will be calculated.  (A version of this function using a spherical
///        footprint is defined elsewhere.)
/// @note  This function has not been optimized for speed.
template<typename Scalar>
void
Median(vector<tuple<int,int,int> > footprint, // a list of (ix,iy,iz,b) values, one for each voxel
      const int image_size[3],                //!< size of the image in x,y,z directions
      Scalar const *const *const *aaafSource, //!< source image array
      Scalar ***aaafDest,                     //!< filter results stored here
      Scalar const *const *const *aaafMask = nullptr,   //!< optional: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0

      ostream *pReportProgress = nullptr)
{
  vector<Scalar> brightnesses(footprint.size());
  for (int iz=0; iz < image_size[2]; iz++) {
    if (pReportProgress)
      *pReportProgress << "  z = " << iz+1 << " of " << image_size[2] << endl;
    #pragma omp parallel for collapse(2)
    for (int iy=0; iy < image_size[1]; iy++) {
      for (int ix=0; ix < image_size[0]; ix++) {
        if ((aaafMask) && (aaafMask[iz][iy][ix] == 0.0))
          continue;
        size_t n_brightnesses = 0;
        aaafDest[iz][iy][ix] = 0.0;
        auto pVoxel=footprint.begin();
        auto pBrightness=brightnesses.begin();
        while(pVoxel != footprint.end()) {
          assert(pBrightness != brightnesses.end());
          int jx = std::get<0>(*pVoxel);
          int jy = std::get<1>(*pVoxel);
          int jz = std::get<2>(*pVoxel);
          int Ix = ix + jx;
          int Iy = iy + jy;
          int Iz = iz + jz;
          if ((Ix < 0) || (Ix >= image_size[0]) ||
              (Iy < 0) || (Iy >= image_size[1]) ||
              (Iz < 0) || (Iz >= image_size[2]))
            continue;
          if ((aaafMask) && (aaafMask[Iz][Iy][Ix] == 0.0))
            continue;
          *pBrightness = aaafSource[Iz][Iy][Ix];
          n_brightnesses++;
          pVoxel++;
          pBrightness++;
        }
        if (n_brightnesses > 0)
          // Now find the median of the first "n_brightnesses"
          // numbers in "brightnesses"
          aaafDest[iz][iy][ix] =
            std::nth_element(std::execution::seq,
                             brightnesses.begin(),
                             n_brightnesses/2,
                             brightnesses.begin() + n_brightnesses);
      }
    }
  }
} // Median(vector<tuple<int,int,int> > footprint,...)





/// @brief Computes the median filter of an image
///        using a spherical footprint, as defined here:
///        https://en.wikipedia.org/wiki/Median_filter
template<typename Scalar>
void
MedianSphere(Scalar radius,              //!< radius of the sphere (in voxels)
             int const image_size[3],    //!< size of the image in x,y,z directions
             Scalar const *const *const *aaafSource, //!< source image array
             Scalar ***aaafDest,             //!< filter results stored here
             Scalar const *const *const *aaafMask = nullptr,   //!< optional: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
             bool smooth_boundary = false, //!< attempt to produce a rounder smoother structure factor that varies between 0 and -1 at its boundary voxels
             ostream *pReportProgress = nullptr //!< report progress to the user
             )
{
  // footprint = a vector of (ix,iy,iz,b) values, one for each voxel
  vector<tuple<int,int,int> > footprint;

  int Ri = ceil(radius);
  for (int iz = -Ri; iz <= Ri; iz++) {
    for (int iy = -Ri; iy <= Ri; iy++) {
      for (int ix = -Ri; ix <= Ri; ix++) {
        bool add_this_voxel = false;
        Scalar r = sqrt(SQR(ix) + SQR(iy) + SQR(iz));
        if (r <= radius)
          footprint.push_back(tuple<int,int,int>(ix,iy,iz));
      } // for (int ix=0; ix < image_size[0]; ix++)
    } // for (int iy=0; iy < image_size[1]; iy++)
  } // for (int iz=0; iz < image_size[2]; iz++)

  if (pReportProgress)
    *pReportProgress << "MedianSphere(r="<<radius<<"(voxels)) progress:" <<endl;

  Median(footprint,
         image_size,
         aaafSource,
         aaafDest,
         aaafMask,
         pReportProgress);

} // MedianSphere()

#endif //#ifndef CXX17_UNSUPPORTED




/// @brief Calculate the fluctuations in intensity around every voxel lying in
///        a Gaussian-weighted ellipsoidal volume with half-width along x,y,z
///        given by sigma[].  (You can approximate a uniform sphere by setting
///        all entries in sigma[] equal to eachother, and by setting
///        template_background_exponent argument to a large number.)
///        Voxels outside the truncation window are not considered.
///        The width of the truncation window is σ*filter_truncate_ratio.
///        If a non-null aaafMask argument is supplied, then voxels from
///        aaafSource will be ignored if aaafMask[iz][iy][ix] is zero
///        (and the resulting filtered output will be normalized accordingly
///         unless the "normalize" argument is set to false).
/// @note  Generalized Gaussians can be computed by setting
///         template_background_exponent to a number which is != 2
///        (Setting it to a high value, for example, will approximate what
///         the intensity fluctuations would be within a (uniformly weighted)
///         spherical volume of radius=sigma.)  This will slow the calculation.

template<typename Scalar>

void
LocalFluctuations(const int image_size[3], //!< number of voxels in x,y,z directions
                  Scalar const *const *const *aaafSource, //!< original image
                  Scalar ***aaafDest, //!< store filtered image here (fluctuation magnitude)
                  Scalar const *const *const *aaafMask, //!< optional: if not nullptr then ignore voxel ix,iy,iz if aaafMask[iz][iy][ix]==0
                  Scalar sigma[3],  //!< radius (=sigma/√3) of neighborhooed over which we search in x,y,z directions (ellipsoidal shaped search area)
                  Scalar template_background_exponent=2, //!< exponent controlling sharpness of the (generalized) Gaussian (slow if != 2)
                  Scalar filter_truncate_ratio=2.5, //!< width over which we search is this many times larger than the gaussian width parameter (sigma)
                  bool normalize = true, //!< normalize the result?
                  ostream *pReportProgress = nullptr //!< report progress to the user?
                  )
{
  // Filter weights, w_i:
  Filter3D<Scalar, int>
    w = GenFilterGenGauss3D(sigma,
                            template_background_exponent,
                            //template_profile.halfwidth,
                            filter_truncate_ratio);

  // GenFilterGenGauss3D() creates normalized gaussians with integral 1.
  // That's not what we want for the weights, w_i:
  // The weights w_i should be 1 in the viscinity we care about, and 0 outside
  // So, we can undo this normalization by dividing all the w_i weights
  // by their maximum value max(w_i)    (at the central peak)
  // This will mean the maximum w_i is 1, and the others decay to 0, as we want
  Scalar wpeak = w.aaafH[0][0][0];
  w.MultiplyScalar(1.0 / wpeak);
  Scalar neighborhood_nvoxels = 1.0 / wpeak;
  // wpeak = (2*π*σ^2)^(-3/2) = height of a 3D Gaussian distribution of width σ.
  // = 1 / integral of an unnormalized 3D Gaussian,exp(-0.5*(x^2+y^2+z^2)/(σ^2))


  if (pReportProgress)
    *pReportProgress << " ------ Calculating the average of nearby voxels: ------\n";


  // P = original image after subtracting average nearby intensities:
  //     (equivalently, the original image after subtracting low frequencies)
  // P_dot_P = fluctuations around this local average value

  if (pReportProgress)
    *pReportProgress
      << " -- Attempting to allocate space for two more images.       --\n"
      << " -- (If this crashes your computer, find a computer with   --\n"
      << " --  more RAM and use \"ulimit\", OR use a smaller image.)   --\n";
  Scalar ***aaafP       = Alloc3D<Scalar>(image_size);
  Scalar ***aaafP_dot_P = Alloc3D<Scalar>(image_size);

  // Initially, set the contents of aaafP equal to the source image.
  for (int iz = 1; iz < image_size[2]-1; iz++)
    for (int iy = 1; iy < image_size[1]-1; iy++)
      for (int ix = 1; ix < image_size[0]-1; ix++)
        aaafP[iz][iy][ix] = aaafSource[iz][iy][ix];

  // First, let's calculate the weighted average voxel intensity in the 
  // source image
  if (template_background_exponent == 2.0) {

    // then do it the fast way with seperable (ordinary) Gaussian filters
    ApplyGauss(image_size,
               aaafSource,
               aaafP,    // <-- save result here
               aaafMask,
               sigma,
               filter_truncate_ratio,
               normalize,
               pReportProgress);
  }
  else {
    w.Apply(image_size,
            aaafSource,
            aaafP,    // <-- save result here
            aaafMask,
            normalize,
            pReportProgress);
  }

  // Subtract the average value from the image intensity, and store in P:
  for(int iz=0; iz<image_size[2]; iz++)
    for(int iy=0; iy<image_size[1]; iy++)
      for(int ix=0; ix<image_size[0]; ix++)
        aaafP[iz][iy][ix] = aaafSource[iz][iy][ix] - aaafP[iz][iy][ix];


  // Now, calculate <P_, P_>

  // Because memory is so precious, we must reuse arrays whenever we can.
  // So we will now use "P" (which used to denote voxel intensity)
  // to store the square of intensity instead
  for(int iz=0; iz<image_size[2]; iz++)
    for(int iy=0; iy<image_size[1]; iy++)
      for(int ix=0; ix<image_size[0]; ix++)
        aaafP[iz][iy][ix] *= aaafP[iz][iy][ix];

  if (pReportProgress)
    *pReportProgress << "\n"
      " ------ Calculating fluctuations around that average ------\n" << endl;

  if (template_background_exponent == 2.0) {
    // then do it the fast way with seperable (ordinary) Gaussian filters
    ApplyGauss(image_size,
               aaafP,
               aaafP_dot_P,    // <-- save result here
               aaafMask,
               sigma,
               filter_truncate_ratio,
               normalize,
               pReportProgress);
  }
  else {
    w.Apply(image_size,
            aaafP,
            aaafP_dot_P,  // <-- store result here
            aaafMask,
            normalize,
            pReportProgress);
  }

  // Now calculate "rms" (sqrt(variance))
  // Save the result in tomo_out

  // Write out RMS variance from the average of nearby voxels:
  for(int iz=0; iz<image_size[2]; iz++) {
    for(int iy=0; iy<image_size[1]; iy++) {
      for(int ix=0; ix<image_size[0]; ix++) {
        //     variance = <P_, P_>

        Scalar variance = aaafP_dot_P[iz][iy][ix];

        //Optional:
        //Compensate for dividing by w.aaafH[][][] by "wpeak" earlier.
        //This enables us to interpret variance as RMSE, (a.k.a. root-
        //mean-squared-error.  The formula above only calculates the
        //"mean" if w_i are normalized, which they were before we divided
        //them all by wpeak.  (Also: This insures that this fast method is
        //equivalent to the "SlowDebug" method commented out above.)
        //(Note: It is important for other reasons that w_i (w.aaafH[][][])
        //       were not normalized until this point.)

        variance *= wpeak;

        if (variance < 0.0)
          variance = 0.0;

        aaafDest[iz][iy][ix] = sqrt(variance);
      } //for(int ix=0; ix<image_size[0]; ix++)
    } //for(int iy=0; iy<image_size[1]; iy++)
  } //for(int iz=0; iz<image_size[2]; iz++)

  Dealloc3D(aaafP);
  aaafP = nullptr;
  Dealloc3D(aaafP_dot_P);
  aaafP_dot_P = nullptr;
} //LocalFluctuations()




/// @brief Calculate the fluctuations in intensity around every voxel lying in
///        a Gaussian-weighted ellipsoidal volume with half-width along x,y,z
///        given by radius[].  (You can approximate a uniform sphere by setting
///        all entries in radius[] equal to eachother, and by setting
///        template_background_exponent argument to a large number.)
///        Voxels outside the truncation window are not considered.
///        The width of the truncation window is σ*filter_truncate_ratio.
///        If a non-null aaafMask argument is supplied, then voxels from
///        aaafSource will be ignored if aaafMask[iz][iy][ix] is zero
///        (and the resulting filtered output will be normalized accordingly
///         unless the "normalize" argument is set to false).
///
/// What's the relationship between the Gaussian width parameter sigma (σ),
/// and the radius (r)?
///
/// The number of voxels belonging to a Gaussian is equal to the integral
/// of the unnormalized Gaussian (which is exp(-0.5*(x^2+y^2+z^2)/(σ^2)))
/// integrated over all space, resulting in (2*π*σ^2)^(3/2).
/// This is effectively the volume of the Gaussian (denoted "v" below).
/// Let's approximate this by a sphere of radius r containing the same volume v
/// v = (4/3)*π*r^3 = the number of voxels in the sphere
/// solving for r:
/// @code
/// <-> r = (3/(4π))^(1/3)  *  v^(1/3)
///       substituting v = (2*π*σ^2)^(3/2)
/// --> r = σ * 3^(1/3) * 2^(1/2-2/3) * π^(1/2 - 1/3)
///       = σ * (9π/2)^(1/6)
///      ~= 1.5549880806696572 σ
/// @endcode
/// @note  Generalized Gaussians can be computed by setting
///         template_background_exponent to a number which is != 2
///        (Setting it to a high value, for example, will approximate what
///         the intensity fluctuations would be within a (uniformly weighted)
///         spherical volume of radius=sigma.
///        If you do this, you will have to multiply your radius arguments by
///        1.5549880806696572 beforehand to compensate.

template<typename Scalar>
void
LocalFluctuationsByRadius(const int image_size[3], //!< number of voxels in x,y,z directions
                          Scalar const *const *const *aaafSource, //!< original image
                          Scalar ***aaafDest, //!< store filtered image here (fluctuation magnitude)
                          Scalar const *const *const *aaafMask, //!< optional: if not nullptr then ignore voxel ix,iy,iz if aaafMask[iz][iy][ix]==0
                          const Scalar radius[3],  //!< radius (=sigma/√3) of neighborhooed over which we search in x,y,z directions (ellipsoidal shaped search area)
                          Scalar template_background_exponent=2, //!< exponent controlling sharpness of the (generalized) Gaussian (slow if != 2)
                          Scalar filter_truncate_ratio=2.5, //!< width over which we search is this many times larger than the gaussian width parameter (sigma)
                          bool normalize = true, //!< normalize the result?
                          ostream *pReportProgress = nullptr //!< report progress to the user?
                        )
{
  Scalar sigma[3];

  Scalar ratio_r_over_sigma = pow((9.0/2)*M_PI, 1.0/6);
  //Scalar ratio_r_over_sigma = sqrt(3.0);  //earlier crude estimate (don't use)

  for (int d = 0; d < 3; d++)
    sigma[d] = (radius[d] / ratio_r_over_sigma);

  LocalFluctuations(image_size,
                    aaafSource,
                    aaafDest,
                    aaafMask,
                    sigma,
                    template_background_exponent,
                    filter_truncate_ratio,
                    normalize,
                    pReportProgress);
}




} //namespace visfd



#endif //#ifndef _FILTER3D_HPP
