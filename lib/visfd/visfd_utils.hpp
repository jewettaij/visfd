#ifndef _VISFD_UTILS_HPP
#define _VISFD_UTILS_HPP

#include <cassert>
#include <limits>
using namespace std;
#include <alloc3d.hpp>


#include <eigen3_simple.hpp>  //defines namespace selfadjoint_eigen3
#include <lin3_utils.hpp>  //defines namespace selfadjoint_eigen3
using namespace visfd::selfadjoint_eigen3;




namespace visfd {




template<class Scalar >
static inline Scalar SQR(Scalar x) { return x*x; }

template<class Scalar >
static inline Scalar MIN(Scalar a, Scalar b) { return ((a<=b) ? a : b); }

template<class Scalar >
static inline Scalar MAX(Scalar a, Scalar b) { return ((a>=b) ? a : b); }

template <class Scalar>
static int SGN(Scalar val) {
  return (static_cast<Scalar>(0) < val) - (val < static_cast<Scalar>(0));
}




/// @brief invert a permutation
template<class T, class Integer>
void
invert_permutation(const vector<Integer>& p,
                   vector<T>& p_inv)
{
  p_inv.resize(p.size());
  for (Integer i = 0; i < p.size(); i++)
    p_inv[p[i]] = i;
}



/// @brief apply a permutation to a std::vector in-place
template<class T, class Integer>
void
apply_permutation(const vector<Integer>& p,
                  vector<T>& v)
{
  assert(p.size() == v.size());
  int n = p.size();
  vector<T> v_copy(v);
  for (Integer i = 0; i < n; i++) {
    Integer j = p[i];
    assert((0 <= j) && (j < n));
    v[i] = v_copy[j];
  }
}




/// @brief  Calculate the hessian of a 3D image at a particular position 
///         ix, iy, iz, from the differences between neighboring voxels.
///         It is the responsibility of the caller to smooth the source image
///         (if necessary) before this function is called.
/// @note   You must insure that ix,iy,iz (and their neighbors!) lie within
///         the boundaries of the image.  THERE IS NO BOUNDS CHECKING.
template<class Scalar>
void
CalcHessianFiniteDifferences(Scalar const *const *const *aaafSource, //!< source image
                             int ix, int iy, int iz, //!< location in the image where you wish to calculate the Hessian
                             Scalar (*hessian)[3]  //!<store resulting 3x3 matrixhere
                             )
{
  assert(aaafSource);
  assert(hessian);

  hessian[0][0] = (aaafSource[iz][iy][ix+1] + 
                   aaafSource[iz][iy][ix-1] - 
                   2*aaafSource[iz][iy][ix]);
  hessian[1][1] = (aaafSource[iz][iy+1][ix] + 
                   aaafSource[iz][iy-1][ix] - 
                   2*aaafSource[iz][iy][ix]);
  hessian[2][2] = (aaafSource[iz+1][iy][ix] + 
                   aaafSource[iz-1][iy][ix] - 
                   2*aaafSource[iz][iy][ix]);

  hessian[0][1] = 0.25 * (aaafSource[iz][iy+1][ix+1] + 
                          aaafSource[iz][iy-1][ix-1] - 
                          aaafSource[iz][iy-1][ix+1] - 
                          aaafSource[iz][iy+1][ix-1]);
  hessian[1][0] = hessian[0][1];

  hessian[1][2] = 0.25 * (aaafSource[iz+1][iy+1][ix] + 
                          aaafSource[iz-1][iy-1][ix] - 
                          aaafSource[iz-1][iy+1][ix] - 
                          aaafSource[iz+1][iy-1][ix]);
  hessian[2][1] = hessian[1][2];

  hessian[2][0] = 0.25 * (aaafSource[iz+1][iy][ix+1] + 
                          aaafSource[iz-1][iy][ix-1] - 
                          aaafSource[iz+1][iy][ix-1] - 
                          aaafSource[iz-1][iy][ix+1]);
  hessian[0][2] = hessian[2][0];
} //CalcHessianFiniteDifferences()



/// @brief  Calculate the hessian of a 3D image at a particular position 
///         ix, iy, iz, from the differences between neighboring voxels.
///         It is the responsibility of the caller to smooth the source image
///         (if necessary) before this function is called.
/// @note   This version of the function provides bounds checking for ix,iy,iz.
///         (At the boundaries of the image, it will substitute the voxel
///          brightnesses from the nearest neighbor.)
/// @note   This function was not intended for public use.

template<class Scalar>
void
CalcHessianFiniteDifferences(Scalar const *const *const *aaafSource, //!< source image
                             int ix, int iy, int iz, //!< location in the image where you wish to calculate the Hessian
                             Scalar (*hessian)[3],  //!< store resulting 3x3 matrixhere
                             const int image_size[3] //!< number of voxels in xyz directions (for bounds checking)
                             )
{
  assert(aaafSource);
  assert(hessian);

  int _ix = ix;
  int _iy = iy;
  int _iz = iz;

  if (image_size) {
    assert(image_size[0] >= 3);
    assert(image_size[1] >= 3);
    assert(image_size[2] >= 3);

    if (_ix == 0)
      _ix++;
    else if (_ix == image_size[0]-1)
      _ix--;

    if (_iy == 0)
      _iy++;
    else if (_iy == image_size[1]-1)
      _iy--;

    if (_iz == 0)
      _iz++;
    else if (_iz == image_size[2]-1)
      _iz--;
  }

  CalcHessianFiniteDifferences(aaafSource,
                               _ix, _iy, _iz,
                               hessian);
}



/// @brief  Calculate the hessian of a 3D image at a particular position 
///         ix, iy, iz, from the differences between neighboring voxels.
///         It is the responsibility of the caller to smooth the source image
///         (if necessary) before this function is called.
/// @note   This version of the function provides bounds checking for ix,iy,iz.
///         (At the boundaries of the image, it will substitute the voxel
///          brightnesses from the nearest neighbor.)
/// @note   This function was not intended for public use.

template<class Scalar>
void
CalcGradientFiniteDifferences(Scalar const *const *const *aaafSource, //!< source image
                              int ix, int iy, int iz, //!< location in the image where you wish to calculate the Hessian
                              Scalar *gradient,  //!< store resulting 3x3 matrixhere
                              const int image_size[3] //!< number of voxels in xyz directions (for bounds checking)
                              )
{
  assert(aaafSource);
  assert(gradient);

  int _ix = ix;
  int _iy = iy;
  int _iz = iz;

  if (image_size) {
    assert(image_size[0] >= 3);
    assert(image_size[1] >= 3);
    assert(image_size[2] >= 3);

    if (_ix == 0)
      _ix++;
    else if (_ix == image_size[0]-1)
      _ix--;
    if (_iy == 0)
      _iy++;
    else if (_iy == image_size[1]-1)
      _iy--;
    if (_iz == 0)
      _iz++;
    else if (_iz == image_size[2]-1)
      _iz--;
  }

  gradient[0]=0.5*(aaafSource[_iz][_iy][_ix+1] - 
                   aaafSource[_iz][_iy][_ix-1]);
  gradient[1]=0.5*(aaafSource[_iz][_iy+1][_ix] - 
                   aaafSource[_iz][_iy-1][_ix]);
  gradient[2]=0.5*(aaafSource[_iz+1][_iy][_ix] - 
                   aaafSource[_iz-1][_iy][_ix]);
}




/// @brief 
/// Compute a weighted average of the entries in a 3-dimensional array.
///         This function was not intended for public use.
/// @param  array_size  an array of 3 integers storing the size of the array
/// @param  aaafH       the array whose entries will be averaged
/// @param  aaafW       an optional array of weights (for weighted averages)
/// @return the weighted average


template<class Scalar, class Integer>
Scalar AverageArr(Integer const array_size[3],
                  Scalar const *const *const *aaafH,
                  Scalar const *const *const *aaafW = NULL) 
{
  Scalar total = 0.0;
  Scalar denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        Scalar h = aaafH[iz][iy][ix];
        if (aaafW) {
          Scalar w = aaafW[iz][iy][ix];
          h *= w;
          denom += w;
        }
        else
          denom += 1.0;
        total += h;
      }
    }
  }
  return static_cast<Scalar>(total / denom);
} //void AverageArr()



/// @brief 
/// Compute a weighted average of the square of entries in a 3D array.
///         This function was not intended for public use.
/// @param  array_size  an array of 3 integers storing the size of the array
/// @param  aaafH       the array whose entries will be squared and averaged
/// @param  aaafW       an optional array of weights (for weighted averages)
/// @return the weighted average of the squared entries in the array

template<class Scalar, class Integer>
static
Scalar _AveSqrArr(Integer const array_size[3],
                  Scalar const *const *const *aaafH,
                  Scalar const *const *const *aaafW = NULL) 
{
  Scalar total = 0.0;
  Scalar denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        Scalar h = aaafH[iz][iy][ix];
        h *= h;

        if (aaafW) {
          Scalar w = aaafW[iz][iy][ix];
          h *= w;
          denom += w;
        }
        else
          denom += 1.0;
        total += h;
      }
    }
  }
  return static_cast<Scalar>(sqrt(total / denom));
} //void _AveSqrArr()




/// @brief 
/// Compute the weighted standard deviation of entries in a 3-dimensional array.
///         This function was not intended for public use.
/// @param  array_size  an array of 3 integers storing the size of the array
/// @param  aaafH       we will compute the standard deviation of entries
///                     in this array
/// @param  aaafW       an optional array of weights (for weighted averages)
/// @return the weighted standard deviation

template<class Scalar, class Integer>
Scalar StdDevArr(Integer const array_size[3],
                 Scalar const *const *const *aaafH,
                 Scalar const *const *const *aaafW = NULL) 
{
  Scalar ave = AverageArr(array_size, aaafH, aaafW);
  Scalar total = 0.0;
  Scalar denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        Scalar h = (aaafH[iz][iy][ix] - ave);
        h *= h;

        if (aaafW) {
          Scalar w = aaafW[iz][iy][ix];
          h *= w;
          denom += w;
        }
        else
          denom += 1.0;
        total += h;
      }
    }
  }
  return static_cast<Scalar>(sqrt(total / denom));
} //void StdDevArr()



/// @brief Compute a (weighted) sum of the entries in a 3D array.
///         This function was not intended for public use.
/// @param array_size  an array of 3 integers storing the size of the array
/// @param  aaafH       we will compute the sum of the squares of the 
///                     entries in this array
/// @param  aaafW       an optional array of weights (which multiply the
///                     the corresponding entry value)
/// @return the (weighted) sum

template<class Scalar, class Integer>
static
Scalar _SumArr(Integer const array_size[3],
               Scalar const *const *const *aaafH,
               Scalar const *const *const *aaafW = NULL) 
{
  Scalar total = 0.0;
  Scalar denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        Scalar h = aaafH[iz][iy][ix];
        total += h;
        if (aaafW)
          h *= aaafW[iz][iy][ix];
      }
    }
  }
  return static_cast<Scalar>(total);
} //void _SumArr()



/// @brief
/// Compute a (weighted) sum of the squares of the entries in a 3D array
///         This function was not intended for public use.
/// @param array_size  an array of 3 integers storing the size of the array
/// @param  aaafH       we will compute the sum of the squares of the 
///                     entries in this array
/// @param  aaafW       an optional array of weights (which multiply the
///                     the corresponding squared entry value)
/// @return the (weighted) sum of squares

template<class Scalar, class Integer>
static
Scalar _SumSqrArr(Integer const array_size[3],
                  Scalar const *const *const *aaafH,
                  Scalar const *const *const *aaafW = NULL) 
{
  Scalar total = 0.0;
  Scalar denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        Scalar h = aaafH[iz][iy][ix];
        h *= h;
        if (aaafW)
          h *= aaafW[iz][iy][ix];
        total += h;
      }
    }
  }
  return static_cast<Scalar>(total);
} //void _SumSqrArr()



/// @brief  Add a scalar offset to all of the entries in a 3D array.
///         This function was not intended for public use.
/// @param  offset       the number that will be added to all array entries
/// @param  array_size  an array of 3 integers storing the size of the array
/// @param  aaafH       the 3D array containing the entries to be modified

template<class Scalar, class Integer>
void AddScalarArr(Scalar offset,
                  Integer const array_size[3],
                  Scalar ***aaafH)
{
  for (Integer iz = 0; iz < array_size[2]; iz++)
    for (Integer iy = 0; iy < array_size[1]; iy++)
      for (Integer ix = 0; ix < array_size[0]; ix++)
        aaafH[iz][iy][ix] += offset;
}



/// @brief  Multiply all the entries in a 3D array by a scalar.
///         This function was not intended for public use.
/// @param  scale       the number that will be multiplied by the array entries
/// @param  array_size  an array of 3 integers storing the size of the array
/// @param  aaafH       the 3D array containing the entries to be modified

template<class Scalar, class Integer>
void MultiplyScalarArr(Scalar scale,
                       Integer const array_size[3],
                       Scalar ***aaafH)
{
  for (Integer iz = 0; iz < array_size[2]; iz++)
    for (Integer iy = 0; iy < array_size[1]; iy++)
      for (Integer ix = 0; ix < array_size[0]; ix++)
        aaafH[iz][iy][ix] *= scale;
}





/// @brief  Find the minimum entries in a 3D array.
///         This function was not intended for public use.
/// @param  array_size  an array of 3 integers storing the size of the array
/// @param  aaafI       the 3D array containing the entries
/// @param  aaafMask    (optional) If aaafMask[i][j][k]==0 ignore this entry
/// @param  afLocation  (optional) report the location of the minima (ix,iy,iz)
/// @return the minimum entry in the aaafI array (not ignored by the mask)
///         (or std::numeric_limits::infinty(), if the aaafMask array has no non-zero entries)

template<class Scalar, class Integer>
static
Scalar _MinArr(Integer const array_size[3],
               Scalar const *const *const *aaafI,
               Scalar const *const *const *aaafMask = NULL,
               Integer *afLocation = NULL)
{
  Scalar min_I = std::numeric_limits<Scalar>::infinity(); //return this suspicious value if mask is empty
  bool first_iter = true;
  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        if ((aaafMask) and (aaafMask[iz][iy][ix] == 0.0))
          continue;
        if ((aaafI[iz][iy][ix] < min_I) || (first_iter)) {
          min_I = aaafI[iz][iy][ix];
          if (afLocation) {
            afLocation[0] = ix;
            afLocation[1] = iy;
            afLocation[2] = iz;
          }
        }
        first_iter = false;
      }
    }
  }
  return min_I;
} //void _MinArr()


/// @brief  Find the maximum and maximum entries in a 3D array.
///         This function was not intended for public use.
/// @param  array_size  an array of 3 integers storing the size of the array
/// @param  aaafI       the 3D array containing the entries
/// @param  aaafMask    (optional) If aaafMask[i][j][k]==0 ignore this entry
/// @param  afLocation  (optional) report the location of the maxima (ix,iy,iz)
/// @return the minimum entry in the aaafI array (not ignored by the mask)
///         (or std::numeric_limits::infinty(), if the aaafMask array has no non-zero entries)

template<class Scalar, class Integer>
static
Scalar _MaxArr(Integer const array_size[3],
               Scalar const *const *const *aaafI,
               Scalar const *const *const *aaafMask = NULL,
               Integer *afLocation = NULL)
{
  Scalar max_I = -std::numeric_limits<Scalar>::infinity(); //return this suspicious value if mask is empty
  bool first_iter = true;
  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        if ((aaafMask) and (aaafMask[iz][iy][ix] == 0.0))
          continue;
        if ((aaafI[iz][iy][ix] < max_I) || (first_iter)) {
          max_I = aaafI[iz][iy][ix];
          if (afLocation) {
            afLocation[0] = ix;
            afLocation[1] = iy;
            afLocation[2] = iz;
          }
        }
        first_iter = false;
      }
    }
  }
  return max_I;
} //void _MaxArr()



/// @brief  Generate a histogram of the values contained in a 3D numeric array
/// @param  paHistX     store entry values corresponding to the center of each bin here
/// @param  paHistY     store the number of entries that fall into that bin here
/// @param  nbins       specify number of bins (if positive)
/// @param  bin_width   alternatively, (if nbins<=0), specify the width of each bin
/// @param  aaafI       the 3D array containing the entries
/// @param  aaafMask    (optional) If aaafMask[i][j][k]==0 ignore this entry

template<class Scalar, class Integer>
void
HistogramArr(Scalar **paHistX,  
             size_t  **paHistY,
             Integer &nbins,
             Scalar &bin_width,
             Integer const array_size[3],
             Scalar ***aaafI,
             Scalar const *const *const *aaafMask = NULL) 
{
  Scalar hmin = 0.0;
  Scalar hmax = -1.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        if (hmin > hmax) {
          hmin = aaafI[iz][iy][ix];
          hmax = aaafI[iz][iy][ix];
        }
        if (aaafI[iz][iy][ix] < hmin)
          hmin = aaafI[iz][iy][ix];
        if (aaafI[iz][iy][ix] > hmax)
          hmax = aaafI[iz][iy][ix];
      }
    }
  }

  if (nbins > 1) {
    bin_width = (hmax - hmin) / (nbins - 1);
  }
  else if (bin_width > 0) {
    nbins = 1 + ceil((hmax - hmin)/bin_width);
  }

  assert(paHistX);
  *paHistX = new Scalar [nbins];
  for (int i=0; i<nbins; i++)
    (*paHistX)[i] = hmin + i*bin_width;

  assert(paHistY);
  *paHistY = new size_t [nbins];
  for (int i=0; i < nbins; i++)
    (*paHistY)[i] = 0;
  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        Scalar h = aaafI[iz][iy][ix];
        if ((! aaafMask) || (aaafMask[iz][iy][ix] != 0.0)) {
          Scalar h = aaafI[iz][iy][ix];
          Integer i = floor( (h - hmin)/bin_width );
          assert((0 <= i) && (i < nbins));
          (*paHistY)[i]++;
        }
      }
    }
  }
} //void HistogramArr()



} //namespace visfd



#endif //#ifndef _VISFD_UTILS_HPP
