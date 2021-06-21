///   @file resample.hpp
///   @brief Some code for resizing (changing the resolution of) an image.
///          As of 2021-6-19, only binning is supported.
///          ("Binning" is a way to reduce the resolution of an image in each
///           direction (by a factor of n).  Each new voxel's brightness is the
///           average of the corresponding n^3 voxels in the original image.)
///          Perhaps later, more sophisticated image resampling methods
///          methods will be offerred.
///   @author Andrew Jewett
///   @date 2021-6-19

#ifndef _RESAMPLE_HPP
#define _RESAMPLE_HPP

#include <cstring>
#include <cassert>
#include <cmath>
#include <ostream>
#include <string>
#include <sstream>
using namespace std;
#include <err_visfd.hpp>  // defines the "VisfdErr" exception type
#include <alloc2d.hpp>    // defines Alloc2D() and Dealloc2D()
#include <alloc3d.hpp>    // defines Alloc3D() and Dealloc3D()
#include <filter1d.hpp>   // defines "Filter1D" (used in ApplySeparable())
#include <filter3d.hpp>   // defines common 3D image filters


namespace visfd {


template<typename Scalar, typename Integer>
/// @brief  Reduce the resolution of the image in the x, y, z directions
///         by a factor of bin_size[0], bin_size[1], bin_size[2],
///         ...where bin_size[d] = floor(size_source[d] / size_dest[d]),
///         (where "d" is an integer from 0 to 2).
///         Each new voxel's brightness is the average brightness of the
///         corresponding bin_size[0]*bin_size[1]*bin_size[2] voxels
///         from the original source image.
/// @note   If the size_dest[d] is not a divisor of size_source[d],
///         (ie. if size_source[d] / size_dest[d] is not an integer), then
///         the voxels from the source located after size_dest[d] * bin_size[d]
///         (in the d'th direction) will be discarded.
/// @note   The caller can shift the location of the binning window in the
///         x,y,z direction by passing an optional "offset" argument.  If not
///         NULL, this argument should be an array of size 3.  The d'th entry
///         in this array should be an integer between 0 and bin_size[d]-1.
/// @param size_source contains size of the source image (in the x,y,z directions)
/// @param size_dest contains size of the source image (in the x,y,z directions)
/// @param aaafSource[][][] is the source array (source image)
/// @param aaafDest[][][] will store the (reduced-size) image after binning
///        This array is assumed to have been preallocated.
void BinArray3D(Integer const size_source[3],
                Integer const size_dest[3],
                Scalar const *const *const *aaafSource,
                Scalar ***aaafDest,
                Integer const *offset = nullptr)
{
  Integer bin_size[3];
  for (int d = 0; d < 3; d++) {
    assert((size_source[d] > 0) && (size_dest[d] > 0));
    bin_size[d] = size_source[d] / size_dest[d];
    if (offset && ((offset[d] >= bin_size[d]) || (offset[d] < 0))) {
      stringstream err_msg;
      err_msg << "Error in BinArray3D(): offset[" << d
              << "] should lie between 0 and floor("
              << size_source[d] << " / " << size_dest[d] << ")\n";
      throw VisfdErr(err_msg.str());
    }
  }
  for (Integer Iz=0; Iz < size_dest[2]; Iz++) {
    for (Integer Iy=0; Iy < size_dest[1]; Iy++) {
      for (Integer Ix=0; Ix < size_dest[0]; Ix++) {
        Scalar sum = 0.0;
        // Average the entries in the voxels from the source image
        // which will be merged into this single voxel in the destination image.
        for (Integer dz=0; dz < bin_size[2]; dz++) {
          for (Integer dy=0; dy < bin_size[1]; dy++) {
            for (Integer dx=0; dx < bin_size[0]; dx++) {
              Integer ix = Ix*bin_size[0] + dx;
              Integer iy = Iy*bin_size[1] + dy;
              Integer iz = Iz*bin_size[2] + dz;
              if (offset) {
                ix += offset[0];
                iy += offset[1];
                iz += offset[2];
              }
              assert((0 <= ix) && (ix < size_source[0]));
              assert((0 <= iy) && (iy < size_source[1]));
              assert((0 <= iz) && (iz < size_source[2]));
              sum += aaafSource[iz][iy][ix];
            }
          }
        }
        Scalar ave = sum / (bin_size[0] * bin_size[1] * bin_size[2]);
        aaafDest[Iz][Iy][Ix] = ave;
      }
    }
  }
}



} //namespace visfd


#endif //#ifndef _RESAMPLE_HPP
