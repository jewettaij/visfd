///   @file alloc3d.hpp
///   @brief functions that allocate and deallocate 3D arrays
///   @author Andrew Jewett
///   @date 2021-7-13


#ifndef _ALLOC3D_HPP
#define _ALLOC3D_HPP

#include <cstddef>

namespace visfd {



/// @brief    A function for allocating C-style 3-dimensional arrays of data
///           which are contiguous in memory and in row-major order.
/// @returns  The array after allocation.  (Suppose A = Alloc3D(image_size);
///           Then the contents of A can be accessed using A[iz][iy][ix].)
/// @details  In addition to the array contents, additional space is also
///           allocated to store the internal arrays of pointers
///           A[i] and A[i][j], for all i,j.  All of this
///           additional memory is deleted when Dealloc3D() is invoked.

template<typename Entry, typename Integer>
Entry ***
Alloc3D(const Integer size[3]  //!< size of the array in x,y,z directions
        )
{
  // We need to allocate space to store both the entries of the 3-D array
  // as well as the internal pointers in the pointer->pointer-pointer structure.
  // aC points to a contiguous block of memory large enough to hold everything.
  std::byte *ac = new std::byte [sizeof(Entry**) * (size[2]) +
                                 sizeof(Entry*)  * (size[2]*size[1]) +
                                 sizeof(Entry)   * (size[2]*size[1]*size[0])];

  // traditional 3-D C pointer-to-pointer-to-pointer
  Entry ***aaaX = reinterpret_cast<Entry***>(ac);

  // Offsets into the "ac" array:
  //
  // 1) bytes  0 ... (sizeof(Entry**)*size[2]) - 1
  //               -> store the a[iz] pointer entries
  //
  // 2) bytes  sizeof(Entry**)*size[2] ...
  //           sizeof(Entry**)*size[2] + sizeof(Entry*)*size[2]*size[1] -1 bytes
  //               -> store the a[iz][iy] pointer entries
  //
  // 3) The remaining bytes store the actual entries in the array

  // The aaafX[][] pointers should point to to the 
  // appropriate locations within the contents array.
  for(Integer iz=0; iz<size[2]; iz++) {
    aaaX[iz] = reinterpret_cast<Entry**>(ac + sizeof(Entry**) * size[2] +
                                              sizeof(Entry*)  * size[1] * iz);
    for(Integer iy=0; iy<size[1]; iy++) {
      aaaX[iz][iy] = reinterpret_cast<Entry*>(ac +
                                              sizeof(Entry**) * size[2] +
                                              sizeof(Entry*)*size[2]*size[1]+
                                              sizeof(Entry)*(iz*size[0]*size[1]+
                                                             iy*size[0]));
    }
  }

  return aaaX;

} // Alloc3d()




/// @brief   This function is the corresponding way to dellocate arrays
///          that were created using Alloc3D().

template<typename Entry>
void Dealloc3D(Entry ***aaaX   //!< a 3-D array created by Alloc3D()
               )
{
  if (aaaX) {
    delete [] reinterpret_cast<std::byte*>(aaaX);
    // Note: This also seems to work:
    //   delete [] aaaX;
  }
}





} //namespace visfd


#endif //#ifndef _ALLOC3D_HPP
