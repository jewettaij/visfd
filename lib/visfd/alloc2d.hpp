///   @file alloc2d.hpp
///   @brief functions that allocate and deallocate 2D arrays
///   @author Andrew Jewett
///   @date 2021-7-13

#ifndef _ALLOC2D_HPP
#define _ALLOC2D_HPP


namespace visfd {



/// @brief    A function for allocating C-style 2-dimensional arrays of data
///           (contiguous in memory, in row-major format.)
/// @returns  The array after allocation.  (Suppose A = Alloc2D(image_size).
///           Then the contents of A can be accessed using A[iy][ix].)
/// @details  This array can be deallocated using Dealloc2D().

template<typename Entry>
void Alloc2D(int size[2]           //!< size of the array (number of rows)
             )
{
  Entry **aaX = new Entry* [nrows];
  aaX[0] = new Entry [nrows * ncols];  // 1D C array (contiguous memory)
  for(size_t iy=0; iy<nrows; iy++)
    aaX[iy] = aaX[0] + iy*ncols;
  return aaX;
}



/// @brief
/// This function is the corresponding way to dellocate arrays
/// that were created using Alloc2D()
template<typename Entry>
void Dealloc2D(Entry **aaX       //!< pointer to a 2D C-style array
               )
{
  if (aaX) {
    delete [] aaX[0];
    delete [] aaX;
  }
}



} //namespace visfd


#endif //#ifndef _ALLOC2D_HPP
