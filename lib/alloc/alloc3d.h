#ifndef _ALLOC3D_H
#define _ALLOC3D_H


// Alloc3D() is a function for allocating 3-dimensional arrays of data
// (contiguous in memory, in row-major format.)
// (The functions in this file are used to allocate arrays for storing
//  tomographic data, and also precomputed filter-weights.)


template<class Entry, class Integer>
void Alloc3D(Integer const size[3], 
             Entry **paX,       // <--pointer to 1-D contiguous-memory array
             Entry ****paaaX) { // <--pointer to 3-D multidimensional array
  if (! paX)
    return;

  // Allocate a 3-dimensional table row-major order
  // Optional: Also allocate a conventional 3-dimensional
  //           pointer-to-a-pointer-to-a-pointer data structure (aaaX), that
  //           you can use to access the contents using aaaX[k][j][i] notation.
  *paX = new Entry [size[0] * 
                    size[1] * 
                    size[2]];

  if (! paaaX)
    return;

  *paaaX = new Entry** [size[2]];

  // (The aaafX[][] pointers should point to to the 
  //  appropriate locations within the afX array.)
  for(Integer iz=0; iz<size[2]; iz++) {
    (*paaaX)[iz] = new Entry* [size[1]];
    for(Integer iy=0; iy<size[1]; iy++) {
      (*paaaX)[iz][iy] = &((*paX)[iz*size[0]
                                    *size[1] + 
                                  iy*size[0]]);
    }
  }
}


template<class Entry, class Integer>
void Dealloc3D(Integer const size[3], 
               Entry **paX, 
               Entry ****paaaX) {
  if (paaaX && *paaaX) {
    for(Integer iz=0; iz<size[2]; iz++)
      if ((*paaaX)[iz])
        delete [] (*paaaX)[iz];
    delete [] (*paaaX);
    *paaaX = NULL;
  }
  if (paX && *paX) {
    delete [] *paX;
    *paX = NULL;
  }
}



#endif //#ifndef _ALLOC3D_H
