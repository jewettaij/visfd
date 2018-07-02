#ifndef _MRC_SIMPLE_H
#define _MRC_SIMPLE_H

// "MrcSimple": a class to read and write tomographic data.
//              (stored in .MRC and .REC files)
//
// Tomographic data is internally stored as arrays of floats, however
// this program can read MRC/REC files using several other numeric formats 
// as well (including signed and unsigned 8-bit and 16-bit integers).
// Note: As of 2015-4-20, file data is assumed to be in row-major format.
//       http://en.wikipedia.org/wiki/Row-major_order
//       (IE The "mapC", "mapR", and "mapS" (mapCRS[]) header data is ignored.)
// Note: As of 2015-4-20, this software does not attempt to detect 
//       or convert big-endian or little-endian numbers.
//       To be on the safe side, compile this software on hardware which is
//       compatible with the hardware which was used to write the MRC file.
//       (IE intel CPUs, or ARM CPUs using the same endian settings.)
//
// Documentation on the MRC file format is available here:
//http://bio3d.colorado.edu/imod/doc/mrc_format.txt
//http://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/image-processing-software/#image
//http://ami.scripps.edu/software/mrctools/mrc_specification.php
//http://bio3d.colorado.edu/imod/betaDoc/libhelp/mrcfiles.html#TOP



#include <cstring>
#include <iostream>
using namespace std;
#include "mrc_header.h"


// "MrcSimple": a class to read and write tomographic data.
//              (stored in .MRC and .REC files)
//
// Tomographic data is internally stored as arrays of floats, however
// this program can read MRC/REC files using several other numeric formats 
// as well (including signed and unsigned 8-bit and 16-bit integers).
// Note: As of 2015-4-20, file data is assumed to be in row-major format.
//       http://en.wikipedia.org/wiki/Row-major_order
//       (IE The "mapC", "mapR", and "mapS" (mapCRS[]) header data is ignored.)
// Note: As of 2015-4-20, this software does not attempt to detect 
//       or convert big-endian or little-endian numbers.
//       To be on the safe side, compile this software on hardware which is
//       compatible with the hardware which was used to write the MRC file.
//       (IE intel CPUs, or ARM CPUs using the same endian settings.)
//
// Documentation on the MRC file format is available here:
//http://bio3d.colorado.edu/imod/doc/mrc_format.txt
//http://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/image-processing-software/#image
//http://ami.scripps.edu/software/mrctools/mrc_specification.php
//http://bio3d.colorado.edu/imod/betaDoc/libhelp/mrcfiles.html#TOP



class MrcSimple {
 public:
  MrcHeader header;     //<-- contains the size and dimensions of the tomogram:
                        //    (header.nvoxels[] stores the number of voxels
                        //                          in the x,y,z directions.
                        //     header.cellA[] stores the physical size of
                        //                       the tomographic reconstruction)

  float *afI;           //<--pointer to a (contiguous) 1-D array containing 
                        //   the entire tomogram. The data at x,y,z is located
                        //    at  afI[ix + iy*sizez + iz*sizex*sizey]

  float ***aaafI;       //<-- Useful if you prefer to use 3-index notation:
                        //    aafI[iz][iy][ix] =
                        //     afI[ix + iy*sizez + iz*sizex*sizey]
 private:
  void Alloc();  // Allocate memory for the tomogram.
                 // Allocates and initializes afI and aaafI. Invoke
                 // after using ReadHeader(), or after modifying the header.
  void Dealloc();

 public:
  void Read(string mrc_file_name,  //Read an .MRC/.REC file
            bool rescale=true,     //Optional: rescale densities from 0 to 1?
            float ***aaafMask=NULL);//(also chooses header.uses_signed_bytes)

  void Write(string mrc_file_name); //Write an .MRC/.REC file

  void Read(istream& mrc_file,      //Read an .MRC/.REC file (input stream)
            bool rescale=true,
	    float ***aaafMask=NULL);
  void Write(ostream& mrc_file); //Write an .MRC/.REC file (output stream)




  void FindMinMaxMean(float ***aaafMask=NULL); 
                         // Read all density values from the tomogram and
                         // determine header.dmin,dmax,dmean
                         // The optional pmask argument allows you to indicate
                         // voxels who you want to exclude fron consideration.
  // We want different versions of the same tomogram to be directly comparable,
  // even if they were saved using different formats ("modes").
  // So we rescale them by dividing their densities values by the 
  // magnitude of the range of densities which are allowed in the file.
  // If an optional pmask argument is provided, then voxels in the mask
  // containing zeros 0 are ignored and are not rescaled.
  void Rescale01(float ***aaafMask=NULL); 

  // Sometimes we want to make the "white" voxels look black,
  // and the black voxels look white.
  // The Invert() will change the density of every voxel using:
  //   new_density = (ave_density - old_density)
  void Invert(float ***aaafMask=NULL); 

  // Print information about the tomogram size and format
  inline void PrintStats(ostream &out) { header.PrintStats(out); }

  // ------ Constructors, destructors, memory-allocation ------

  inline void Init() {
    afI = NULL;
    aaafI = NULL;
    header.mode = MrcHeader::MRC_MODE_FLOAT;
  }

  MrcSimple() {
    Init();
  }

  ~MrcSimple() {
    Dealloc();
  }


  void Resize(int const set_nvoxels[3]) {
    // Copy the contents of the image from an existing 3-D array
    Dealloc();
    for (int d=0; d < 3; d++) {
      header.nvoxels[d] = set_nvoxels[d];
      header.mvoxels[d] = set_nvoxels[d]; //default, can be overriden later
      header.cellA[d]   = set_nvoxels[d]; //default, can be overriden later
    }
    Alloc();
  }

  MrcSimple(int const set_nvoxels[3],
            float ***aaaf_contents);

  inline MrcSimple& operator = (const MrcSimple& source) {
    header = source.header;
    Dealloc(); // (just in case)
    Alloc();   // allocates and initializes afI and aaafI
    //for(Int iz=0; iz<header.nvoxels[2]; iz++)
    //  for(Int iy=0; iy<header.nvoxels[1]; iy++)
    //    for(Int ix=0; ix<header.nvoxels[0]; ix++)
    //      aaafI[iz][iy][ix] = source.aaafI[iz][iy][ix];
    // Use memcpy() instead:
    memcpy(afI, source.afI, 
           header.nvoxels[0]*header.nvoxels[1]*header.nvoxels[2]
           *sizeof(float));
  }

  inline MrcSimple(const MrcSimple& source) { //(not tested yet)
    Init();
    operator=(source);
  }

 private:

  // If you want to read the header and array separately, you can do that too:
  // using header.Read(mrc_file), and header.Write(mrc_file)
  // After that, you can read and write the rest of the file using:
 
  //Invoke only after header.Read()
  void ReadArray(istream& mrc_file, int const *axis_order);
  //Invoke after header.Write()
  void WriteArray(ostream& mrc_file) const; 

}; //MrcSimple




// Although not necessarily recommended, you can use << or >> to read/write
inline istream& operator >> (istream& mrc_file, 
                             MrcSimple& tomo)
{
  tomo.Read(mrc_file);
  return mrc_file;
}


inline ostream& operator << (ostream& mrc_file,
                             MrcSimple& tomo)
{
  tomo.Write(mrc_file);
  return mrc_file;
}




#endif //#ifndef _MRC_SIMPLE_H


