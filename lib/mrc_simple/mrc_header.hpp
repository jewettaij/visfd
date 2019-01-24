#ifndef MRC_HEADER_HPP
#define MRC_HEADER_HPP

#include <cstdint>
using namespace std;



// Documentation on the MRC file format is available here:
//http://bio3d.colorado.edu/imod/doc/mrc_format.txt
//http://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/image-processing-software/#image
//http://ami.scripps.edu/software/mrctools/mrc_specification.php
//http://bio3d.colorado.edu/imod/betaDoc/libhelp/mrcfiles.html#TOP




//typedef int Int;
typedef int32_t Int;  //An integer type gauranteed to be exactly 32 bits wide.
                      //(All integers in the header are 32 bits wide)
                      //(As of 2015-3, if you use gcc/g++, you must compile 
                      // using the -std=c++11 flag to support "int32_t")






class MrcHeader {

public:

  Int nvoxels[3];  // How many voxels are in the tomogram? (in x,y,z direction)
                   // (total number = nvoxels[0]*nvoxels[1]*nvoxels[2])

  Int mode; // Format of the numbers representing the density at each voxel 
            // location in the array (as stored in an MRC/REC file):
  const static Int  MRC_MODE_UNDEFINED     = -1;
  const static Int  MRC_MODE_BYTE          = 0;
  const static Int  MRC_MODE_SHORT         = 1;
  const static Int  MRC_MODE_FLOAT         = 2;
  const static Int  MRC_MODE_COMPLEX_SHORT = 3;  // (IGNORED as of 2015-4-16)
  const static Int  MRC_MODE_COMPLEX_FLOAT = 4;  // (IGNORED as of 2015-4-16)
  const static Int  MRC_MODE_USHORT        = 6;
  const static Int  MRC_MODE_RGB           = 16; // (IGNORED as of 2015-4-16)
            // Note: This is only useful when reading or writing MRC files.
            // (Internally, densities are represented using arrays of floats.)


  // ----- (The "nstart" and "mvoxels" entries are not useful for us, but -----
  // -----  they appear early in the header file so I include them here.) -----
  Int nstart[3];   // starting-point for sub-images (IGNORED)             -----
  Int mvoxels[3];  // not sure what this is (IGNORED)                     -----
  // --------------------------------------------------------------------------


  float cellA[3]; // physical size (in Angstroms) of the rectangular box 
                  // containing the tomogram.  Useful to calculate voxel size:
                  // (voxel-width in x direction = cellA[0]/nvoxels[0]
                  //  voxel-width in y direction = cellA[1]/nvoxels[1]
                  //  voxel-width in z direction = cellA[2]/nvoxels[2])

  void Read(istream& mrc_file); //Read an .MRC/.REC file header
  void Write(ostream& mrc_file) const; //Write an .MRC/.REC file header





  // ----------------------------------------------------------
  // ------------------- Ignore the rest ----------------------
  // ----------------------------------------------------------





  // (The remaining entries in a header file are not currently useful to me.)

  float cellB[3];  // unit cell angles (IGNORED, assume 90 degrees)

  Int mapCRS[3];        // which axis (x,y,z) corresponds to which index (i,j,k)
                        // in the array of voxel intensities (aaaI[k][j][i]) ? 
                        // The default (1,2,3) corresponds to x<->i, y<->j, z<->k

  float dmin;           // Minimum-possible numeric value of entry in the array
  float dmax;           // Maximum-possible numeric value of entry in the array
  float dmean;          // Average value of numbers in the array
  Int ispg;             // no idea what this is (IGNORED)
  Int nsymbt;           // no idea what this is (IGNORED)

  static const Int SIZE_HEADER = 1024;
  static const Int SIZE_PER_FIELD = 4;  //Each Int or float requires 4 bytes
  static const Int NUM_USED_FIELDS = 52; //Number of entries read from header
  static const Int SIZE_HEADER_USED = NUM_USED_FIELDS*SIZE_PER_FIELD;//in bytes
  static const Int SIZE_REMAINING_HEADER = SIZE_HEADER - SIZE_HEADER_USED;

  char extra_raw_data[25*SIZE_PER_FIELD]; // (IGNORED)

  float origin[3];      // I'm guessing this is only relevant to reconstructions
                        // I'm guessing this might be the location of a point
                        // within the pivot axis (axes?) around which the sample
                        // was rotated when generating the original tilt series.

  char remaining_raw_data[SIZE_REMAINING_HEADER];

  bool use_signed_bytes;// signed bytes? (only relevant when mode=MRC_MODE_BYTE)
  void PrintStats(ostream &out);//Print information about tomogram size & format


  // constructor sets the tomogram size to 0,0,0
  MrcHeader() {
    mapCRS[0] = 1;
    mapCRS[1] = 2;
    mapCRS[2] = 3;
    nvoxels[0] = 0; 
    nvoxels[1] = 0;
    nvoxels[2] = 0;
    cellA[0] = 0.0; 
    cellA[1] = 0.0; 
    cellA[2] = 0.0;
    dmin = 0.0;   //impossible values
    dmax = -1.0;  //impossible values
    dmean = -1.0;
    mode = MRC_MODE_UNDEFINED;
    cellB[0] = 90.0; 
    cellB[1] = 90.0; 
    cellB[2] = 90.0; 
    mvoxels[0] = 0; 
    mvoxels[1] = 0;
    mvoxels[2] = 0;
    origin[0] = -1;
    origin[1] = -1;
    origin[2] = -1;
    nstart[0] = -1;
    nstart[1] = -1;
    nstart[2] = -1;
    use_signed_bytes = true;
  }

};


// Although not necessarily recommended, you can use << or >> to read/write
static inline
istream& operator>>(istream& mrc_file, 
                    MrcHeader& mrc_header)
{
  mrc_header.Read(mrc_file);
  return mrc_file;
}


static inline
ostream& operator<<(ostream& mrc_file, 
                    const MrcHeader& mrc_header)
{
  mrc_header.Write(mrc_file);
  return mrc_file;
}




#endif //#ifndef MRC_HEADER_HPP


