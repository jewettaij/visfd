#ifndef _MRC_SIMPLE_HPP
#define _MRC_SIMPLE_HPP

#include <cstring>
#include <iostream>
using namespace std;
#include "mrc_header.hpp"


/// @brief "MrcSimple": a class to read and write tomographic data.
///              (stored in .MRC and .REC files)
///
/// Tomographic data is internally stored as arrays of floats, however
/// this program can read MRC/REC files using several other numeric formats 
/// as well (including signed and unsigned 8-bit and 16-bit integers).
///
/// --- PLEASE REPORT BUGS ---
///
/// @note As of 2019-4-08, support for MRC files which are not in row-major
///       format has been implemented but has not been tested.
///       http://en.wikipedia.org/wiki/Row-major_order
///       (IE This program has only been tested on MRC files whose "mapC",
///        "mapR", and "mapS" (mapCRS[]) header data equals 1,2,3 respectively.)
///        
/// @note As of 2019-4-08, the interpretation of signed bytes (mode 0) MRC files
///       is different in IMOD, compared to other software tools, including this
///       software.  When reading files using signed bytes, IMOD adds 128 to
///       all the entries (so that the resulting numbers are strictly positive).
///       This does not occur when reading signed-byte files using this library.
///       (This software may also fail to realize when signed bytes are in use.)
///        
/// @note As of 2019-4-08, this software does not attempt to detect 
///       or convert big-endian or little-endian numbers.
///       To be on the safe side, compile this software on hardware which is
///       compatible with the hardware which was used to write the MRC file.
///       (IE intel CPUs, or ARM CPUs using the same endian settings.)
///
/// @note  Documentation on the MRC file format is available here:
/// http://bio3d.colorado.edu/imod/doc/mrc_format.txt
/// http://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/image-processing-software/#image
/// http://ami.scripps.edu/software/mrctools/mrc_specification.php
/// http://bio3d.colorado.edu/imod/betaDoc/libhelp/mrcfiles.html#TOP
///
/// @note  Most of this code was pillaged from IMOD
///        (by David Mastronarde and coworkers)


class MrcSimple {
public:
  MrcHeader header;     ///< contains the size and dimensions of the tomogram:
                        ///<   (header.nvoxels[] stores the number of voxels
                        ///<                         in the x,y,z directions.
                        ///<    header.cellA[] stores the physical size of
                        ///<                    the tomographic reconstruction)

  float *afI;           ///< pointer to a (contiguous) 1-D array containing 
                        ///< the entire tomogram. The data at x,y,z is located
                        ///<   at  afI[ix + iy*sizez + iz*sizex*sizey]

  float ***aaafI;       ///< Useful if you prefer to use 3-index notation:
                        ///<    aafI[iz][iy][ix] =
                        ///<     afI[ix + iy*sizez + iz*sizex*sizey]
private:

  /// @brief  Allocate memory for the tomogram (image).
  ///   Allocates and initializes afI and aaafI. Invoke
  ///   after using ReadHeader(), or after modifying the header.
  void Alloc();

  /// @brief  Deallocate memory for the tomogram (image).
  void Dealloc();

public:

  /// @brief  Read an .MRC/.REC file
  void Read(string mrc_file_name,  //!<name of the file
            bool rescale=false,     //!<Optional: rescale densities from 0 to 1?
            float ***aaafMask=NULL//!<Optional: ignore zero-valued voxels in the aaafMask[][][] array
            );

  /// @brief  Write an .MRC/.REC file
  void Write(string mrc_file_name); //Write an .MRC/.REC file

  /// @brief  Read an .MRC/.REC file
  void Read(istream& mrc_file,      //!< Read an .MRC/.REC file (input stream)
            bool rescale=false,      //!<Optional: rescale densities from 0 to 1?
	    float ***aaafMask=NULL  //!<Optional: ignore zero-valued voxels in the aaafMask[][][] array
            );
  void Write(ostream& mrc_file); //Write an .MRC/.REC file (output stream)




  // @brief Read all density values from the tomogram and
  //  determine header.dmin,dmax,dmean
  //  The optional pmask argument allows you to indicate
  //  voxels who you want to exclude fron consideration.
  //  We want different versions of the same tomogram to be directly comparable,
  //  even if they were saved using different formats ("modes").
  //  So we rescale them by dividing their densities values by the 
  //  magnitude of the range of densities which are allowed in the file.
  //  If an optional pmask argument is provided, then voxels in the mask
  //  containing zeros 0 are ignored and are not rescaled.
  void FindMinMaxMean(float ***aaafMask=NULL); 

  /// @brief  linearly rescale the voxel intensities in the image so that the
  ///         new voxel intensities range from outA to outB
  void Rescale01(float ***aaafMask,
                 float outA=0.0,
                 float outB=1.0);

  /// @brief Sometimes we want to make the "white" voxels look black,
  /// and the black voxels look white.
  /// The Invert() will change the density of every voxel using:
  ///   new_density = (ave_density - old_density)
  void Invert(float ***aaafMask=NULL); 

  /// @brief  Print information about the tomogram size and format
  void PrintStats(ostream &out) { header.PrintStats(out); }

  // ------ Constructors, destructors, memory-allocation ------

private:
  void Init() {
    afI = NULL;
    aaafI = NULL;
    header.mode = MrcHeader::MRC_MODE_FLOAT;
  }

public:
  MrcSimple() {
    Init();
  }

  ~MrcSimple() {
    Dealloc();
  }


  /// @brief  Change the size (number of voxels) in the image.

  void Resize(int const set_nvoxels[3]//!<number of voxels in the xyz directions
              )
  {
    // Copy the contents of the image from an existing 3-D array
    Dealloc();
    for (int d=0; d < 3; d++) {
      header.nvoxels[d] = set_nvoxels[d];
      header.mvoxels[d] = set_nvoxels[d]; //default, can be overriden later
      //header.cellA[d]   = set_nvoxels[d]; //default, can be overriden later
      //header.origin[d]  = 0.0;
    }
    Alloc();
  }


  MrcSimple(int const set_nvoxels[3],//!< number of voxels in the xyz directions
            float ***aaaf_contents //!< a 3D array of voxel values
            );


  MrcSimple(const MrcSimple& source):
    MrcSimple(source.header.nvoxels,
              source.aaafI)
  {
    header = source.header;
  }

  void swap(MrcSimple &other) {
    std::swap(header, other.header);
    std::swap(afI, other.afI);
    std::swap(aaafI, other.aaafI);
    // (do I need to do something fancier for multidimensional arrays (aaafI)?)
  }

  MrcSimple&
    operator = (MrcSimple source) {
    this->swap(source);
    return *this;
  }


private:

  /// @brief
  /// Invoke only after header.Read()
  /// If you want to read the header and array separately, you can do that too:
  /// using header.Read(mrc_file)
  /// After that, you can read the rest of the file using:
  void ReadArray(istream& mrc_file, int const *axis_order);

  /// @brief
  /// Invoke only after header.Read()
  /// If you want to write the header and array separately, you can do that too:
  /// using header.Write(mrc_file)
  /// After that, you can write the rest of the file using:
  void WriteArray(ostream& mrc_file) const; 

}; //MrcSimple




/// @brief you can also use >> to write the file
static inline
istream& operator >> (istream& mrc_file, 
                      MrcSimple& tomo)
{
  tomo.Read(mrc_file);
  return mrc_file;
}


/// @brief you can also use << to read the file
static inline
ostream& operator << (ostream& mrc_file,
                      MrcSimple& tomo)
{
  tomo.Write(mrc_file);
  return mrc_file;
}



/// @brief  This simple function prints a warning message
///         whenever the image uses signed bytes.
///         (Since my users are frequently IMOD users,
///          I find myself using this silly function frequently.)

inline void
WarnMRCSignedBytes(const MrcSimple &image,
                   string file_name,
                   ostream &report_warning)
{
  if (image.header.use_signed_bytes &&
      image.header.mode == MrcHeader::MRC_MODE_BYTE) {
    report_warning <<
      "#####################################################################\n"
      "WARNING: File \""<< file_name <<"\"\n"
      "         uses signed bytes.\n"
      "         Do not compare voxel brightnesses in this file with voxel\n"
      "         brightnesses reported by IMOD or 3dmod.  (To avoid this message\n"
      "         in the future, convert this file to an unambiguous format.\n"
      "         You can do this using \"convert_to_float file1 file2\",\n"
      "         or \"newstack -mode 2 -in file1 -ou file2\", for example.)\n"
      "#####################################################################"
                   << endl;
  }
}


#endif //#ifndef _MRC_SIMPLE_HPP


