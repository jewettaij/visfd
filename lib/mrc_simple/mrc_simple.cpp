// "MrcSimple": a class to read and write tomographic data.
//              (stored in .MRC and .REC files)
//
// Tomographic data is internally stored as arrays of floats, however
// this program can read MRC/REC files using several other numeric formats 
// as well (including signed and unsigned 8-bit and 16-bit integers).
// Note: As of 2018-6-30, this software does not attempt to detect 
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
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;
#include <alloc3d.hpp>
#include <err_report.hpp>
#include "mrc_simple.hpp"


template<class RealNum >
static inline RealNum ABS(RealNum x) { return ((x<0.0) ? -x: x); }


MrcSimple::
MrcSimple(int const set_nvoxels[3],
          float ***aaaf_contents)
{
  Init();
  Resize(set_nvoxels);
  for(Int iz=0; iz<header.nvoxels[2]; iz++)
    for(Int iy=0; iy<header.nvoxels[1]; iy++)
      for(Int ix=0; ix<header.nvoxels[0]; ix++)
        aaafI[iz][iy][ix] = aaaf_contents[iz][iy][ix];
  FindMinMaxMean();
}







void MrcSimple::Alloc() {
  Alloc3D(header.nvoxels,
          &afI,
          &aaafI);
}



void MrcSimple::Dealloc() {
  Dealloc3D(header.nvoxels,
            &afI,
            &aaafI);
}



template<class Int, class Entry>
void PermuteCArray(Int n,
		   Entry *target,
		   Int const *p) {
  assert(target && p);
  Int *tmp = new Int[n];
  for (Int i=0; i<n; i++)
    tmp[i] = target[i];
  for (Int i=0; i<n; i++)
    target[i] = tmp[ p[i] ];
  delete [] tmp;
}


template<class Int>
void PermuteInverseCArrayA(Int n,
                           Int *p) {
  assert(p);
  Int *tmp = new Int[n];
  for (Int i=0; i<n; i++)
    tmp[i] = p[i];
  for (Int i=0; i<n; i++)
    p[ tmp[i] ] = i;
  delete [] tmp;
}



void MrcSimple::Read(istream& mrc_file,
                     bool rescale,
                     float ***aaafMask) {
  header.Read(mrc_file);

  //    Row-major format ??
  Int *axis_order = NULL;
  if ((header.mapCRS[0] != 1) ||
      (header.mapCRS[1] != 2) ||
      (header.mapCRS[2] != 3)) {
    cerr <<
      "--------------------------------------------------------------------\n"
      "--------------------------------------------------------------------\n"
      "WARNING: The data in the file is not stored in row-major order.\n"
      "         This program has not been tested carefully on these kinds\n"
      "         of files.  Any calculated results will be written to a file\n"
      "         in row-major order, but the file's header entries may be\n"
      "         incorrect.  Use with caution.  Please report bugs.\n"
      "--------------------------------------------------------------------\n"
      "--------------------------------------------------------------------\n";
    // If the image data is not stored row-major format (header.mapCRS[]),
    // we must rearrange it. (This makes my once pretty code much uglier.)
    // Make a copy of "mapCRS[]" now so that later we can keep track of 
    // which axis (x,y,z) was stored in each index (i,j,k) of the array.
    axis_order = new Int[3];
    axis_order[0] = header.mapCRS[0] - 1;
    axis_order[1] = header.mapCRS[1] - 1;
    axis_order[2] = header.mapCRS[2] - 1;
    // Internally I decided to keep the image in row-major format,
    // so, to be consistent, I change the header entries accordingly.
    header.mapCRS[0] = 1;
    header.mapCRS[1] = 2;
    header.mapCRS[2] = 3;
    // Re-arrange the array storing the size of the image.
    PermuteCArray(3, header.nvoxels, axis_order);
    PermuteCArray(3, header.mvoxels, axis_order);
    PermuteCArray(3, header.origin, axis_order);
    // Do we need to re-arrange the other header entries as well?
    PermuteCArray(3, header.cellA, axis_order);
    // (Did I forget anything?)

  } // deal with non-row-major axis order


  // I'm not sure what some of these header entries mean, so I
  // just tried to fill them in with reasonable values:

  header.mvoxels[0] = header.nvoxels[0];
  header.mvoxels[1] = header.nvoxels[1];
  header.mvoxels[2] = header.nvoxels[2];

  // Read the image
  ReadArray(mrc_file, axis_order);


  if (rescale)
    Rescale01(aaafMask);


  // clean up
  if (axis_order)
    delete [] axis_order;

} //MrcSimple::Read()



void MrcSimple::Read(string in_file_name,
                     bool rescale,
                     float ***aaafMask) {
  Int len_in_file_name = in_file_name.size();
  fstream mrc_file;
  mrc_file.open(in_file_name.c_str(), ios::binary | ios::in);
  if (! mrc_file) 
    throw InputErr("Error: unable to open \""+ in_file_name +"\" for reading.\n");
  // Try to infer signed-vs-unsigned integers from the file name:
  //http://www.cgl.ucsf.edu/pipermail/chimera-users/2010-June/005245.html
  if ((len_in_file_name > 4)
      && 
      (in_file_name.substr(len_in_file_name-4, len_in_file_name) == ".rec")) {
    header.use_signed_bytes = false; //(note: this could be changed later by Read())
  }
  Read(mrc_file, rescale, aaafMask);
  mrc_file.close();
}





void MrcSimple::ReadArray(istream& mrc_file,
                          int const *axis_order=NULL) {

  Dealloc(); //free up any space you may have allocated earlier
  Alloc();   //allocate space for the array

  // Something like this might be faster:
  //mrc_file.read((char*)afI, sizeof(float)*num_voxels);
  // (I don't remember why I did not use this approach.
  //  Hopefully the code here fast enough in any case.)

  size_t num_voxels = header.nvoxels[0] *
                      header.nvoxels[1] *
                      header.nvoxels[2];

  Int iX = 0;
  Int iY = 0;
  Int iZ = 0;
  Int NX = header.nvoxels[0];
  Int NY = header.nvoxels[1];
  Int NZ = header.nvoxels[2];
  Int inv_axis_order[3];
  if (axis_order) {
    inv_axis_order[0] = axis_order[0];
    inv_axis_order[1] = axis_order[1];
    inv_axis_order[2] = axis_order[2];
    PermuteInverseCArrayA(3, inv_axis_order);
    NX = header.nvoxels[ inv_axis_order[0] ];
    NY = header.nvoxels[ inv_axis_order[1] ];
    NZ = header.nvoxels[ inv_axis_order[2] ];
  }
  

  // REMOVE THIS CRUFT:
  //for (size_t i = 0; i < num_voxels; i++) {

  for(Int iZ=0; iZ<NZ; iZ++) {
    for(Int iY=0; iY<NY; iY++) {
      for(Int iX=0; iX<NX; iX++) {
        Int ix = iX;
        Int iy = iY;
        Int iz = iZ;
        if (axis_order) {
          int ixyz[3];
          ixyz[0] = iX;
          ixyz[1] = iY;
          ixyz[2] = iZ;
          ix = ixyz[ inv_axis_order[0] ];
          iy = ixyz[ inv_axis_order[1] ];
          iz = ixyz[ inv_axis_order[2] ];
        }

        switch (header.mode) {

        case MrcHeader::MRC_MODE_BYTE:
          if (header.use_signed_bytes)
          {
            int8_t entry;
            mrc_file.read((char*)&entry, sizeof(char));

            #ifdef ENABLE_IMOD_COMPATIBILITY
            // MRC files using signed integers (eg. "mode 0") are
            // interpreted differently by IMOD than they are by other software.
            // In IMOD, voxel brightness values are adjusted to be so that 
            // they lie in the range from 0..255.
            // However files that use SIGNED bytes store their
            // voxel brightnesses in the range from -128..127.
            // For these files, IMOD automatically adds 128 to each voxel
            // brightness value to insure the result lies between 0..255.
            // Note: Neither UCSF Chimera, nor CCP-EM do this.
            //  (CCP-EM is the distributor of the "mrcfile" python module.)
            // To be compatible with the IMOD ecosystem, we must add 128:
            entry += 128;
            // Note: Later on when writing these files (in signed byte
            //       format), remember to subtract this offset before writing.
            // Note: Surprisingly, IMOD does not seem to add an offset to voxel
            //       brightnesses from files using "mode 1" (signed int16 format).
            #endif // #ifdef ENABLE_IMOD_COMPATIBILITY

            aaafI[iz][iy][ix] = static_cast<float>(entry);
          }
          else
          {
            uint8_t entry;
            mrc_file.read((char*)&entry, sizeof(char));
            aaafI[iz][iy][ix] = static_cast<float>(entry);
          }
          break;

        case MrcHeader::MRC_MODE_SHORT:
          {
            int16_t entry;
            mrc_file.read((char*)&entry, sizeof(int16_t));
            aaafI[iz][iy][ix] = static_cast<float>(entry);
          }
          break;

        case MrcHeader::MRC_MODE_USHORT:
          {
            uint16_t entry;
            mrc_file.read((char*)&entry, sizeof(uint16_t));
            aaafI[iz][iy][ix] = static_cast<float>(entry);
          }
          break;

        case MrcHeader::MRC_MODE_FLOAT:
          {
            float entry;
            mrc_file.read((char*)&entry, sizeof(float));
            //aaafI[iz][iy][ix] = static_cast<float>(entry);
            aaafI[iz][iy][ix] = entry;  //(already type float)
          }
          break;
        default:
          throw InputErr("UNSUPPORTED MODE in MRC file (unsupported MRC format)");
          exit(-1);
          break;
        } // switch (header.mode)


	// REMOVE THIS CRUFT:
        // find indices for the next entry in the array:
        //iX++;
        //if (iX >= NX) {
        //  iX = 0;
        //  iY++;
        //  if (iY >= NY) {
        //    iY = 0;
        //    iZ++;
        //    if (iZ >= NZ) {
        //      iZ = 0;
        //    }
        //  }
        //}


      } // for(Int iZ=0; iZ<NZ; iZ++)
    } // for(Int iY=0; iY<NY; iY++) {
  } // for(Int iX=0; iX<NX; iX++) {



  // REMOVE THIS CRUFT:
  //} //for (size_t i = 0; i < num_voxels; i++)



} //MrcSimple::ReadArray()




void MrcSimple::Write(string out_file_name) {
  fstream mrc_file;
  mrc_file.open(out_file_name.c_str(), ios::binary | ios::out);
  if (! mrc_file) 
    throw InputErr("Error: unable to open \""+ out_file_name+"\" for writing.\n");
  Write(mrc_file);  // You can also use "mrc_file << tomo;"
  mrc_file.close();
}


void MrcSimple::Write(ostream& mrc_file) {
  // First, write the MRC file header:
  //header.Write(mrc_file); <-- OOPS, that does not work.  Commenting out.
  // Note that any changes made to the tomogram density data will effect the
  // header. So we make a copy of the original header and modify it accordingly:
  // As of 2015-4-16, regardless of the original header, I save the list
  // of numbers in "float" format.  This means I must change mrc.mode.

  FindMinMaxMean(); // calculates the "dmin", "dmax", and "dmean" member values
  MrcHeader new_header = header;
  new_header.mode = MrcHeader::MRC_MODE_FLOAT;
  new_header.Write(mrc_file);

  // Finally, write the array of densities:
  WriteArray(mrc_file);

} //MrcSimple::Write(ostream& mrc_file) const






void MrcSimple::WriteArray(ostream& mrc_file) const {
  // Write the MRC contents:
  for(Int iz=0; iz<header.nvoxels[2]; iz++)
    for(Int iy=0; iy<header.nvoxels[1]; iy++)
      for(Int ix=0; ix<header.nvoxels[0]; ix++)
        mrc_file.write((char*)&(aaafI[iz][iy][ix]), sizeof(float));
}



void MrcSimple::FindMinMaxMean(float ***aaafMask) {
  double density_total = 0.0;
  double dmin = 0.0;  //impossible initial value
  double dmax = -1.0; //impossible initial value
  long long nvoxels_kept = 0;
  for(Int iz=0; iz<header.nvoxels[2]; iz++) {
    for(Int iy=0; iy<header.nvoxels[1]; iy++) {
      for(Int ix=0; ix<header.nvoxels[0]; ix++) {
        if (aaafMask && (aaafMask[iz][iy][ix] == 0))
          continue; //if a mask is supplied, ignore voxels when mask=0
        density_total += aaafI[iz][iy][ix];
        if (dmin > dmax) {
          dmin = aaafI[iz][iy][ix];
          dmax = aaafI[iz][iy][ix];
        }
        else {
          if (aaafI[iz][iy][ix] > dmax)
            dmax = aaafI[iz][iy][ix];
          if (aaafI[iz][iy][ix] < dmin)
            dmin = aaafI[iz][iy][ix];
        }
        nvoxels_kept++;
      }
    }
  }
  header.dmin = dmin;
  header.dmax = dmax;
  header.dmean = density_total  /  nvoxels_kept;
}


void MrcSimple::Rescale01(float ***aaafMask,
                          float outA,
                          float outB)
{
  FindMinMaxMean(aaafMask); // find the min, max voxel intensities, considering
                            // only voxels where aaafMask[z][y][x] != 0
  float dmin = header.dmin;
  float dmax = header.dmax;
  for(Int iz=0; iz<header.nvoxels[2]; iz++) {
    for(Int iy=0; iy<header.nvoxels[1]; iy++) {
      for(Int ix=0; ix<header.nvoxels[0]; ix++) {
        aaafI[iz][iy][ix] = outA +
          (outB-outA) * (aaafI[iz][iy][ix] - dmin) / (dmax - dmin);
      }
    }
  }
  FindMinMaxMean(NULL); // find min, max, mean, and update "header"
}



void MrcSimple::Invert(float ***aaafMask)
{
  double sum = 0.0;
  long n = 0;
  for(Int iz=0; iz<header.nvoxels[2]; iz++) {
    for(Int iy=0; iy<header.nvoxels[1]; iy++) {
      for(Int ix=0; ix<header.nvoxels[0]; ix++) {
        if (aaafMask && (aaafMask[iz][iy][ix] == 0))
          continue; //if a mask is supplied, ignore voxels when mask=0
        sum += aaafI[iz][iy][ix];
        n += 1;
      }
    }
  }
  double ave = sum / n;
  float dmin = ave;
  float dmax = ave;
  for(Int iz=0; iz<header.nvoxels[2]; iz++) {
    for(Int iy=0; iy<header.nvoxels[1]; iy++) {
      for(Int ix=0; ix<header.nvoxels[0]; ix++) {
        if (aaafMask && (aaafMask[iz][iy][ix] == 0))
          continue; //if a mask is supplied, ignore voxels when mask=0
        aaafI[iz][iy][ix] = 2.0*ave - aaafI[iz][iy][ix];
        if (aaafI[iz][iy][ix] < dmin)
          dmin = aaafI[iz][iy][ix];
        if (dmax < aaafI[iz][iy][ix])
          dmax = aaafI[iz][iy][ix];
      }
    }
  }
  header.dmean = ave;
  header.dmin = dmin;
  header.dmax = dmax;
}


