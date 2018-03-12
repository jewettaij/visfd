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
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;
#include <alloc3d.h>
#include <err_report.h>
#include "mrc_simple.h"


template<class RealNum >
inline RealNum ABS(RealNum x) { return ((x<0.0) ? -x: x); }


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
  header.mvoxels[0] = header.nvoxels[0];
  header.mvoxels[1] = header.nvoxels[1];
  header.mvoxels[2] = header.nvoxels[2];
  header.origin[0] = ceil(header.nvoxels[0]*0.5);  //default location of
  header.origin[1] = ceil(header.nvoxels[1]*0.5);  //origin is at the
  header.origin[2] = ceil(header.nvoxels[2]*0.5);  //center of the image
}



void MrcSimple::Dealloc() {
  Dealloc3D(header.nvoxels,
            &afI,
            &aaafI);
  header.origin[0] = -1;
  header.origin[1] = -1;
  header.origin[2] = -1;
}





void MrcSimple::Read(istream& mrc_file,
                     bool rescale,
		     float ***aaafMask) {
  header.Read(mrc_file);
  if ((header.mapCRS[0] != 1) ||
      (header.mapCRS[1] != 2) ||
      (header.mapCRS[2] != 3))
    throw InputErr("Error: This program currently only supports .MRC/.REC files in row-major\n"
      "       format, and whose mapC, mapR, mapS numbers (in the file header)\n"
      "       are 1,2,3, respectively.  Please use an alternate program to convert\n"
      "       to this format.  For more information on MRC file format, see:\n"
      "       http://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/image-processing-software/#image"
      "       http://bio3d.colorado.edu/imod/doc/mrc_format.txt\n"
                 "       http://ami.scripps.edu/software/mrctools/mrc_specification.php\n");
  ReadArray(mrc_file);
  if (rescale)
    Rescale01(aaafMask);
}



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
    header.use_signed_bytes = false;
  }
  Read(mrc_file, rescale, aaafMask);
  mrc_file.close();
}


void MrcSimple::ReadArray(istream& mrc_file) {
  Dealloc(); //free up any space you may have allocated earlier
  Alloc();   //allocate space for the array
  //mrc_file.read((char*)afI, sizeof(float)*num_voxels);
  for(Int iz=0; iz<header.nvoxels[2]; iz++) {
    //aaafI[iz] = new float* [header.nvoxels[1]];
    for(Int iy=0; iy<header.nvoxels[1]; iy++) {
      //aaafI[iz][iy] = &(afI[iz*header.nvoxels[0]*header.nvoxels[1] + 
      //                                    iy*header.nvoxels[0]]);
      for(Int ix=0; ix<header.nvoxels[0]; ix++) {

        switch (header.mode) {

        case MrcHeader::MRC_MODE_BYTE:
          if (header.use_signed_bytes)
          {
            int8_t entry;
            mrc_file.read((char*)&entry, sizeof(char));
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

        // Debugging check:
        //if ((abs(ix-100) <= 2) && (iy==100) && (iz==100))
        //  cerr << "aaafI["<<iz<<"]["<<iy<<"]["<<ix<<"] = "
        //       << aaafI[iz][iy][ix] << endl;

      } //for(Int ix=0; ix<header.nvoxels[0]; ix++) {
    } //for(Int iy=0; iy<header.nvoxels[1]; iy++) {
  } //for(Int iz=0; iz<header.nvoxels[2]; iz++) {
} //MrcSimple::ReadArray(istream& mrc_file)



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
  MrcHeader new_header = header;
  new_header.mode = MrcHeader::MRC_MODE_FLOAT;

  FindMinMaxMean(); // calculates the "dmin", "dmax", and "dmean" member values
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


void MrcSimple::Rescale01(float ***aaafMask)
{
  FindMinMaxMean(aaafMask); // find the min, max voxel intensities, considering
                            // only voxels where aaafMask[z][y][x] != 0
  float dmin = header.dmin;
  float dmax = header.dmax;
  for(Int iz=0; iz<header.nvoxels[2]; iz++) {
    for(Int iy=0; iy<header.nvoxels[1]; iy++) {
      for(Int ix=0; ix<header.nvoxels[0]; ix++) {
        aaafI[iz][iy][ix] = 
          (aaafI[iz][iy][ix] - dmin) / (dmax - dmin);
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


