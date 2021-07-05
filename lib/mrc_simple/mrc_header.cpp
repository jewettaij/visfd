#include <cstring>
#include <cassert>
#include <fstream>
#include <iostream>
using namespace std;
#include "err_mrcfile.hpp"
#include "mrc_header.hpp"




void MrcHeader::Read(istream& mrc_file) {
  // Read the first "SIZE_HEADER" bytes into a temporary array "header_data"
  char header_data[MrcHeader::SIZE_HEADER]; //make sure we read the correct number of bytes
  mrc_file.read(header_data, MrcHeader::SIZE_HEADER);

  // copy the bytes we dont use into the remaining_raw_data array
  memcpy(remaining_raw_data,  
         header_data + MrcHeader::SIZE_HEADER_USED, 
         MrcHeader::SIZE_REMAINING_HEADER);
          
  // Convert the binary data into numbers (Ints and floats)

  // 0    NX       number of columns (fastest changing in map)
  // 1    NY       number of rows
  // 2    NZ       number of sections (slowest changing in map)

  // Commenting out the lines with "mrc_file.read()", although they also work.
  //mrc_file.read((char*)&nvoxels[0], sizeof(Int));
  //mrc_file.read((char*)&nvoxels[1], sizeof(Int));
  //mrc_file.read((char*)&nvoxels[2], sizeof(Int));
  //nvoxels[0] = *((Int*)(pheader_data)+0);
  nvoxels[0] = *(reinterpret_cast<Int*>(header_data)+0);
  nvoxels[1] = *(reinterpret_cast<Int*>(header_data)+1);
  nvoxels[2] = *(reinterpret_cast<Int*>(header_data)+2);


  // 3    MODE     data type :
  //       0        image : signed 8-bit bytes range -128 to 127
  //       1        image : 16-bit halfwords
  //       2        image : 32-bit reals
  //       3        transform : complex 16-bit integers
  //       4        transform : complex 32-bit reals
  //       6        image : unsigned 16-bit range 0 to 65535
  //mrc_file.read((char*)&mode, sizeof(Int));
  mode = *(reinterpret_cast<Int*>(header_data)+3);


  if (mode == MrcHeader::MRC_MODE_BYTE) {
    // Try to determine if 8-bit integers are signed or unsigned.
    // Apparently, there are no general guidelines indicating whether 8-bit 
    // integers are signed or unsigned in MRC files.  
    // If either of these statements are true, bytes are signed.
    //   1) check the filename to see if it ends in ".rec" instead of ".mrc"
    //     http://www.cgl.ucsf.edu/pipermail/chimera-users/2010-June/005245.html
    //   2) check for a bit in the "imodFlags" section
    //      of the header.  (That's what we do here.)
    //     http://bio3d.colorado.edu/imod/doc/mrc_format.txt
    // Otherwise, the are (by default) assumbe to be unsigned.
    //
    //152 230  4   int     imodStamp   1146047817 indicates that file was created by IMOD or 
    //                                 other software that uses bit flags in the following field
    //156 234  4   int     imodFlags   Bit flags:
    //                                 1 = bytes are stored as signed
    //                                 2 = pixel spacing was set from size in extended header 
    //                                 4 = origin is stored with sign inverted from definition
    //                                     below
    Int imodStamp = *(reinterpret_cast<Int*>(header_data)+38); //(byte 152-155)
    if (imodStamp == 1146047817) { 
      // Then the file was created by IMOD or other software using "bit flags"
      // See: http://bio3d.colorado.edu/imod/doc/mrc_format.txt
      Int imodFlags = *(reinterpret_cast<Int*>(header_data)+39);//(byte 156-159)
      use_signed_bytes = (imodFlags & 1);
    }
  }



  // 4    NXSTART number of first column in map (Default = 0)
  // 5    NYSTART number of first row in map
  // 6    NZSTART number of first section in map
  //mrc_file.read((char*)&nstart[0], sizeof(Int));
  //mrc_file.read((char*)&nstart[1], sizeof(Int));
  //mrc_file.read((char*)&nstart[2], sizeof(Int));
  nstart[0] = *(reinterpret_cast<Int*>(header_data)+4);
  nstart[1] = *(reinterpret_cast<Int*>(header_data)+5);
  nstart[2] = *(reinterpret_cast<Int*>(header_data)+6);


  // 7     MX       number of intervals along X
  // 8     MY       number of intervals along Y
  // 9     MZ       number of intervals along Z
  //mrc_file.read((char*)&mvoxels[0], sizeof(Int));
  //mrc_file.read((char*)&mvoxels[1], sizeof(Int));
  //mrc_file.read((char*)&mvoxels[2], sizeof(Int));
  mvoxels[0] = *(reinterpret_cast<Int*>(header_data)+7);
  mvoxels[1] = *(reinterpret_cast<Int*>(header_data)+8);
  mvoxels[2] = *(reinterpret_cast<Int*>(header_data)+9);


  // 10-12    CELLA    cell dimensions in angstroms
  //mrc_file.read((char*)&(cellA[0]), sizeof(float));
  //mrc_file.read((char*)&(cellA[1]), sizeof(float));
  //mrc_file.read((char*)&(cellA[2]), sizeof(float));
  cellA[0] = *(reinterpret_cast<float*>(header_data)+10);
  cellA[1] = *(reinterpret_cast<float*>(header_data)+11);
  cellA[2] = *(reinterpret_cast<float*>(header_data)+12);


  // 13-15    CELLB    cell angles in degrees
  //mrc_file.read((char*)&(cellB[0]), sizeof(float));
  //mrc_file.read((char*)&(cellB[1]), sizeof(float));
  //mrc_file.read((char*)&(cellB[2]), sizeof(float));
  cellB[0] = *(reinterpret_cast<float*>(header_data)+13);
  cellB[1] = *(reinterpret_cast<float*>(header_data)+14);
  cellB[2] = *(reinterpret_cast<float*>(header_data)+15);


  // 16    MAPC     axis corresp to cols (1,2,3 for X,Y,Z)
  // 17    MAPR     axis corresp to rows (1,2,3 for X,Y,Z)
  // 18    MAPS     axis corresp to sections (1,2,3 for X,Y,Z)
  //mrc_file.read((char*)&(mapCRS[0]), sizeof(Int));
  //mrc_file.read((char*)&(mapCRS[1]), sizeof(Int));
  //mrc_file.read((char*)&(mapCRS[2]), sizeof(Int));
  mapCRS[0] = *(reinterpret_cast<Int*>(header_data)+16);
  mapCRS[1] = *(reinterpret_cast<Int*>(header_data)+17);
  mapCRS[2] = *(reinterpret_cast<Int*>(header_data)+18);


  // 19    DMIN     minimum brightness value
  // 20    DMAX     maximum brightness value
  // 21    DMEAN    mean brightness value
  //mrc_file.read((char*)&(dmin), sizeof(float));
  //mrc_file.read((char*)&(dmax), sizeof(float));
  //mrc_file.read((char*)&(dmean), sizeof(float));
  dmin  = *(reinterpret_cast<float*>(header_data)+19);
  dmax  = *(reinterpret_cast<float*>(header_data)+20);
  dmean = *(reinterpret_cast<float*>(header_data)+21);


  #ifdef ENABLE_IMOD_COMPATIBILITY
  if (use_signed_bytes) {
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
    dmin += 128;
    dmax += 128;
    dmean += 128;
    // Note: Later on when writing these files (in signed byte
    //       format), remember to subtract this offset before writing.
    // Note: Surprisingly, IMOD does not seem to add an offset to voxel
    //       brightnesses from files using "mode 1" (signed int16 format).
  }
  #endif // #ifdef ENABLE_IMOD_COMPATIBILITY


  // 22    ISPG     space group number 0 or 1 (default=0)
  //mrc_file.read((char*)&(ispg), sizeof(Int));
  ispg = *(reinterpret_cast<Int*>(header_data)+22);

  // 23    NSYMBT   number of bytes used for symmetry data (0 or 80)
  //mrc_file.read((char*)&(nsymbt), sizeof(Int));
  nsymbt = *(reinterpret_cast<Int*>(header_data)+23);

  // 24-49    EXTRA    extra space used for anything   - 0 by default
  //mrc_file.read(extra_raw_data, 25*MrcHeader::SIZE_PER_FIELD)
  // copy the bytes we dont use into the remaining_raw_data array
  memcpy(extra_raw_data,
         header_data + (24*MrcHeader::SIZE_PER_FIELD),
         25*MrcHeader::SIZE_PER_FIELD);

  // 50-52    ORIGIN   origin in X,Y,Z used for transforms
  //mrc_file.read((char*)&(origin[0]), sizeof(float));
  //mrc_file.read((char*)&(origin[1]), sizeof(float));
  //mrc_file.read((char*)&(origin[2]), sizeof(float));
  origin[0] = *(reinterpret_cast<float*>(header_data)+49);
  origin[1] = *(reinterpret_cast<float*>(header_data)+50);
  origin[2] = *(reinterpret_cast<float*>(header_data)+51);

} //MrcHeader::Read(istream& mrc_file)




void MrcHeader::Write(ostream& mrc_file) const {

  // 0    NX       number of columns (fastest changing in map)
  // 1    NY       number of rows
  // 2    NZ       number of sections (slowest changing in map)

  // Commenting out the lines with "mrc_file.read()", although they also work.
  mrc_file.write((char*)&nvoxels[0], sizeof(Int));
  mrc_file.write((char*)&nvoxels[1], sizeof(Int));
  mrc_file.write((char*)&nvoxels[2], sizeof(Int));

  // 3    MODE     data type :
  //       0        image : signed 8-bit bytes range -128 to 127
  //       1        image : 16-bit halfwords
  //       2        image : 32-bit reals
  //       3        transform : complex 16-bit integers
  //       4        transform : complex 32-bit reals
  //       6        image : unsigned 16-bit range 0 to 65535
  mrc_file.write((char*)&mode, sizeof(Int));


  // 4    NXSTART number of first column in map (Default = 0)
  // 5    NYSTART number of first row in map
  // 6    NZSTART number of first section in map
  mrc_file.write((char*)&nstart[0], sizeof(Int));
  mrc_file.write((char*)&nstart[1], sizeof(Int));
  mrc_file.write((char*)&nstart[2], sizeof(Int));


  // 7    MX       number of intervals along X
  // 8    MY       number of intervals along Y
  // 9    MZ       number of intervals along Z
  mrc_file.write((char*)&mvoxels[0], sizeof(Int));
  mrc_file.write((char*)&mvoxels[1], sizeof(Int));
  mrc_file.write((char*)&mvoxels[2], sizeof(Int));


  // 10-12    CELLA    cell dimensions in angstroms
  mrc_file.write((char*)&(cellA[0]), sizeof(float));
  mrc_file.write((char*)&(cellA[1]), sizeof(float));
  mrc_file.write((char*)&(cellA[2]), sizeof(float));


  // 13-15    CELLB    cell angles in degrees
  mrc_file.write((char*)&(cellB[0]), sizeof(float));
  mrc_file.write((char*)&(cellB[1]), sizeof(float));
  mrc_file.write((char*)&(cellB[2]), sizeof(float));


  // 16    MAPC     axis corresp to cols (1,2,3 for X,Y,Z)
  // 17    MAPR     axis corresp to rows (1,2,3 for X,Y,Z)
  // 18    MAPS     axis corresp to sections (1,2,3 for X,Y,Z)
  mrc_file.write((char*)&(mapCRS[0]), sizeof(Int));
  mrc_file.write((char*)&(mapCRS[1]), sizeof(Int));
  mrc_file.write((char*)&(mapCRS[2]), sizeof(Int));

  // 19    DMIN     minimum brightness value
  // 20    DMAX     maximum brightness value
  // 21    DMEAN    mean brightness value
  mrc_file.write((char*)&(dmin), sizeof(float));
  mrc_file.write((char*)&(dmax), sizeof(float));
  mrc_file.write((char*)&(dmean), sizeof(float));


  // 22    ISPG     space group number 0 or 1 (default=0)
  mrc_file.write((char*)&(ispg), sizeof(Int));

  // 23    NSYMBT   number of bytes used for symmetry data (0 or 80)
  mrc_file.write((char*)&(nsymbt), sizeof(Int));

  // 24-49    EXTRA    extra space used for anything   - 0 by default
  mrc_file.write(extra_raw_data,
                 25*sizeof(Int));

  // 50-52    ORIGIN   origin in X,Y,Z used for transforms
  mrc_file.write((char*)&(origin[0]), sizeof(float));
  mrc_file.write((char*)&(origin[1]), sizeof(float));
  mrc_file.write((char*)&(origin[2]), sizeof(float));

  // I don't care about the remaining junk in the header
  mrc_file.write(remaining_raw_data,
                 MrcHeader::SIZE_REMAINING_HEADER);

} //MrcHeader::Write(ostream& mrc_file) const





void MrcHeader::PrintStats(ostream& out) {
  out << "  mrc file stats:\n"
      << "    number of voxels: "
      << nvoxels[0] << " x " 
      << nvoxels[1] << " x "
      << nvoxels[2] << "" 
      << endl;
  //out << "  mrc file physical dimensions ("<< cellA[0] << "," << cellA[1] << "," << cellA[2] << ")" << endl;
  out << "    voxel size in file header: "
      << cellA[0]/nvoxels[0] << " x " 
      << cellA[1]/nvoxels[1] << " x " 
      << cellA[2]/nvoxels[2] << "" << endl;
  out << "    table axis order: "
      << mapCRS[0] << " " 
      << mapCRS[1] << " "
      << mapCRS[2] << ""
      << endl;
  out << "    mode: " << mode << endl;
  out << "    minimum brightness: "<< dmin << endl;
  out << "    maximum brightness: "<< dmax << endl;
  out << "    mean brightness: "<< dmean << endl;
  out << "    origin: "
      << origin[0] << " " 
      << origin[1] << " " 
      << origin[2] << "" << endl;
}
