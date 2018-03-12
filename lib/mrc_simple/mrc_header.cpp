#include <cstring>
#include <cassert>
#include <fstream>
#include <iostream>
using namespace std;
#include <err_report.h>
#include "mrc_header.h"




void MrcHeader::Read(istream& mrc_file) {
  // Read the first "SIZE_HEADER" bytes into a temporary array "header_data"
  char header_data[MrcHeader::SIZE_HEADER]; //make sure we read the correct number of bytes
  mrc_file.read(header_data, MrcHeader::SIZE_HEADER);

  // copy the bytes we dont use into the remaining_raw_data array
  memcpy(remaining_raw_data,  
         header_data + MrcHeader::SIZE_HEADER_USED, 
         MrcHeader::SIZE_REMAINING_HEADER);
          
  // Convert the binary data into numbers (Ints and floats)

  // 1    NX       number of columns (fastest changing in map)
  // 2    NY       number of rows   
  // 3    NZ       number of sections (slowest changing in map)

  // Commenting out the lines with "mrc_file.read()", although they also work.
  //mrc_file.read((char*)&nvoxels[0], sizeof(Int));
  //mrc_file.read((char*)&nvoxels[1], sizeof(Int));
  //mrc_file.read((char*)&nvoxels[2], sizeof(Int));
  //nvoxels[0] = *((Int*)(pheader_data)+0);
  nvoxels[0] = *(reinterpret_cast<Int*>(header_data)+0);
  nvoxels[1] = *(reinterpret_cast<Int*>(header_data)+1);
  nvoxels[2] = *(reinterpret_cast<Int*>(header_data)+2);


  // 4    MODE     data type :
  //       0        image : signed 8-bit bytes range -128 to 127
  //       1        image : 16-bit halfwords
  //       2        image : 32-bit reals
  //       3        transform : complex 16-bit integers
  //       4        transform : complex 32-bit reals
  //       6        image : unsigned 16-bit range 0 to 65535
  //mrc_file.read((char*)&mode, sizeof(Int));
  mode = *(reinterpret_cast<Int*>(header_data)+3);


  // 5    NXSTART number of first column in map (Default = 0)
  // 6    NYSTART number of first row in map
  // 7    NZSTART number of first section in map
  //mrc_file.read((char*)&nstart[0], sizeof(Int));
  //mrc_file.read((char*)&nstart[1], sizeof(Int));
  //mrc_file.read((char*)&nstart[2], sizeof(Int));
  nstart[0] = *(reinterpret_cast<Int*>(header_data)+4);
  nstart[1] = *(reinterpret_cast<Int*>(header_data)+5);
  nstart[2] = *(reinterpret_cast<Int*>(header_data)+6);


  // 8     MX       number of intervals along X
  // 9     MY       number of intervals along Y
  // 10    MZ       number of intervals along Z
  //mrc_file.read((char*)&mvoxels[0], sizeof(Int));
  //mrc_file.read((char*)&mvoxels[1], sizeof(Int));
  //mrc_file.read((char*)&mvoxels[2], sizeof(Int));
  mvoxels[0] = *(reinterpret_cast<Int*>(header_data)+7);
  mvoxels[1] = *(reinterpret_cast<Int*>(header_data)+8);
  mvoxels[2] = *(reinterpret_cast<Int*>(header_data)+9);


  // 11-13    CELLA    cell dimensions in angstroms
  //mrc_file.read((char*)&(cellA[0]), sizeof(float));
  //mrc_file.read((char*)&(cellA[1]), sizeof(float));
  //mrc_file.read((char*)&(cellA[2]), sizeof(float));
  cellA[0] = *(reinterpret_cast<float*>(header_data)+10);
  cellA[1] = *(reinterpret_cast<float*>(header_data)+11);
  cellA[2] = *(reinterpret_cast<float*>(header_data)+12);


  // 14-16    CELLB    cell angles in degrees
  //mrc_file.read((char*)&(cellB[0]), sizeof(float));
  //mrc_file.read((char*)&(cellB[1]), sizeof(float));
  //mrc_file.read((char*)&(cellB[2]), sizeof(float));
  cellB[0] = *(reinterpret_cast<float*>(header_data)+13);
  cellB[1] = *(reinterpret_cast<float*>(header_data)+14);
  cellB[2] = *(reinterpret_cast<float*>(header_data)+15);


  // 17    MAPC     axis corresp to cols (1,2,3 for X,Y,Z)
  // 18    MAPR     axis corresp to rows (1,2,3 for X,Y,Z)
  // 19    MAPS     axis corresp to sections (1,2,3 for X,Y,Z)
  //mrc_file.read((char*)&(mapCRS[0]), sizeof(Int));
  //mrc_file.read((char*)&(mapCRS[1]), sizeof(Int));
  //mrc_file.read((char*)&(mapCRS[2]), sizeof(Int));
  mapCRS[0] = *(reinterpret_cast<Int*>(header_data)+16);
  mapCRS[1] = *(reinterpret_cast<Int*>(header_data)+17);
  mapCRS[2] = *(reinterpret_cast<Int*>(header_data)+18);


  // 20    DMIN     minimum density value
  // 21    DMAX     maximum density value
  // 22    DMEAN    mean density value
  //mrc_file.read((char*)&(dmin), sizeof(float));
  //mrc_file.read((char*)&(dmax), sizeof(float));
  //mrc_file.read((char*)&(dmean), sizeof(float));
  dmin  = *(reinterpret_cast<float*>(header_data)+19);
  dmax  = *(reinterpret_cast<float*>(header_data)+20);
  dmean = *(reinterpret_cast<float*>(header_data)+21);


  // 23    ISPG     space group number 0 or 1 (default=0)
  //mrc_file.read((char*)&(ispg), sizeof(Int));
  ispg = *(reinterpret_cast<Int*>(header_data)+22);

  // 24    NSYMBT   number of bytes used for symmetry data (0 or 80)
  //mrc_file.read((char*)&(nsymbt), sizeof(Int));
  nsymbt = *(reinterpret_cast<Int*>(header_data)+23);

  // 25-49    EXTRA    extra space used for anything   - 0 by default
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

} //MrcHeader::Read(istream& mrc_file)




void MrcHeader::Write(ostream& mrc_file) const {

  // 1    NX       number of columns (fastest changing in map)
  // 2    NY       number of rows   
  // 3    NZ       number of sections (slowest changing in map)

  // Commenting out the lines with "mrc_file.read()", although they also work.
  mrc_file.write((char*)&nvoxels[0], sizeof(Int));
  mrc_file.write((char*)&nvoxels[1], sizeof(Int));
  mrc_file.write((char*)&nvoxels[2], sizeof(Int));

  // 4    MODE     data type :
  //       0        image : signed 8-bit bytes range -128 to 127
  //       1        image : 16-bit halfwords
  //       2        image : 32-bit reals
  //       3        transform : complex 16-bit integers
  //       4        transform : complex 32-bit reals
  //       6        image : unsigned 16-bit range 0 to 65535
  mrc_file.write((char*)&mode, sizeof(Int));


  // 5    NXSTART number of first column in map (Default = 0)
  // 6    NYSTART number of first row in map
  // 7    NZSTART number of first section in map
  mrc_file.write((char*)&nstart[0], sizeof(Int));
  mrc_file.write((char*)&nstart[1], sizeof(Int));
  mrc_file.write((char*)&nstart[2], sizeof(Int));


  // 8    MX       number of intervals along X
  // 9    MY       number of intervals along Y
  // 10    MZ       number of intervals along Z
  mrc_file.write((char*)&mvoxels[0], sizeof(Int));
  mrc_file.write((char*)&mvoxels[1], sizeof(Int));
  mrc_file.write((char*)&mvoxels[2], sizeof(Int));


  // 11-13    CELLA    cell dimensions in angstroms
  mrc_file.write((char*)&(cellA[0]), sizeof(float));
  mrc_file.write((char*)&(cellA[1]), sizeof(float));
  mrc_file.write((char*)&(cellA[2]), sizeof(float));


  // 14-16    CELLB    cell angles in degrees
  mrc_file.write((char*)&(cellB[0]), sizeof(float));
  mrc_file.write((char*)&(cellB[1]), sizeof(float));
  mrc_file.write((char*)&(cellB[2]), sizeof(float));


  // 17    MAPC     axis corresp to cols (1,2,3 for X,Y,Z)
  // 18    MAPR     axis corresp to rows (1,2,3 for X,Y,Z)
  // 19    MAPS     axis corresp to sections (1,2,3 for X,Y,Z)
  mrc_file.write((char*)&(mapCRS[0]), sizeof(Int));
  mrc_file.write((char*)&(mapCRS[1]), sizeof(Int));
  mrc_file.write((char*)&(mapCRS[2]), sizeof(Int));

  // 20    DMIN     minimum density value
  // 21    DMAX     maximum density value
  // 22    DMEAN    mean density value
  mrc_file.write((char*)&(dmin), sizeof(float));
  mrc_file.write((char*)&(dmax), sizeof(float));
  mrc_file.write((char*)&(dmean), sizeof(float));


  // 23    ISPG     space group number 0 or 1 (default=0)
  mrc_file.write((char*)&(ispg), sizeof(Int));

  // 24    NSYMBT   number of bytes used for symmetry data (0 or 80)
  mrc_file.write((char*)&(nsymbt), sizeof(Int));

  // 25-49    EXTRA    extra space used for anything   - 0 by default
  mrc_file.write(extra_raw_data, (1+(49-25))*sizeof(Int));

  // 50-52    ORIGIN   origin in X,Y,Z used for transforms
  mrc_file.write((char*)&(origin[0]), sizeof(float));
  mrc_file.write((char*)&(origin[1]), sizeof(float));
  mrc_file.write((char*)&(origin[2]), sizeof(float));

  // I don't care about the remaining junk in the header
  mrc_file.write(remaining_raw_data,
                 MrcHeader::SIZE_REMAINING_HEADER);

} //MrcHeader::Write(ostream& mrc_file) const





void MrcHeader::PrintStats(ostream& out) {
  out << "tomogram number of voxels ("<< nvoxels[0] << ", " 
                                      << nvoxels[1] << ", "
                                      << nvoxels[2] << ")" 
                                      << endl;
  //out << "tomogram physical dimensions ("<< cellA[0] << "," << cellA[1] << "," << cellA[2] << ")" << endl;
  out << "tomogram voxel size ("
      << cellA[0]/nvoxels[0] << ", " 
      << cellA[1]/nvoxels[1] << ", " 
      << cellA[2]/nvoxels[2] << ")" << endl;
  out << "tomogram table axis order (" << mapCRS[0] << ", " 
                                       << mapCRS[1] << ", "
                                       << mapCRS[2] << ")"
                                       << endl;
  out << "tomogram mode " << mode << endl;
  out << "tomogram minimum, maximum density: " << dmin << " "
                                               << dmax<<endl;
  out << "tomogram origin ("<< origin[0] << ", " 
                            << origin[1] << ", " 
                            << origin[2] << ")" << endl;
}
