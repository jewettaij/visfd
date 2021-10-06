#include <string>
#include <iostream>
using namespace std;
#include "err.hpp"
#include "mrc_simple.hpp"



string g_program_name = "crop_mrc";

int main(int argc, char **argv) {
  if (argc <= 8) {
    cerr << "This program expects 8 arguments.  Syntax:\n"
	 << "  " << g_program_name << " input_file output_file xmin xmax ymin ymax zmin zmax\n"
         << "Exiting.\n" << endl;
    exit(0);
  }

  try {
    string in_file_name(argv[1]);
    string out_file_name(argv[2]);
    int xmin = atoi(argv[3]);
    int xmax = atoi(argv[4]);
    if (xmax < xmin) throw InputErr("Error xmin > xmax\n");
    int ymin = atoi(argv[5]);
    int ymax = atoi(argv[6]);
    if (ymax < ymin) throw InputErr("Error ymin > ymax\n");
    int zmin = atoi(argv[7]);
    int zmax = atoi(argv[8]);
    if (zmax < zmin) throw InputErr("Error zmin > zmax\n");

    // Read the file
    cerr << "Reading tomogram \""<<in_file_name<<"\"" << endl;
    MrcSimple tomo;
    tomo.Read(in_file_name, false);//You can also use "tomo.Read(cin);" or "cin>>tomo;"
    tomo.PrintStats(cerr);  //Optional (display the tomogram size & format)

    if (xmin < 0) {
      xmin = 0;
      cerr << "####################################################\n"
           << "######### WARNING: xmin changed from " << xmin << " to 0\n"
           << "####################################################" << endl;
    }
    if (xmax >= tomo.header.nvoxels[0]) {
      cerr << "####################################################\n"
           << "######### WARNING: xmax changed from " << xmax << " to "
           << tomo.header.nvoxels[0]-1 << "\n"
           << "####################################################" << endl;
      xmax = tomo.header.nvoxels[0] - 1;
    }
    if (ymin < 0) {
      ymin = 0;
      cerr << "####################################################\n"
           << "######### WARNING: ymin changed from " << ymin << " to 0\n"
           << "#####################################################" << endl;
    }
    if (ymax >= tomo.header.nvoxels[1]) {
      cerr << "#####################################################\n"
           << "######### WARNING: xmax changed from " << ymax << " to "
           << tomo.header.nvoxels[1]-1 << "\n"
           << "#####################################################" << endl;
      ymax = tomo.header.nvoxels[1] - 1;
    }
    if (zmin < 0) {
      zmin = 0;
      cerr << "#####################################################\n"
           << "######### WARNING: zmin changed from " << zmin << " to 0\n"
           << "#####################################################" << endl;
    }
    if (zmax >= tomo.header.nvoxels[2]) {
      cerr << "##################################################\n"
           << "######### WARNING: xmax changed from " << zmax << " to "
           << tomo.header.nvoxels[2]-1 << "\n"
           << "##################################################" << endl;
      zmax = tomo.header.nvoxels[2] - 1;
    }

    MrcSimple cropped_tomo;
    cropped_tomo.header = tomo.header;

    float voxel_size[3];
    voxel_size[0] = cropped_tomo.header.cellA[0] / cropped_tomo.header.nvoxels[0];
    voxel_size[1] = cropped_tomo.header.cellA[1] / cropped_tomo.header.nvoxels[1];
    voxel_size[2] = cropped_tomo.header.cellA[2] / cropped_tomo.header.nvoxels[2];

    // Optional: reduce the cellA[0], cellA[1], cellA[2] size due to cropping
    int nvoxels[3];
    cropped_tomo.header.cellA[0] *= (1.0 + xmax-xmin) / cropped_tomo.header.nvoxels[0];
    cropped_tomo.header.cellA[1] *= (1.0 + ymax-ymin) / cropped_tomo.header.nvoxels[1];
    cropped_tomo.header.cellA[2] *= (1.0 + zmax-zmin) / cropped_tomo.header.nvoxels[2];
    nvoxels[0] = 1 + xmax-xmin;
    nvoxels[1] = 1 + ymax-ymin;
    nvoxels[2] = 1 + zmax-zmin;
    cropped_tomo.Resize(nvoxels);
    cropped_tomo.header.nvoxels[0] = nvoxels[0];
    cropped_tomo.header.nvoxels[1] = nvoxels[1];
    cropped_tomo.header.nvoxels[2] = nvoxels[2];
    //           and shift origin[0],origin[1],origin[2] accordingly
    cropped_tomo.header.origin[0]= tomo.header.origin[0]+(xmin-1)*voxel_size[0];
    cropped_tomo.header.origin[1]= tomo.header.origin[1]+(ymin-1)*voxel_size[1];
    cropped_tomo.header.origin[2]= tomo.header.origin[2]+(zmin-1)*voxel_size[2];

    for (int iz = 0; iz <= zmax-zmin; iz++)
      for (int iy = 0; iy <= ymax-ymin; iy++)
	for (int ix = 0; ix <= xmax-xmin; ix++)
	  cropped_tomo.aaafI[iz][iy][ix] = 
	    tomo.aaafI[iz+zmin][iy+ymin][ix+xmin];

    cropped_tomo.FindMinMaxMean(); //update the min, max, mean header params

    // Write the file
    cerr << "writing cropped tomogram (in float mode)" << endl;
    cerr << "  new tomogram after cropping:" << endl;
    cropped_tomo.PrintStats(cerr);
    cropped_tomo.Write(out_file_name); //(You can also use cout<<cropped_tomo;)

  } //try {
  catch (const std::exception& e) {
    cerr << "\n" << e.what() << endl;
    exit(1);
  }
} // main()

