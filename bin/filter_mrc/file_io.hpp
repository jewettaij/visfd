#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

#include <err_report.hpp>


template<class Scalar>
static void
WriteOrientedPointCloudPLY(string filename,
                           vector<array<Scalar,3> > coords,
                           vector<array<Scalar,3> > norms)
{
  assert(coords.size() == norms.size());
  
  size_t n = coords.size();
  fstream ply_file;
  ply_file.open(filename.c_str(), ios::out);
  ply_file <<
    "ply\n"
    "format ascii 1.0\n"
    "comment  created by visfd\n"
    "element vertex " << n << "\n"
    "property float x\n"
    "property float y\n"
    "property float z\n"
    "property float nx\n"
    "property float ny\n"
    "property float nz\n"
    //"property list uchar int vertex_index\n" <- needed for "faces" only
    "end_header\n";
           
  for (size_t i=0; i < n; i++)
    ply_file << coords[i][0] << " "<< coords[i][1] << " " << coords[i][2] <<" "
             << norms[i][0] << " " << norms[i][1] << " " << norms[i][2] <<"\n";
  ply_file.close();
}




template<class Scalar>
static void
WriteOrientedPointCloudOBJ(string filename,
                           vector<array<Scalar,3> > coords,
                           vector<array<Scalar,3> > norms)
{
  assert(coords.size() == norms.size());
  
  size_t n = coords.size();
  fstream obj_file;
  obj_file.open(filename.c_str(), ios::out);
  obj_file << "# WaveFront *.obj file created by visfd\n";
  obj_file << "\n"
           << "g obj1_\n"
           << "\n";
           
  for (size_t i=0; i < n; i++)
    obj_file << "v "
             << coords[i][0] <<" "<< coords[i][1] <<" "<< coords[i][2] <<"\n";

  obj_file << "\n";

  for (size_t i=0; i < n; i++)
    obj_file << "vn "
             << norms[i][0] <<" "<< norms[i][1] <<" "<< norms[i][2] <<"\n";

  obj_file.close();
}



template<class Scalar>
static void
WriteOrientedPointCloudBNPTS(string filename,
                             vector<array<Scalar,3> > coords,
                             vector<array<Scalar,3> > norms)
{
  assert(coords.size() == norms.size());
  throw InputErr("The WriteOrientedPointCloudBNPTS() function is not working correctly.\n"
                 "No code should be invoking this function until the bug(s) is fixed.\n"
                 "If you are seeing this message, the error is the fault of the programmer.\n"
                 "Please contact the developer.");
  size_t n = coords.size();
  fstream bnpts_file;
  bnpts_file.open(filename.c_str(), ios::out | ios::binary);
  for (size_t i=0; i < n; i++) {
    float xyz[3];           //(convert from "Scalar" to float)
    xyz[0] = coords[i][0];  //(convert from "Scalar" to float)
    xyz[1] = coords[i][1];
    xyz[2] = coords[i][2];
    float norm[3];          //(convert from "Scalar" to float)
    norm[0] = norms[i][0];  //(convert from "Scalar" to float)
    norm[1] = norms[i][1];
    norm[2] = norms[i][2];
    bnpts_file.write((char*)&(xyz[0]), sizeof(float));
    bnpts_file.write((char*)&(xyz[1]), sizeof(float));
    bnpts_file.write((char*)&(xyz[2]), sizeof(float));
    bnpts_file.write((char*)&(norm[0]), sizeof(float));
    bnpts_file.write((char*)&(norm[1]), sizeof(float));
    bnpts_file.write((char*)&(norm[2]), sizeof(float));
  }
  bnpts_file.close();
}




template<class Scalar, class VectorContainer>
static void
WriteOrientedPointCloud(string pointcloud_file_name,
                        const int image_size[3],
                        VectorContainer const *const *const *aaaafVector,
                        Scalar const *const *const *aaafMask = NULL,
                        const Scalar *voxel_width=NULL)
{
  assert(aaaafVector);
  vector<array<Scalar,3> > coords;
  vector<array<Scalar,3> > norms;

  // Did the caller specify the physical width of each voxel?
  Scalar _voxel_width[3] = {1.0, 1.0, 1.0};
  if (voxel_width) {
    _voxel_width[0] = voxel_width[0];
    _voxel_width[1] = voxel_width[1];
    _voxel_width[2] = voxel_width[2];
  }

  for (int iz=0; iz < image_size[2]; iz++) {
    for (int iy=0; iy < image_size[1]; iy++) {
      for (int ix=0; ix < image_size[0]; ix++) {
        if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
          continue;

        array<Scalar,3> xyz;
        xyz[0] = ix * _voxel_width[0];
        xyz[1] = iy * _voxel_width[1];
        xyz[2] = iz * _voxel_width[2];
        coords.push_back(xyz);

        array<Scalar,3> norm;
        norm[0] = aaaafVector[iz][iy][ix][0];
        norm[1] = aaaafVector[iz][iy][ix][1];
        norm[2] = aaaafVector[iz][iy][ix][2];
        norms.push_back(norm);
      }
    }
  }

  //WriteOrientedPointCloudBNPTS(pointcloud_file_name + ".bnpts",
  //                             coords,
  //                             norms);

  //WriteOrientedPointCloudOBJ(pointcloud_file_name,
  //                           coords,
  //                           norms);

  WriteOrientedPointCloudPLY(pointcloud_file_name,
                             coords,
                             norms);

} //WriteOrientedPointCloud()

