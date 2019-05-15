#ifndef _FILE_IO_HPP
#define _FILE_IO_HPP

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;
#include <err_visfd.hpp>
#include <mrc_simple.hpp>
#include "err.hpp"
#include "settings.hpp"


static vector<string> &split(const string &s, char delim, vector<string> &elems){
  stringstream ss(s);
  string item;
  while (getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}


static vector<string> split(const string &s, char delim) {
  vector<string> elems;
  split(s, delim, elems);
  return elems;
}


/// @brief   Read a multi-column text file and store the words on each line
///          in a vector of vectors.  This function can also be used to
///          read a file containing multiple columns of numbers, for example
///          a list of 3D coordinates.  ("Entry" can be numeric.)
template<typename Entry>

static void
ReadMulticolumnFile(istream &f,  //!< the file to be read
                    vector<vector<Entry> > &vvDest, //!< store the results here
                    char comment_char='#' //!<ignore text following this character ('\0' disables)
                    )
{
  vvDest.clear();
  string strLine;
  size_t i_line = 1;
  while (getline(f, strLine))
  {
    // ignore text after comments
    size_t ic = strLine.find(comment_char);
    if (ic != string::npos)
      strLine = strLine.substr(0, ic);
    stringstream ssLine(strLine);
    vector<Entry> row;
    Entry x;
    try {
      while (ssLine >> x) {
        row.push_back(x);
      }
    } // try {
    catch ( ... ) {
      stringstream err_msg;
      err_msg << "Error: File read error (invalid entry?) on line: "
              << i_line << "\n";
      throw VisfdErr(err_msg.str().c_str());
    }
    vvDest.push_back(row);
    i_line++;
  }
} // ReadMulticolumnFile()



/// @brief   Read a multi-column text file and store the words on each line
///          in a vector of vectors.  This function can also be used to
///          read a file containing multiple columns of numbers, for example
///          a list of 3D coordinates.  ("Entry" can be numeric.)
template<typename Entry>

static void
ReadMulticolumnFile(string file_name,  //!< the name of the file to be read
                    vector<vector<Entry> > &vvDest, //!< store the results here
                    char comment_char='#' //!<ignore text following this character ('\0' disables)
                    )
{
  fstream f;
  f.open(file_name.c_str(), ios::in);
  if (! f)
    throw VisfdErr("Error: unable to open \""+
                    file_name +"\" for reading.\n");
  ReadMulticolumnFile(f, vvDest, comment_char);
  f.close();
}



/// @brief  Convert a vector-of-vectors of strings
///         (which typically represent the words on each line of a file)
///         with a vector-of-vectors of numbers.
///         The text on each line might be expressed using using IMOD notation:
///         "Pixel [91.3, 118.7, 231] = 134"
///         or using ordinary whitespace-delimited notation:
///         90.3 117.7 230.
///         If it is expressed in IMOD notation, subtract 1
///         and return true (informing the caller we don't need to divide
///         the coordinates later by the voxel_width).  Otherwise return false.

template<typename Coordinate>

static bool
ConvertStringsToCoordinates(const vector<vector<string> > &vvWords_orig, //!< words on each line of a file
                            vector<vector<Coordinate> > &vvCoords, //!< convert these words to a matrix of numbers
                            int num_columns = 3 //!< number of columns in the matrix (same for each row)
                            )
{
  vvCoords.clear();

  bool is_output_from_imod = false;
  vector<vector<string> > vvWords = vvWords_orig;
  vector<Coordinate> linked_set;

  for (size_t i = 0; i < vvWords.size(); i++)
  {
    if ((vvWords[i].size() > 0) && (vvWords[i][0] == "Pixel")) {
      vvWords[i].erase(vvWords[i].begin(), vvWords[i].begin() + 1);
      is_output_from_imod = true;
    }
    vector<Coordinate> xyz;

    if (vvWords[i].size() == 0) {
      vvCoords.push_back(vector<Coordinate>(0));
      continue;
    }

    for (int d=0; d < num_columns; d++) {
      // Get rid if weird characters enclosing or separating the words
      // such as '(', ')', or ','
      if ((vvWords[i][d].size() >= 0) && (vvWords[i][d][0] == '(')) {
        // Delete the '(' character from the beginning of the string.
        // (IMOD surrounds the coordinates with '(' and ')' characters.)
        vvWords[i][d] = vvWords[i][d].substr(1, string::npos);
        is_output_from_imod = true;
      }
      int n_chars;
      n_chars = vvWords[i][d].size();
      if ((vvWords[i][d].size()>=0) && (vvWords[i][d][n_chars-1]==')')) {
        // Delete the ')' character from the end of the string.
        // (IMOD surrounds the coordinates with '(' and ')' characters.)
        vvWords[i][d] = vvWords[i][d].substr(0, n_chars-1);
        is_output_from_imod = true;
      }
      n_chars = vvWords[i][d].size();
      if ((vvWords[i][d].size()>=0) && (vvWords[i][d][n_chars-1]==','))
        vvWords[i][d] = vvWords[i][d].substr(0, n_chars-1);

      // Now convert the remaining string to a number
      try {
        double x = stof(vvWords[i][d]); //<-- coordinate (x, y, or z)
      }
      catch ( ... ) {
        stringstream err_msg;
        err_msg << "Error: File read error (invalid entry?) on line: "
                << i+1 << "\n";
        throw VisfdErr(err_msg.str().c_str());
      }

      // The way we interpret the number depends on whether or not
      // we suspect the number was printed by IMOD, which uses
      // a somewhat unconventional coordinate system
      Coordinate x = stod(vvWords[i][d]);
      if (is_output_from_imod)
        x = floor(x) - 1.0;

      // Finally, store the coordinate in "xyz"
      xyz.push_back(x);

    } //for (int d=0; d < vvWords[i].size(); d++)

    vvCoords.push_back(xyz);

  } //for (size_t i = 0; i < vvWords.size(); i++)

  assert(vvCoords.size() == vvWords.size());


  return is_output_from_imod;

} // ConvertStringsToCoordinates()



/// @brief  Read a file containing a list of coordinates.
///         The text on each line might be expressed using using IMOD notation:
///         "Pixel [91.3, 118.7, 231] = 134"
///         or using ordinary whitespace-delimited notation:
///         90.3 117.7 230.
///         If it is expressed in IMOD notation, subtract 1, 
///         but do not divide by the voxel_width.  (IMOD coordinates are in
///         units of voxels, not physical distance, however indices start at 1).
///         Otherwise, simply divide the numbers by the voxel width.

template<typename Coordinate>
bool
ReadCoordinates(string filename,  //!< the name of the file to be read
                vector<array<Coordinate, 3>> &vaCoords, //!< store the coordinates here
                char comment_char='#' //!<ignore text following this character ('\0' disables)
                )
{
  vector<vector<string> > vvWords;  // words on each line
  vector<vector<Coordinate> > vvCoords; // replace each word with a number

  ReadMulticolumnFile(filename, vvWords, comment_char);

  for (size_t i = 0; i < vvWords.size(); i++) {
    if (vvWords.size() == 0)
      vvWords.erase(vvWords.begin()+i, vvWords.begin()+i+1);
    else if (vvWords.size() < 3) {
      stringstream err_msg;
      err_msg << "Format error near line "<<i+1<<" of file \""+
                     filename+"\"\n";
      throw InputErr(err_msg.str());
    }
  }

  bool is_in_voxels =
    ConvertStringsToCoordinates(vvWords,   //convert words to numbers
                                vvCoords); // store coordinates here

  vaCoords.clear();
  for (size_t i = 0; i < vvCoords.size(); i++) {
    if (vvCoords[i].size() >= 3) {
      array<float, 3> xyz = {vvCoords[i][0], vvCoords[i][1], vvCoords[i][2]};
      vaCoords.push_back(xyz);
    }
  }
  return is_in_voxels;
} //ReadCoordinates()





/// @brief   Read a list of blob coordinates, diameters, and scores
///          from a text file.

template<typename Scalar, typename Coordinate>

static void
ReadBlobCoordsFile(string in_coords_file_name, //!< name of file we will read
                   vector<array<Coordinate, 3> > *pCrds=nullptr, //!< store the blob coordinates here (if !=nullptr)
                   vector<Scalar> *pDiameters=nullptr, //!< store the blob diameters here (if !=nullptr)
                   vector<Scalar> *pScores=nullptr, //!< store blob scores here (if !=nullptr)
                   Scalar distance_scale=1.0, //!< divide all distances and coordinates by this value
                   Scalar diameter_override=-1.0, //!< use this diameter (useful if no 4th column is present)
                   Scalar score_default=0.0, //!< default "score" (if no 5th column is present)
                   Scalar diameter_factor=1.0 //!< multiply all diameters in the file by this number
                   )
{
  fstream coords_file;
  coords_file.open(in_coords_file_name.c_str(), ios::in);
  if (! coords_file)
    throw VisfdErr("Error: unable to open \""+
                   in_coords_file_name +"\" for reading.\n");

  bool custom_diameters = false;
  while (coords_file) {
    string strLine;
    getline(coords_file, strLine);
    if (strLine.size() == 0)
      continue;
    stringstream ssLine(strLine);
    double x, y, z;
    ssLine >> x;
    ssLine >> y;
    ssLine >> z;
    double ix, iy, iz;
    ix = floor((x / distance_scale) + 0.5);
    iy = floor((y / distance_scale) + 0.5);
    iz = floor((z / distance_scale) + 0.5);
    array<Coordinate, 3> ixiyiz;
    ixiyiz[0] = ix;
    ixiyiz[1] = iy;
    ixiyiz[2] = iz;

    Scalar diameter = -1.0;
    Scalar score = score_default;
    if (ssLine) { // Does the file contain a 4th column? (the diameter)
      Scalar _diameter;
      ssLine >> _diameter;
      // convert from physical distance to # of voxels:
      _diameter /= distance_scale;
      if (ssLine)
        diameter = _diameter;
      custom_diameters = true;
    }

    if (diameter < 0) { //If file does not contain a 4th column, or if < 0
      diameter = diameter_override;
      if (diameter_override < 0)
        diameter = 1.0;   // (sphere will be 1 voxel wide by default)
    }

    if (diameter_override >= 0) //override the diameter ?
      diameter = diameter_override;
    else
      diameter *= diameter_factor;

    if (ssLine) {
      Scalar _score;
      ssLine >> _score;
      if (ssLine)
        score = _score;
    }

    if (pCrds)
      (*pCrds).push_back(ixiyiz);
    if (pDiameters)
      (*pDiameters).push_back(diameter);
    if (pScores)
      (*pScores).push_back(score);

  } //while (coords_file) {...
} //ReadBlobCoordsFile()




template<typename Scalar>

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




template<typename Scalar>

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



template<typename Scalar>

static void
WriteOrientedPointCloudBNPTS(string filename,
                             vector<array<Scalar,3> > coords,
                             vector<array<Scalar,3> > norms)
{
  assert(coords.size() == norms.size());
  throw VisfdErr("The WriteOrientedPointCloudBNPTS() function is not working correctly.\n"
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




template<typename Scalar, typename VectorContainer>

static void
WriteOrientedPointCloud(string pointcloud_file_name,
                        const int image_size[3],
                        VectorContainer const *const *const *aaaafVector,
                        Scalar const *const *const *aaafMask = nullptr,
                        const Scalar *voxel_width=nullptr)
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





/// @brief   Parse the list of voxel coordinates provided 
///          by the user to use in "link" constraints.
/// @return  whether or not the coordinates are in units of voxels.
///          (If not, they will need to be converted to these units later
///           by dividing by the physical width represented by each voxel.)

template<typename Coordinate>
static bool
ProcessLinkConstraints(string must_link_filename,
                       vector<vector<array<Coordinate, 3> > > &must_link_constraints)
{
  vector<vector<string> > vvWords;  // split file into lines and words
  vector<vector<Coordinate> > vvCoords; // replace each word with a number
  vector<array<Coordinate, 3> > linked_group;

  ReadMulticolumnFile(must_link_filename, vvWords); //read the words

  bool is_in_voxels = 
    ConvertStringsToCoordinates(vvWords, //convert words to numbers
                                vvCoords);

  // Now generate a vector of linked groups
  for (size_t i = 0; i < vvCoords.size(); i++) {
    if (vvCoords[i].size() == 0) {
      if (linked_group.size() > 0)
        must_link_constraints.push_back(linked_group);
      linked_group.clear();
    }
    else if (vvCoords[i].size() == 3) {
      array<Coordinate, 3> voxel_location;
      for (int d = 0; d < 3; d++)
        voxel_location[d] = vvCoords[i][d];
      linked_group.push_back(voxel_location);
    }
    else {
      stringstream err_msg;
      err_msg << "Error: Each line of file \""
              << must_link_filename << "\"\n"
              <<"       should contain either 3 numbers or 0 numbers.\n";
      throw VisfdErr(err_msg.str());
    }
  } //for (size_t i = 0; i < vvCoords.size(); i++)

  // Deal with any left-over linked_groups that we haven't added to the list yet
  if (linked_group.size() > 0)
    must_link_constraints.push_back(linked_group);

  if (must_link_constraints.size() == 0) {
    stringstream err_msg;
    err_msg << "Error: Format error in file \""<<must_link_filename<<"\".\n"
            << "       File contains no voxel coordinates.\n";
    throw VisfdErr(err_msg.str());
  }

  for (size_t i = 0; i < must_link_constraints.size(); i++) {
    linked_group = must_link_constraints[i];
    if ((linked_group.size() < 2) || (linked_group[0] == linked_group[1])) {
      stringstream err_msg;
      err_msg
        << "Error: Format error in file \""<<must_link_filename<<"\".\n"
        << "       Each broup must contain at least 2 voxels.  (Voxels appear on different\n"
        << "       lines, so blank-line delimters must not separate SINGLE non-blank lines)\n"
        << "       Furthermore, the voxels in each set must be unique.\n";

      throw VisfdErr(err_msg.str());
    }
  }

  return is_in_voxels;

} //ProcessLinkConstraints()



#endif //#ifndef _FILE_IO_HPP
