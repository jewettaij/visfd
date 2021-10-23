#ifndef _FILE_IO_HPP
#define _FILE_IO_HPP

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;
#include <connect.hpp>
using namespace visfd;
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



/// @brief split the text in "line" into white-space or comma delimited words
static void
Str2Words(string source,          //!< line of text containing words
          vector<string> &vDest, //!< store the vector of strings here
          char comment_char='#'   //!<ignore text following this character ('\0' disables)
          )
{
  vDest.clear();
  stringstream ssSource(source);
  char c;
  string word;
  bool preceeded_by_comma = false;
  while (ssSource.get(c)) {
    if (isspace(c) || (c == ',') || (c == comment_char)) {
      if ((word.size() > 0) ||
          ((c == ',') && preceeded_by_comma)) {
        vDest.push_back(word);
        word = string(""); //clear word
        preceeded_by_comma = false;
        if (c == ',')
          preceeded_by_comma = true;
      }
      if (c == comment_char)
        break;
    }
    else
      word.push_back(c); // add c to the end of this word
  }
  // now deal with text in the last word in the source
  if ((word.size() > 0) || preceeded_by_comma)
    vDest.push_back(word);
} //Str2Words()



/// @brief  Convert a vector of strings
///         (which typically represent the words on one line in a file)
///         with a vector of numbers.
///         The text on each line might be expressed using using IMOD notation:
///         "Pixel (91.3, 118.7, 231) = 134"   or
///         "(91.3, 118.7, 231)"
///         ...or using ordinary whitespace-delimited notation:
///         "90.3 117.7 230"
///         If any of the rows of strings appear to use IMOD notation
///         (eg. "Pixel (241, 315, 82) = 12.3"))
///         then extract the 3 numbers in parenthesis (and subtract 1 from each)
///         and discard the rest.
/// @return If IMOD notation was used (eg, if the first word is "Pixel")
///         then eturn true.  Otherwise return false.

template<typename Coordinate>
static bool
IMODWords2Crds(const vector<string> &vWords, //!< words containing x,y,z coordinates
               vector<Coordinate> &xyz,  //!< store the numbers here
               char comment_char='#' //!<ignore text following this character ('\0' disables)
               )
{
  // (I really hate this redundant ugly code.
  //  I would reorganize this if I had time, but it's working at the moment. -A)

  xyz.clear();

  bool is_output_from_imod = false;
  bool contains_parens = false;
  vector<string> vWords_cpy = vWords;

  // First, if the line of text contains a comment character,
  // skip over the remaining text in that line.
  for (size_t i = 0; i < vWords_cpy.size(); i++) {
    string word = vWords_cpy[i];
    // ignore text after comments
    size_t ic = word.find(comment_char);
    if (ic != string::npos) {
      word = word.substr(0, ic);  //discard remaining characters in this word
      vWords_cpy[i] = word;
      if (i < vWords_cpy.size())  //discard the remaining words on this line
        vWords_cpy.erase(vWords_cpy.begin()+i+1, vWords_cpy.end());
    }
  }

  if ((vWords_cpy.size() > 0) && (vWords_cpy[0] == "Pixel")) {
      vWords_cpy.erase(vWords_cpy.begin(), vWords_cpy.begin() + 1);
      is_output_from_imod = true;
      contains_parens = true;
  }

  if (vWords_cpy.size() == 0) {
    //xyz.push_back(vector<Coordinate>(0));
    xyz.clear();
    return is_output_from_imod;
  }

  int d = 0;
  while (d < vWords_cpy.size()) {
    // Get rid if weird characters enclosing or separating the words
    // such as '(', ')', or ','
    if ((vWords_cpy[d].size() >= 0) && (vWords_cpy[d][0] == '(')) {
      // Delete the '(' character from the beginning of the string.
      // (IMOD surrounds the coordinates with '(' and ')' characters.)
      contains_parens = true;
      vWords_cpy[d] = vWords_cpy[d].substr(1, string::npos);
      if (vWords_cpy[d].size() == 0) {
        vWords_cpy.erase(vWords_cpy.begin()+d, vWords_cpy.begin()+d+1);
        continue;
      }
    }
    int n_chars;
    n_chars = vWords_cpy[d].size();
    if ((vWords_cpy[d].size()>=0) && (vWords_cpy[d][n_chars-1]==')')) {
      // Delete the ')' character from the end of the string.
      // (IMOD surrounds the coordinates with '(' and ')' characters.)
      contains_parens = true;
      vWords_cpy[d] = vWords_cpy[d].substr(0, n_chars-1);
      if (vWords_cpy[d].size() == 0) {
        vWords_cpy.erase(vWords_cpy.begin()+d, vWords_cpy.begin()+d+1);
        continue;
      }
    }

    if ((n_chars>0) && (vWords_cpy[d][n_chars-1]==',')) {
      vWords_cpy[d] = vWords_cpy[d].substr(0, n_chars-1);
      if (vWords_cpy[d].size() == 0) {
        vWords_cpy.erase(vWords_cpy.begin()+d, vWords_cpy.begin()+d+1);
        continue;
      }
    }

    // Is this "word" actually a comma-separated list of numbers?
    // If so, split this string into multiple words.
    {
      stringstream word_ss(vWords_cpy[d]);
      string token;
      vector<string> tokens;
      while (getline(word_ss, token, ',')) {
        //if (token.size() > 0)
        tokens.push_back(token);
      }
      if (tokens.size() > 1) {
        vWords_cpy[d] = tokens[0];
        vWords_cpy.insert(vWords_cpy.begin()+d+1,
                          tokens.begin()+1, tokens.end());
      }
    }

    // Now convert the remaining string to a number
    Coordinate x; //<-- shorthand for all 3 coordinates (x, y, z)
    // The way we interpret the number depends on whether or not
    // we suspect the number was printed by IMOD, which uses
    // a coordinate system with integer coordinates beginning at 1, not 0.
    if ((d < 3) || (! is_output_from_imod))
    {
      stringstream word_ss(vWords_cpy[d]);
      try {
        word_ss >> x;
        // To be as general as possible, I use >> instead of stod(vWords_cpy[d])
        // (I don't assume "Coordinate" is of type "double".)
      }
      catch ( ... ) {
        stringstream err_msg;
        err_msg << "Error: File read error (invalid entry?) on line:\n"
                << "     ";
        for (int i = 0; i < vWords.size(); i++)
          err_msg << " " << vWords[i];
        err_msg << "\n";
        throw VisfdErr(err_msg.str());
      }
      if (contains_parens || is_output_from_imod)
        // we want x,y,z indices to begin at 0 not 1,
        // but in IMOD, x,y,z indices begin at 1.  So subtract 1.
        x = floor(x) - 1.0;
      // Finally, store the coordinate in "xyz"
      xyz.push_back(x);
    }

    d++;
  } //for (int d=0; d < vWords_cpy.size(); d++)

  return contains_parens;

} // IMODWords2Crds()



/// @brief  Convert a string into a 1-D table of numbers.
///         (The string typically corresponds to a line from a text file.)
///         If the string appears to use IMOD notation
///         (eg. "Pixel (241, 315, 82) = 12.3"))
///         then extract the 3 numbers in parenthesis, and discard the rest.
/// @return Returns true if the string uses IMOD notation.

template<typename Coordinate>
static bool
IMODStr2Crds(string line,             //!< line of text containing x,y,z coordinates
             vector<Coordinate> &xyz, //!< convert these words to a vector of numbers
             char comment_char='#'    //!<ignore text following this character ('\0' disables)
             )
{
  // split the text in "line" into white-space delimited words
  vector<string> vWords;
  Str2Words(line, vWords, comment_char);
  return IMODWords2Crds(vWords, xyz, comment_char);
}



/// @brief  Convert a 2-dimensional table of strings into a table of numbers.
///         (The words on each row in the table are typically parsed from
///          different linse in a multi-column text file.)
///         If any of the rows of strings appear to use IMOD notation
///         (eg. "Pixel (241, 315, 82) = 12.3"))
///         then extract the 3 numbers in parenthesis, and discard the rest.
/// @return Returns true if any of the vectors of strings use IMOD notation.
template<typename Coordinate>
static bool
IMODWords2Crds(const vector<vector<string> > &vvWords_orig, //!< 2D table of strings
               vector<vector<Coordinate> > &vvCoords,  //!< store the resulting 2D table of numbers here
               char comment_char='#' //!<ignore text following this character ('\0' disables)
               )
{
  vvCoords.clear();

  bool is_output_from_imod = false;
  vector<vector<string> > vvWords = vvWords_orig;

  for (size_t i = 0; i < vvWords.size(); i++)
  {
    vector<Coordinate> xyz;
    if (IMODWords2Crds(vvWords_orig[i], xyz))
      is_output_from_imod = true;
    if (xyz.size() > 0)
      vvCoords.push_back(xyz);
  } //for (size_t i = 0; i < vvWords.size(); i++)

  //assert(vvCoords.size() == vvWords.size());

  return is_output_from_imod;

} // IMODWords2Crds()




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
  string line;
  size_t i_line = 1;
  while (getline(f, line))
  {
    // ignore text after comments <--COMMENTING OUT. WE'LL HANDLE COMMENTS LATER
    //size_t ic = line.find(comment_char);
    //if (ic != string::npos)
    //  line = line.substr(0, ic);

    // Split this line of text into words separated by commas and/or spaces.
    vector<string> vWords;
    Str2Words(line, vWords, comment_char);
 
    if (vWords.size() == 0)
      continue;

    vector<Entry> vDest;

    try {
      for (int i = 0; i < vWords.size(); i++) {
        stringstream word_ss(vWords[i]);
        Entry x;
        word_ss >> x;
        vDest.push_back(x);
      }
    } // try {
    catch ( ... ) {
      stringstream err_msg;
      err_msg << "Error: File read error (invalid entry?) on line: "
              << i_line << "\n";
      throw VisfdErr(err_msg.str());
    }

    vvDest.push_back(vDest);
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
  f.open(file_name, ios::in);
  if (! f)
    throw VisfdErr("Error: unable to open \""+
                    file_name +"\" for reading.\n");
  ReadMulticolumnFile(f, vvDest, comment_char);
  f.close();
}





/// @brief  Read a file containing a list of coordinates.
///         The text on each line might be expressed using using IMOD notation:
///         "Pixel (91.3, 118.7, 231) = 134"
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
    if (vvWords[i].size() == 0)
      vvWords.erase(vvWords.begin()+i, vvWords.begin()+i+1);
    else if (vvWords[i].size() < 3) {
      stringstream err_msg;
      err_msg << "Format error near line "<<i+1<<" of file \""
              << filename+"\"\n";
      throw InputErr(err_msg.str());
    }
  }

  bool is_in_voxels =
    IMODWords2Crds(vvWords,  //convert words to numbers
                   vvCoords, // store coordinates here
                   comment_char);

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
/// @returns true if any of the file coordinates were surrounded by parenthesis.
///          (For example: "(23 51 49)" or "Pixel (411, 196, 171) = 63".
///          In that case, it is assumed the numbers were printed by IMOD
///          and are in units of voxels instead of units of physical distance.)

template<typename Scalar, typename Coordinate>
static bool
ReadBlobCoordsFile(string in_coords_file_name, //!< name of file we will read
                   vector<array<Coordinate, 3> > *pCrds=nullptr, //!< store the blob coordinates here (if !=nullptr)
                   vector<Scalar> *pDiameters=nullptr, //!< store the blob diameters here (if !=nullptr)
                   vector<Scalar> *pScores=nullptr, //!< store blob scores here (if !=nullptr)
                   Scalar diameter_override=-1.0, //!< use this diameter (useful if no 4th column is present)
                   Scalar score_default=0.0, //!< default "score" (if no 5th column is present)
                   Scalar diameter_factor=1.0, //!< multiply all diameters in the file by this number
                   char comment_char='#' //!<ignore text following this character ('\0' disables)
                   )
{
  bool has_parens = false;
  fstream coords_file;
  coords_file.open(in_coords_file_name, ios::in);
  if (! coords_file)
    throw VisfdErr("Error: unable to open \""+
                   in_coords_file_name +"\" for reading.\n");

  size_t i_line = 0;

  //bool custom_diameters = false;
  while (coords_file) {

    string line;
    getline(coords_file, line);

    // Convert this line of text into an array of numbers
    vector<Coordinate> vCoords;
    if (IMODStr2Crds(line,    //convert words to numbers
                     vCoords, // store coordinates here
                     comment_char))
      has_parens = true;
    //(is_in_voxels warns us that IMOD notation was used so units are in voxels)

    if (vCoords.size() == 0)
      continue;

    if (! ((vCoords.size() == 3) ||
           (vCoords.size() == 4) ||
           (vCoords.size() == 5)))
    {
      stringstream err_msg;
      err_msg << "Error: Error on line " << i_line+1 << " of file \""
              << in_coords_file_name << "\"\n"
              << "       Each line should contain either 3-5 numbers, or 0 numbers (blank).\n";
      throw VisfdErr(err_msg.str());
    }

    array<Coordinate, 3> xyz;
    xyz[0] = vCoords[0];
    xyz[1] = vCoords[1];
    xyz[2] = vCoords[2];

    Scalar diameter = -1.0;
    if (vCoords.size() > 3) {
      diameter = vCoords[3];
      //custom_diameters = true;
    }
    if (diameter < 0) { //If file does not contain a 4th column, or if < 0
      diameter = diameter_override;
    }

    if (diameter_override >= 0) //override the diameter ?
      diameter = diameter_override;
    else
      diameter *= diameter_factor;

    Scalar score = score_default;
    if (vCoords.size() > 4)
      score = vCoords[4];

    if (pCrds)
      (*pCrds).push_back(xyz);
    if (pDiameters)
      (*pDiameters).push_back(diameter);
    if (pScores)
      (*pScores).push_back(score);

    i_line++;
  } //while (coords_file) {...
  return has_parens;
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
  ply_file.open(filename, ios::out);
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
  obj_file.open(filename, ios::out);
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
  bnpts_file.open(filename, ios::out | ios::binary);
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
                       vector<vector<array<Coordinate, 3> > > &must_link_constraints,
                       vector<vector<DirectionPairType> > &must_link_directions)
{
  vector<vector<string> > vvWords;  // split file into lines and words
  vector<vector<Coordinate> > vvCoords; // replace each word with a number
  vector<array<Coordinate, 3> > linked_group;
  vector<DirectionPairType> linked_group_directions;

  ReadMulticolumnFile(must_link_filename, vvWords); //read the words

  bool is_in_voxels = 
    IMODWords2Crds(vvWords, //convert words to numbers
                   vvCoords);

  // Now generate a vector of linked groups
  for (size_t i = 0; i < vvCoords.size(); i++) {
    if (vvCoords[i].size() == 0) {
      if (linked_group.size() > 0) {
        assert(linked_group.size() == linked_group_directions.size());
        must_link_constraints.push_back(linked_group);
        must_link_directions.push_back(linked_group_directions);
      }
      linked_group.clear();
      linked_group_directions.clear();
    }
    else if ((vvCoords[i].size() == 3) ||
             (vvCoords[i].size() == 4)) {
      array<Coordinate, 3> voxel_location;
      for (int d = 0; d < 3; d++)
        voxel_location[d] = vvCoords[i][d];
      linked_group.push_back(voxel_location);
      DirectionPairType flip_direction = AUTO;
      if (vvCoords[i].size() == 4) {
        if (vvCoords[i][3] > 0)
          flip_direction = SAME_DIRECTION;
        else if (vvCoords[i][3] < 0)
          flip_direction = OPPOSITE_DIRECTION;
      }
      linked_group_directions.push_back(flip_direction);
    }
    else {
      stringstream err_msg;
      err_msg << "Error: Each line of file \""
              << must_link_filename << "\"\n"
              <<"       should contain either 3 numbers, 4 numbers, or 0 numbers.\n";
      throw VisfdErr(err_msg.str());
    }
  } //for (size_t i = 0; i < vvCoords.size(); i++)

  // Deal with any left-over linked_groups that we haven't added to the list yet
  if (linked_group.size() > 0) {
    assert(linked_group.size() == linked_group_directions.size());
    must_link_constraints.push_back(linked_group);
    must_link_directions.push_back(linked_group_directions);
  }

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
        << "       Each group must contain at least 2 voxels.  (Voxels appear on different\n"
        << "       lines, so blank-line delimters must not separate SINGLE non-blank lines)\n"
        << "       Furthermore, the voxels in each set must be unique.\n";

      throw VisfdErr(err_msg.str());
    }
  }

  return is_in_voxels;

} //ProcessLinkConstraints()



#endif //#ifndef _FILE_IO_HPP
