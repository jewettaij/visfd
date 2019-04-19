///   @file feature.hpp
///   @brief  THIS CODE IS NOT INTENDED FOR PUBLIC USE
///           (The functions defined here are invoked in "feature.hpp")
///   @author Andrew Jewett
///   @date 2019-4-17

#ifndef _DO_NOT_USE_FEATURE_HPP
#define _DO_NOT_USE_FEATURE_HPP

#include <cstring>
#include <cassert>
#include <cmath>
#include <limits>
#include <ostream>
#include <vector>
#include <tuple>
#include <set>
#include <queue>
using namespace std;
#include <err_visfd.hpp> // defines the "VisfdErr" exception type
#include <eigen3_simple.hpp>  // defines matrix diagonalizer (DiagonalizeSym3())
#include <visfd_utils.hpp>    // defines invert_permutation(), AveArray(), ...
#include <alloc2d.hpp>    // defines Alloc2D() and Dealloc2D()
#include <alloc3d.hpp>    // defines Alloc3D() and Dealloc3D()
#include <filter1d.hpp>   // defines "Filter1D" (used in ApplySeparable())
#include <filter3d.hpp>   // defines common 3D image filters



namespace visfd {




/// @brief  The following function finds EITHER local minima OR local maxima
///         in an image.
///
/// For this version, the location where each minima or maxima occurs (ix,iy,iz)
/// is represented by an integer (called an "index") which equals:
/// @code
///          index = ix + iy*image_size[0] + iz*image_size[0]*image_size[1]
/// @endcode
/// where ix,iy,iz are the coordinates of the corresponding voxel in the image,
/// and image_size[] stores the size of the image in the x,y,z directions.
/// If pv_minima_indices!=nullptr, then *pv_minima_indices will store a list
/// of the indices corresponding to the locations of the local minima.
/// If pv_maxima_indices!=nullptr, then *pv_maxima_indices will store a list
/// of the indices corresponding to the locations of the local maxima.
/// The corresponding voxel intensities (brightness values) will be stored in
/// *pv_minima_scores and *pv_maxima_scores (assuming they are != nullptr).
/// Thresholds can be used to discard minima or maxima whose corresponding
/// voxel intensities are not sufficiently low or high, respectively.
/// If the aaafMask[][][] is not equal to nullptr, then local minima and maxima
/// will be ignored if the corresponding entry in aaafMask[][][] equals 0.
///
/// @note: THIS VERSION OF THE FUNCTION WAS NOT INTENDED FOR PUBLIC USE.

template<typename Scalar, typename IntegerIndex, typename Label>
static void
_FindExtrema(int const image_size[3],          //!< size of the image in x,y,z directions
             Scalar const *const *const *aaafSource, //!< image array aaafSource[iz][iy][ix]
             Scalar const *const *const *aaafMask, //!< optional: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
             vector<IntegerIndex> *pv_minima_indices, //!< a list of integers uniquely identifying the location of each minima
             vector<IntegerIndex> *pv_maxima_indices, //!< a list of integers uniquely identifying the location of each maxima
             vector<Scalar> *pv_minima_scores, //!< store corresponding minima aaafSource[iz][iy][ix] values here (if not nullptr)
             vector<Scalar> *pv_maxima_scores, //!< store corresponding maxima aaafSource[iz][iy][ix] values here (if not nullptr)
             vector<IntegerIndex> *pv_minima_nvoxels, //!< store number of voxels in each minima (usually 1)
             vector<IntegerIndex> *pv_maxima_nvoxels, //!< store number of voxels in each maxima (usually 1)
             Scalar minima_threshold = std::numeric_limits<Scalar>::infinity(),  // Ignore minima which are not sufficiently low
             Scalar maxima_threshold = -std::numeric_limits<Scalar>::infinity(), // Ignore maxima which are not sufficiently high
             int connectivity=3,                      //!< square root of search radius around each voxel (1=nearest_neighbors, 2=diagonal2D, 3=diagonal3D)
             bool allow_borders=true,  //!< if true, plateaus that touch the image border (or mask boundary) are valid extrema
             Label ***aaaiDest=nullptr,  //!< optional: create an image showing where the extrema are?
             ostream *pReportProgress=nullptr)  //!< print progress to the user?
{
  assert(aaafSource);

  bool find_minima = (pv_minima_indices != nullptr);
  bool find_maxima = (pv_maxima_indices != nullptr);

  vector<IntegerIndex> minima_indices; //store minima indices here
  vector<IntegerIndex> maxima_indices; //store maxima indices here
  vector<Scalar>  minima_scores;  //store the brightness of the minima here
  vector<Scalar>  maxima_scores;  //store the brightness of the maxima here
  vector<IntegerIndex> minima_nvoxels; //store number of voxels in each minima here
  vector<IntegerIndex> maxima_nvoxels; //store number of voxels in each maxima here
  if (pv_minima_indices == nullptr)
    pv_minima_indices = &minima_indices;
  if (pv_maxima_indices == nullptr)
    pv_maxima_indices = &maxima_indices;
  if (pv_minima_scores == nullptr)
    pv_minima_scores = &minima_scores;
  if (pv_maxima_scores == nullptr)
    pv_maxima_scores = &maxima_scores;
  if (pv_minima_nvoxels == nullptr)
    pv_minima_nvoxels = &minima_nvoxels;
  if (pv_maxima_nvoxels == nullptr)
    pv_maxima_nvoxels = &maxima_nvoxels;


  ptrdiff_t ***aaaiExtrema; // 3-D array storing the plateau to which each voxel belongs (if any)
  ptrdiff_t *aiExtrema;     // the same array, stored contiguously in 1-D
  Alloc3D(image_size,
          &aiExtrema,
          &aaaiExtrema);

  // We will assign voxels in aaaiExtrema[][][] to the following values:
  ptrdiff_t NEITHER = 0;//(means this voxel is neither a local minima or maxima)
  ptrdiff_t UNDEFINED = image_size[0] * image_size[1] * image_size[2] + 1;
  ptrdiff_t QUEUED    = image_size[0] * image_size[1] * image_size[2] + 2;
  // ... or we will assign the voxel to an integer denoting which
  //     local minima or maxima it belongs to.
  //Note: Both "UNDEFINED" and "QUEUED" are impossible values (in this respect).
  //      There cannot be that many local maxima or minima in the image.


  for (int iz = 0; iz < image_size[2]; iz++)
    for (int iy = 0; iy < image_size[1]; iy++)
      for (int ix = 0; ix < image_size[0]; ix++)
        aaaiExtrema[iz][iy][ix] = UNDEFINED;


  if (pReportProgress)
    *pReportProgress
      << " -- Attempting to allocate space for one more image.       --\n"
      << " -- (If this crashes your computer, find a computer with   --\n"
      << " --  more RAM and use \"ulimit\", OR use a smaller image.)   --\n"
      << "\n";

  if (pReportProgress)
    *pReportProgress << "---- searching for local minima & maxima ----\n";

  // Figure out which neighbors to consider when searching neighboring voxels
  int (*neighbors)[3] = nullptr; //a pointer to a fixed-length array of 3 ints
  int num_neighbors = 0;
  {
    // How big is the search neighborhood around each minima?
    int r_neigh_sqd = connectivity;
    int r_neigh = floor(sqrt(r_neigh_sqd));
    
    vector<array<int, 3> > vNeighbors;
    for (int jz = -r_neigh; jz <= r_neigh; jz++) {
      for (int jy = -r_neigh; jy <= r_neigh; jy++) {
        for (int jx = -r_neigh; jx <= r_neigh; jx++) {
          array<int, 3> j_xyz;
          j_xyz[0] = jx;
          j_xyz[1] = jy;
          j_xyz[2] = jz;
          if ((jx == 0) && (jy == 0) && (jz == 0))
            continue;
          else if (jx*jx+jy*jy+jz*jz > r_neigh_sqd)
            continue;
          else
            vNeighbors.push_back(j_xyz);
        }
      }
    }
    // Convert the list of neighbors vNeigh to a contiguous 2D C-style array:
    num_neighbors = vNeighbors.size();
    neighbors = new int[num_neighbors][3];
    for (int j = 0; j < num_neighbors; j++)
      for (int d = 0; d < 3; d++)
        neighbors[j][d] = vNeighbors[j][d];
    // We will use neighbors[][] from now on..
  } // ...done figuring out the list of voxel neighbors


  // Now search the image for minima, maxima.
  // Do not assume that the local minima or maxima are point-like.
  // They could be "plateaus".
  // So at every voxel location ix,iy,iz, search for neighboring voxels
  // with identical intensities (brightness).  These voxels collectively are
  // a "plateau".  Find the neighbors of these voxels to determine whether
  // this "plateau" is also a local minima or maxima.
  for (int iz0 = 0; iz0 < image_size[2]; iz0++) {

    if (pReportProgress)
      *pReportProgress << "  searching for extrema at z="<<iz0+1<<" (of "<<image_size[2]<<")" << endl;

    for (int iy0 = 0; iy0 < image_size[1]; iy0++) {

      for (int ix0 = 0; ix0 < image_size[0]; ix0++) {

        if (aaafMask && aaafMask[iz0][iy0][ix0] == 0.0)
          continue;

        if (aaaiExtrema[iz0][iy0][ix0] != UNDEFINED)
          continue;

        bool is_minima = true;
        bool is_maxima = true;


        //#ifndef NDEBUG
        //set<array<int, 3> > visited;
        //#endif


        // It's possible that a local minima or maxima consists of more than one
        // adjacent voxels with identical intensities (brightnesses).
        // These voxels belong to a "plateau" of voxels of the same height
        // which are above (or below) all of the surrounding voxels.
        // (Note: Most of the time this plateau consists of only one voxel.)
        // The "q_plateau" variable is a data structure to keep track of which
        // voxels belong to the current plateau.
        queue<array<int, 3> > q_plateau;
        array<int, 3> i_xyz0;
        i_xyz0[0] = ix0;//note:this can be shortened to 1 line with make_array()
        i_xyz0[1] = iy0;
        i_xyz0[2] = iz0;
        q_plateau.push(i_xyz0);
        size_t n_plateau = 0; // the number of voxels in this plateau

        // We also need a reverse lookup table to check whether a given voxel
        // is in the queue (or has already been previously assigned).
        // That's what aaaiExtrema[][][] is.
        // If aaaiExtrema[iz][iy][ix] == QUEUED, then we know it's in the queue.
        // We need to come up with an impossible value for QUEUED:
        aaaiExtrema[iz0][iy0][ix0] = QUEUED; //make sure we don't visit it twice

        while (! q_plateau.empty())
        {
          array<int, 3> p = q_plateau.front();
          q_plateau.pop();
          int ix = p[0];
          int iy = p[1];
          int iz = p[2];
          n_plateau++;
          assert(aaaiExtrema[iz][iy][ix] == QUEUED);


          //#ifndef NDEBUG
          //assert(visited.find(p) == visited.end());
          //visited.insert(p);
          //#endif


          for (int j = 0; j < num_neighbors; j++) {
            int jx = neighbors[j][0];
            int jy = neighbors[j][1];
            int jz = neighbors[j][2];
            int iz_jz = iz + jz;
            int iy_jy = iy + jy;
            int ix_jx = ix + jx;
            if (((iz_jz < 0) || (image_size[2] <= iz_jz))
                ||
                ((iy_jy < 0) || (image_size[1] <= iy_jy))
                ||
                ((ix_jx < 0) || (image_size[0] <= ix_jx))
                ||
                (aaafMask && (aaafMask[iz_jz][iy_jy][ix_jx] == 0.0)))
            {
              if (! allow_borders) {
                is_minima = false;
                is_maxima = false;
              }
              continue;
            }



            if (aaafSource[iz_jz][iy_jy][ix_jx] == aaafSource[iz][iy][ix])
            {
              if (aaaiExtrema[iz_jz][iy_jy][ix_jx] == UNDEFINED)//don't visit twice
              {
                // then add this voxel to the q_plateau
                array<int, 3> ij_xyz;
                ij_xyz[0] = ix_jx;
                ij_xyz[1] = iy_jy;
                ij_xyz[2] = iz_jz;
                q_plateau.push(ij_xyz);
                aaaiExtrema[iz_jz][iy_jy][ix_jx] = QUEUED;

                //#ifndef NDEBUG
                //assert(visited.find(ij_xyz) == visited.end());
                //#endif
              }
            }
            else {
              if (aaafSource[iz+jz][iy+jy][ix+jx] < aaafSource[iz][iy][ix])
                is_minima = false;
              else if (aaafSource[iz+jz][iy+jy][ix+jx] > aaafSource[iz][iy][ix])
                is_maxima = false;
              else
                assert(false);
            }
          } //for (int j = 0; j < num_neighbors; j++)
        } // while (! q_plateau.empty())


        // If this voxel is either a minima or a maxima, add it to the list.
        if (is_minima && find_minima) {
          if (((! aaafMask) || (aaafMask[iz0][iy0][ix0] != 0)) &&
              ((aaafSource[iz0][iy0][ix0] <= minima_threshold)
               // || (maxima_threshold < minima_threshold)
               ))
          {
            // convert from a 3D location (ix,iy,iz) to a 1D index ("index")
            IntegerIndex index = ix0 + image_size[0]*(iy0 + image_size[1]*iz0);
            //minima_crds.push_back(ixiyiz);
            pv_minima_indices->push_back(index);
            pv_minima_scores->push_back(aaafSource[iz0][iy0][ix0]);
            pv_minima_nvoxels->push_back(n_plateau);
          }
        }
        if (is_maxima && find_maxima) {
          if (((! aaafMask) || (aaafMask[iz0][iy0][ix0] != 0)) &&
              ((aaafSource[iz0][iy0][ix0] >= maxima_threshold)
               // || (maxima_threshold < minima_threshold)
               ))
          {
            // convert from a 3D location (ix,iy,iz) to a 1D index ("index")
            IntegerIndex index = ix0 + image_size[0]*(iy0 + image_size[1]*iz0);
            //maxima_crds.push_back(ixiyiz);
            pv_maxima_indices->push_back(index);
            pv_maxima_scores->push_back(aaafSource[iz0][iy0][ix0]);
            pv_maxima_nvoxels->push_back(n_plateau);
          }
        }



        //#ifndef NDEBUG
        //visited.clear();
        //#endif



        //now loop through those same voxels again to reset the aaaiExtrema array

        // What to store in aaaiExtrema[][][] for voxels in this plateau?
        ptrdiff_t assigned_plateau_label;
        if (is_maxima)
          assigned_plateau_label = pv_maxima_scores->size();
        else if (is_minima)
          assigned_plateau_label = -pv_minima_scores->size();
        else
          assigned_plateau_label = NEITHER;
        aaaiExtrema[iz0][iy0][ix0] = assigned_plateau_label;


        assert(q_plateau.empty());
        q_plateau.push(i_xyz0);


        while (! q_plateau.empty())
        {
          array<int, 3> p = q_plateau.front();
          q_plateau.pop();
          int ix = p[0];
          int iy = p[1];
          int iz = p[2];
          assert(aaaiExtrema[iz][iy][ix] != QUEUED);

          //#ifndef NDEBUG
          //assert(visited.find(p) == visited.end());
          //visited.insert(p);
          //#endif


          for (int j = 0; j < num_neighbors; j++) {
            int jx = neighbors[j][0];
            int jy = neighbors[j][1];
            int jz = neighbors[j][2];
            int iz_jz = iz + jz;
            int iy_jy = iy + jy;
            int ix_jx = ix + jx;
            if (((iz_jz < 0) || (image_size[2] <= iz_jz))
                ||
                ((iy_jy < 0) || (image_size[1] <= iy_jy))
                ||
                ((ix_jx < 0) || (image_size[0] <= ix_jx))
                ||
                (aaafMask && (aaafMask[iz_jz][iy_jy][ix_jx] == 0.0)))
              continue;
            if (aaaiExtrema[iz_jz][iy_jy][ix_jx] == QUEUED) //don't visit twice
            {
              assert(aaafSource[iz_jz][iy_jy][ix_jx] == aaafSource[iz][iy][ix]);

              // then add this voxel to the q_plateau
              array<int, 3> ij_xyz;
              ij_xyz[0] = ix_jx;
              ij_xyz[1] = iy_jy;
              ij_xyz[2] = iz_jz;
              q_plateau.push(ij_xyz);

              // Commenting out:
              //assert(! (is_minima && is_maxima));   <- wrong
              // Actually this is possible.  For example in an image where all
              // the voxels are identical, all voxels are BOTH minima a maxima.
              // This is such a rare, pathelogical case, I don't worry about it.
              // (The aaaiExtrema[][][] array is only used to prevent visiting
              //  the same voxel twice.  It doesn't matter what is stored there
              //  as long as it fulfills this role.)

              aaaiExtrema[iz_jz][iy_jy][ix_jx] = assigned_plateau_label;
            }
          } //for (int j = 0; j < num_neighbors; j++)
        } // while (! q_plateau.empty())
      } // for (int ix0 = 0; ix0 < image_size[0]; ix0++)
    } // for (int iy0 = 0; iy0 < image_size[1]; iy0++)
  } // for (int iz0 = 0; iz0 < image_size[2]; iz0++)


  // Sort the minima and maxima in increasing and decreasing order, respectively
  IntegerIndex n_minima, n_maxima;
  bool descending_order;

  if (pv_minima_indices) {
    // Sort the minima in increasing order
    n_minima = pv_minima_indices->size();
    if (n_minima > 0) {
      descending_order = false;
      vector<tuple<Scalar, IntegerIndex> > score_index(n_minima);
      for (IntegerIndex i = 0; i < n_minima; i++)
        score_index[i] = make_tuple((*pv_minima_scores)[i], i);
      if (pReportProgress)
        *pReportProgress << "-- Sorting minima according to their scores... ";
      if (descending_order)
        sort(score_index.rbegin(),
             score_index.rend());
      else
        sort(score_index.begin(),
             score_index.end());
      vector<IntegerIndex> permutation(n_minima);
      for (IntegerIndex i = 0; i < score_index.size(); i++)
        permutation[i] = get<1>(score_index[i]);
      score_index.clear();
      apply_permutation(permutation, *pv_minima_indices);
      apply_permutation(permutation, *pv_minima_scores);
      apply_permutation(permutation, *pv_minima_nvoxels);
      // Optional: The minima in the image are not in sorted order either. Fix?
      vector<IntegerIndex> perm_inv;
      invert_permutation(permutation, perm_inv);
      for (int iz = 0; iz < image_size[2]; iz++) {
        for (int iy = 0; iy < image_size[1]; iy++) {
          for (int ix = 0; ix < image_size[0]; ix++) {
            if (aaafMask && (aaafMask[iz][iy][ix] == 0))
              continue;
            if (aaaiExtrema[iz][iy][ix] < 0)
              aaaiExtrema[iz][iy][ix]=-perm_inv[(-aaaiExtrema[iz][iy][ix])-1]-1;
          }
        }
      }
      if (pReportProgress)
        *pReportProgress << "done --" << endl;
    }
  }

  if (pv_maxima_indices) {
    // Sort the maxima in decreasing order
    n_maxima = pv_maxima_indices->size();
    if (n_maxima > 0) {
      descending_order = true;
      vector<tuple<Scalar, IntegerIndex> > score_index(n_maxima);
      for (IntegerIndex i = 0; i < n_maxima; i++)
        score_index[i] = make_tuple((*pv_maxima_scores)[i], i);
      if (pReportProgress)
        *pReportProgress << "-- Sorting maxima according to their scores... ";
      if (descending_order)
        sort(score_index.rbegin(),
             score_index.rend());
      else
        sort(score_index.begin(),
             score_index.end());
      vector<IntegerIndex> permutation(n_maxima);
      for (IntegerIndex i = 0; i < score_index.size(); i++)
        permutation[i] = get<1>(score_index[i]);
      score_index.clear();
      apply_permutation(permutation, *pv_maxima_indices);
      apply_permutation(permutation, *pv_maxima_scores);
      apply_permutation(permutation, *pv_maxima_nvoxels);
      if (pReportProgress)
        *pReportProgress << "done --" << endl;
      // Optional: The maxima in the image are not in sorted order either. Fix?
      vector<IntegerIndex> perm_inv;
      invert_permutation(permutation, perm_inv);
      for (int iz = 0; iz < image_size[2]; iz++) {
        for (int iy = 0; iy < image_size[1]; iy++) {
          for (int ix = 0; ix < image_size[0]; ix++) {
            if (aaafMask && (aaafMask[iz][iy][ix] == 0))
              continue;
            if (aaaiExtrema[iz][iy][ix] > 0)
              aaaiExtrema[iz][iy][ix] = perm_inv[aaaiExtrema[iz][iy][ix]-1]+1;
          }
        }
      }
    }
  }

  if (aaaiDest) {
    for (int iz = 0; iz < image_size[2]; iz++) {
      for (int iy = 0; iy < image_size[1]; iy++) {
        for (int ix = 0; ix < image_size[0]; ix++) {
          if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
            continue; // don't modify voxels outside the mask
          aaaiDest[iz][iy][ix] = aaaiExtrema[iz][iy][ix];
          // Until now, minima in the image are represented by negative integers
          // and maxima are represented by positive integers.
          if (((! find_minima) || (! find_maxima)) &&
              (aaaiExtrema[iz][iy][ix] < 0))
            // If we only asked for minima OR maxima (not both)
            // then just use positive integers only.
            aaaiExtrema[iz][iy][ix] = -aaaiExtrema[iz][iy][ix];
        }
      }
    }
  }

  Dealloc3D(image_size,
            &aiExtrema,
            &aaaiExtrema);

  delete [] neighbors;

} // _FindExtrema()




/// @brief  The following function finds EITHER local minima OR local maxima
///         in an image.
///
/// For this version, the location where each minima or maxima occurs (ix,iy,iz)
/// is represented by an integer (called an "index") which equals:
/// @code
///          index = ix + iy*image_size[0] + iz*image_size[0]*image_size[1]
/// @endcode
/// where ix,iy,iz are the coordinates of the corresponding voxel in the image.
///
/// @note: THIS VERSION OF THE FUNCTION WAS NOT INTENDED FOR PUBLIC USE.

template<typename Scalar, typename IntegerIndex, typename Label>
static void
_FindExtrema(int const image_size[3],          //!< size of the image in x,y,z directions
             Scalar const *const *const *aaafI,    //!< image array aaafI[iz][iy][ix]
             Scalar const *const *const *aaafMask, //!< optional: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
             vector<IntegerIndex> &extrema_indices, //!< a list of integers uniquely identifying the location of each minima or maxima
             vector<Scalar> &extrema_scores, //!< corresponding voxel brightness at that location
             vector<IntegerIndex> &extrema_nvoxels, //!< how many voxels belong to each minima or maxima?
             bool seek_minima=true,    //!< search for minima or maxima?
             Scalar threshold=std::numeric_limits<Scalar>::infinity(), // Ignore minima or maxima which are not sufficiently low or high
             int connectivity=3,       //!< square root of search radius around each voxel (1=nearest_neighbors, 2=diagonal2D, 3=diagonal3D)
             bool allow_borders=true,  //!< if true, plateaus that touch the image border (or mask boundary) are valid extrema
             Label ***aaaiDest=nullptr,  //!< optional: create an image showing where the extrema are?
             ostream *pReportProgress=nullptr)  //!< print progress to the user?
{
  // NOTE:
  // C++ will not allow us to supply nullptr to a function that expects a pointer 
  // to a template expression: Template argument deduction/substitution fails.
  // We need to re-cast "nullptr" as a pointer with the correct type.
  // One way to do that is to define these new versions of nullptr:
  vector<IntegerIndex> *null_vI = nullptr;  
  vector<IntegerIndex> *null_vi = nullptr;  
  vector<Scalar> *null_vf = nullptr;  
  if (seek_minima) {
    _FindExtrema(image_size,
                 aaafI,
                 aaafMask,
                 &extrema_indices, // store maxima locations here
                 null_vI, // <-- don't search for maxima
                 &extrema_scores, // store minima values here
                 null_vf, // <-- don't search for maxima
                 &extrema_nvoxels, // store number of voxels in each minima
                 null_vi, // <-- don't search for maxima
                 threshold,
                 -std::numeric_limits<Scalar>::infinity(),
                 connectivity,
                 allow_borders,
                 aaaiDest,
                 pReportProgress);
  }
  else {
    if (threshold == std::numeric_limits<Scalar>::infinity())
      threshold = -std::numeric_limits<Scalar>::infinity();
    _FindExtrema(image_size,
                 aaafI,
                 aaafMask,
                 null_vI, // <-- don't search for minima_crds,
                 &extrema_indices, // store maxima locations here
                 null_vf, // <-- don't search for minima_scores,
                 &extrema_scores, // store maxima values here
                 null_vi, // <-- don't search for minima
                 &extrema_nvoxels, // store number of voxels in each maxima
                 std::numeric_limits<Scalar>::infinity(),
                 threshold,
                 connectivity,
                 allow_borders,
                 aaaiDest,
                 pReportProgress);
  }
  
} //_FindExtrema()



} //namespace visfd



#endif //#ifndef _DO_NOT_USE_FEATURE_HPP
