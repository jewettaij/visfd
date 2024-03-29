///   @file segmentation.hpp
///   @brief a collection of functions relevant for image segmentation
///   @author Andrew Jewett
///   @date 2019-4-15

#ifndef _SEGMENTATION_HPP
#define _SEGMENTATION_HPP

#include <cstring>
#include <cassert>
#include <cmath>
#include <limits>
#include <ostream>
#include <vector>
#include <map>
#include <tuple>
#include <set>
#include <queue>
using namespace std;
#include <err_visfd.hpp> // defines the "VisfdErr" exception type
#include <eigen3_simple.hpp>  // defines matrix diagonalizer (DiagonalizeSym3())
#include <lin3_utils.hpp>     // defines "Normalize3()"
#include <visfd_utils.hpp>    // defines invert_permutation(), AveArray(), ...
#include <alloc2d.hpp>    // defines Alloc2D() and Dealloc2D()
#include <alloc3d.hpp>    // defines Alloc3D() and Dealloc3D()
#include <filter1d.hpp>   // defines "Filter1D" (used in ApplySeparable())
#include <filter3d.hpp>   // defines common 3D image filters
#include <morphology.hpp> // defines FindExtrema()


namespace visfd {


/// @brief  Perform the Meyer's flood-fill algorithm to segment a 3D image.
///         After this function has completed, every voxel in aaaiDest[][][]
///         will be assigned an integer which is one of 4 types:
///      1) an integer between 1 and N_BASINS (the number of local minima in
///         the image, indicating the basin to which the voxel belongs)
///      2) label_boundary (which is 0 by default)
///         (if show_boundaries == true, and if
///          the voxel is located at the boundary between 2 or more basins)
///      3) label_undefined (which is -1 by default)
///         (when the voxel's intensity passes the threshold set by the caller)
///      4) left unmodified (if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0)
///
///      If aaafMarkers != nullptr, then the voxels in this image (if > 0)
///      are assumed to correspond to the labels that the caller wants
///      these voxels to have.
///
/// @returns  The number of watershed basins found.
/// @note  Users can specify the number, locations, and labels of the watershed
///        basins in the aaaiMarkers[][][] array.  This array is an image whose
///        voxel brightnesses should be filled with nonnegative (>=0) integers.
///        Entries in the array >= 1 are the initial watershed basins and their
///        labels.  (Entries in the array <= 0 will be ignored.)
///        The number of distinct basins found by this function will equal
///        the number of distinct entries found in the aaaiMarkers[][][] array.
///        If this array is unspecified (ie. if aaaiMarkers == nullptr),
///        then the local minima or maxima will be used
///        as the initial watershed basins by default,
///        (depending on the "start_from_minima" argument).
/// @note  If the "halt_threshold" argument is not supplied by the caller,
///        then std::numeric_limits::infinity is used by default.

template<typename Scalar, typename Label, typename Coordinate>

size_t
Watershed(int const image_size[3],                   //!< #voxels in xyz
          Scalar const *const *const *aaafSource,    //!< intensity of each voxel
          Label ***aaaiDest,                         //!< watershed segmentation results will go here.
          Scalar const *const *const *aaafMask,      //!< OPTIONAL: Ignore voxels whose mask value is 0
          Label const *const *const *aaaiMarkers,    //!< OPTIONAL: If not NULL, users can specify locations and labels of basins. Entries <= 0 are ignored.
          Scalar halt_threshold=std::numeric_limits<Scalar>::infinity(), //!< OPTIONAL: voxels with brightness exceeding this threshold are assigned "label_undefined"
          bool start_from_minima=true,               //!< OPTIONAL: start from local minima? (if false, maxima will be used)
          int connectivity=1,                        //!< OPTIONAL: square root of the search radius around each voxel (1=nearest_neighbors, 2=2D_diagonal, 3=3D_diagonal)
          bool show_boundaries=true,                 //!< OPTIONAL: should voxels on the boundary between two basins be labelled differently?
          Label label_boundary=0,                   //!< OPTIONAL: if so, what label should boundary voxels be assigned to?
          Label label_undefined=-1,                 //!< OPTIONAL: assign this label to masked voxels or to voxels with brightness exceeding the halt_threshold
          vector<array<Coordinate, 3> > *pv_basin_locations=nullptr, //!< OPTIONAL: store the location of each minima or maxima
          vector<Scalar> *pv_basin_scores=nullptr, //!< OPTIONAL: store the voxel intensities (brightnesses) at these locations here
          ostream *pReportProgress=nullptr           //!< OPTIONAL: print progress to the user?
          )
{
  assert(image_size);
  assert(aaafSource);
  assert(aaaiDest);

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

  assert(neighbors);

  Scalar SIGN_FACTOR = 1.0;
  if (! start_from_minima) {
    SIGN_FACTOR = -1.0;
    if (halt_threshold == std::numeric_limits<Scalar>::infinity())
      // Then the caller did not bother to specify a threshold.  Note: the
      // default threshold assumes we are searching for minima, not maxima.
      // If the caller did not specify a halting threshold, AND if we ARE
      // looking for maxima (not minima), then we have to reverse the sign of
      // the default threshold as well (this is a little confusing to me too):
      halt_threshold = -std::numeric_limits<Scalar>::infinity();
  }

  vector<array<Coordinate, 3> > basin_locations;
  if (pv_basin_locations == nullptr)
    pv_basin_locations = &basin_locations;

  vector<Scalar> basin_scores;
  if (pv_basin_scores == nullptr)
    pv_basin_scores = &basin_scores;

  vector<size_t> basin_nvoxels;
  vector<size_t> *pv_basin_nvoxels = &basin_nvoxels;


  ptrdiff_t WATERSHED_BOUNDARY = 0; //an impossible value
  ptrdiff_t UNDEFINED = -1; //an impossible value
  size_t num_basins; // keep track of the number of watershed basins found
  ptrdiff_t max_label = 0;

  // Did the caller already specify the labels they want us to use?
  // If so, start with these labels (instead of the minima or maxima).
  if (aaaiMarkers) {
    num_basins = 0;
    set<ptrdiff_t> labels_so_far;
    // NOT NEEDED (pv_basin_nvoxels), so commenting out:
    // map<ptrdiff_t, size_t> labels_2_basin_id;
    for (int iz=0; iz<image_size[2]; iz++) {
      for (int iy=0; iy<image_size[1]; iy++) {
        for (int ix=0; ix<image_size[0]; ix++) {
          if (aaafMask && aaafMask[iz][iy][ix] == 0.0)
            continue;
          ptrdiff_t label = aaaiMarkers[iz][iy][ix];
          if (label > 0) {
            //aaaiDest[iz][iy][ix] = label;
            if (labels_so_far.find(label) == labels_so_far.end()) {
              num_basins++;
              labels_so_far.insert(label);
              max_label = std::max(max_label, label);
              if (pv_basin_locations) {
                array<Coordinate,3> ixyz;
                ixyz[0] = ix;
                ixyz[1] = iy;
                ixyz[2] = iz;
                (*pv_basin_locations).push_back(ixyz);
              }
              if (pv_basin_scores)
                (*pv_basin_scores).push_back(aaafSource[iz][iy][ix]);
              // NOT NEEDED (pv_basin_nvoxels), so commenting out:
              //if (pv_basin_nvoxels) {
              //  labels_2_basin_id[label] = (*pv_basin_nvoxels).size();
              //  (*pv_basin_nvoxels).push_back(1);
              //}
            } //if (labels_so_far.find(label) != labels_so_far.end()) {
            // NOT NEEDED (pv_basin_nvoxels), so commenting out:
            //else if (pv_basin_nvoxels) {
            //  assert(labels_2_basin_id.find(label)!=labels_2_basin_id.end());
            //  size_t label_id = labels_2_basin_id[label];
            //  pv_basin_nvoxels[label_id]++;
            //}
          } // if (label > 0)
          //else
          //  aaaiDest[iz][iy][ix] = UNDEFINED;
        } // for (int ix=0; ix<image_size[0]; ix++)
      } // for (int iy=0; iy<image_size[1]; iy++)
    } // for (int iz=0; iz<image_size[2]; iz++)
  } // if (aaaiMarkers)

  else // if (! aaaiMarkers)
  {
    // Find all the local minima (or maxima?) in the image.
    num_basins =
      _FindExtrema(image_size,
                   aaafSource,
                   aaafMask,
                   *pv_basin_locations,
                   *pv_basin_scores,
                   *pv_basin_nvoxels,
                   start_from_minima, //<-- minima or maxima?
                   halt_threshold,
                   connectivity,
                   true,
                   static_cast<Label***>(nullptr),
                   pReportProgress);
    max_label = num_basins;
  } // if (! aaaiMarkers)

  ptrdiff_t QUEUED = max_label + 1; //an impossible value
  //OLD VERSION: ptrdiff_t UNDEFINED = max_label + 1; //an impossible value
  //OLD VERSION: ptrdiff_t QUEUED = max_label + 2; //an impossible value

  //initialize aaaiDest[][][]
  for (int iz=0; iz<image_size[2]; iz++) {
    for (int iy=0; iy<image_size[1]; iy++) {
      for (int ix=0; ix<image_size[0]; ix++) {
        //if (aaafMask && aaafMask[iz][iy][ix] == 0.0)
        //  continue;
        aaaiDest[iz][iy][ix] = UNDEFINED;
      }
    }
  }

  // Define "q", the set of voxels which are adjacent to the voxels we have
  // processed so far.  These are the voxels to be processed in the next
  // iteration.  This is implemented as a priority-queue, sorted by intensity.
  // Each voxel in "q" has an intensity, a basin-ID (Label), and a location.

  priority_queue<tuple
                  <
                   Scalar,    // the "height" of that voxel (brightness)
                   ptrdiff_t,     // the basin-ID to which the voxel belongs
                   array<Coordinate, 3>  // location of the voxel
                  >
                > q;
               

  if (pReportProgress)
    *pReportProgress <<
      " ---- Watershed segmentation algorithm ----\n"
      "starting from " << num_basins << " different local "
                     << (start_from_minima ? "minima" : "maxima") << endl;

  // Initialize the queue so that it only contains voxels which lie
  // at the location of a local minima (or maxima).
  // Also keep track of which minima (which "basin") they belong to.

  for (size_t i=0; i < num_basins; i++) {
    // Create an entry in q for each of the local minima (or maxima)

    // Assign a different integer to each of these minima, starting at 1
    ptrdiff_t which_basin = i;

    int ix = (*pv_basin_locations)[i][0];
    int iy = (*pv_basin_locations)[i][1];
    int iz = (*pv_basin_locations)[i][2];

    // These entries in the priority queue will be sorted by "score"
    // which is the intensity of the image at this voxel location.
    Scalar score = (*pv_basin_scores)[i];
    assert(score == aaafSource[iz][iy][ix]);
    score *= SIGN_FACTOR; //(enable search for either local minima OR maxima)

    // Note:FindExtrema() should avoid minima above the halt_threshold,
    //     or maxima below the halt_threshold. We check for that with an assert:
    assert(score <= halt_threshold * SIGN_FACTOR);

    // copy the ix,iy,iz coordinates into an array<Coordinate, 3>
    array<Coordinate, 3> icrds;
    // (It seems like there should be a way to do this in 1 line, but I'm
    //  unfamiliar with the new C++ array initializer syntax and I'm on a plane 
    //  without internet at the moment.  So I'll just do it the obvious way.)
    icrds[0] = ix;
    icrds[1] = iy;
    icrds[2] = iz;

    q.push(make_tuple(-score, // <-- entries sorted lexicographically by -score
                      which_basin,
                      icrds));

    assert(aaaiDest[iz][iy][ix] == UNDEFINED);
    aaaiDest[iz][iy][ix] = QUEUED;

  } // for (size_t i=0; i < num_basins; i++)


  // Count the number of voxels in the image we will need to consider.
  // (This will be used for printing out progress estimates later.)
  size_t n_voxels_processed = 0;
  //size_t n_voxels_image = image_size[0] * image_size[1] * image_size[2];
  size_t n_voxels_image = 0;
  for (int iz=0; iz<image_size[2]; iz++) {
    for (int iy=0; iy<image_size[1]; iy++) {
      for (int ix=0; ix<image_size[0]; ix++) {
        if (aaafMask && aaafMask[iz][iy][ix] == 0.0)
          continue;
        if (aaafSource[iz][iy][ix]*SIGN_FACTOR > halt_threshold*SIGN_FACTOR)
          continue;
        n_voxels_image++;
      }
    }
  }


  // Loop over all the voxels on the periphery
  // of the voxels with lower intensity (brightness).
  // This is a priority-queue of voxels which lie adjacent to 
  // the voxels which have been already considered.
  // As we consider new voxels, add their neighbors to this queue as well
  // until we process all voxels in the image.

  //#ifndef NDEBUG
  //set<array<Coordinate, 3> > visited;
  //#endif

  while (! q.empty())
  {
    tuple<Scalar, ptrdiff_t, array<Coordinate, 3> > p = q.top();
    q.pop();

    // Figure out the properties of that voxel (position, intensity/score,basin)
    Scalar i_score = -std::get<0>(p); //abs(i_score) = voxel intensity = aaafSource[iz][iy][ix]
    Scalar i_which_basin = std::get<1>(p); // the basin to which that voxel belongs (tentatively)
    int ix = std::get<2>(p)[0]; // voxel location
    int iy = std::get<2>(p)[1]; //   "      "
    int iz = std::get<2>(p)[2]; //   "      "

    // Should we ignore this voxel?

    if (i_score > halt_threshold * SIGN_FACTOR) {
      // stop when the voxel brightness(*SIGN_FACTOR) exceeds the halt_threshold
      aaaiDest[iz][iy][ix] = UNDEFINED;
      continue;
    }

    if (aaafMask && aaafMask[iz][iy][ix] == 0.0) {
      // ignore voxel if the user specified a mask and the voxel does not belong
      aaaiDest[iz][iy][ix] = UNDEFINED;
      continue;
    }

    assert(aaaiDest[iz][iy][ix] == QUEUED);
    // Now we assign this voxel to the basin
    aaaiDest[iz][iy][ix] = i_which_basin + 1;
    // (Note: This will prevent the voxel from being visited again.)


    if (pReportProgress) {
      // Every time the amount of progress increases by 1%, inform the user:
      n_voxels_processed++;
      size_t percentage = (n_voxels_processed*100) / n_voxels_image;
      size_t percentage_previous = ((n_voxels_processed-1)*100)/n_voxels_image;
      if (percentage != percentage_previous)
        *pReportProgress << " percent complete: " << percentage << endl;
    }


    //#ifndef NDEBUG
    //array<Coordinate, 3> ixiyiz;
    //ixiyiz[0] = ix;
    //ixiyiz[1] = iy;
    //ixiyiz[2] = iz;
    //assert(visited.find(ixiyiz) == visited.end());
    //visited.insert(ixiyiz);
    //#endif //#ifndef NDEBUG



    // ---------- Meyer (and Beucher's?) inter-pixel flood algorith: ---------
    // Find all the local minima (or maxima) in the image.
    // Each one of these minima (or maxima) is the center (seed) of a "basin"
    // (a group of voxels that surround that minima or maxima).
    // Initially these basins contain only the
    // voxel(s) which are located at that local minima (or maxima).
    // so we assign each of the voxels to the basin to which it belongs.
    // Then we add nearby voxels (adjacent) voxels to each basin until
    // a collision between basins occurs.
    //
    // Details:  Iterate over every voxel in the priority queue:
    //
    // Add each of the local minima (or maxima) voxels to a priority queue
    // of voxels to be considered later.
    //
    // 1) Retrieve the next lowest (or highest) intensity voxel from this queue.
    //    By definition, it will be adjacent to at least one basin.
    //    (The basin to which the voxel belonged that added it to the queue.)
    // 2) If that voxel is adjacent to only one basin assign it to that basin.
    //    If it is adjacent to multiple basins, assign it to WATERSHED_BOUNDARY.
    // 3) Search over that voxel's neighbors.  Add any neighboring voxels
    //    which have not yet been assigned to a basin already to the queue.
    // 4) Iterate until there are no more voxels left in the queue.


    // ---- loop over neighbors ----
    // Let (jx,jy,jz) denote the index offset for the neighbor voxels
    // (located at ix+jx,iy+jy,iz+jz)

    for (int j = 0; j < num_neighbors; j++) {
      int jx = neighbors[j][0];
      int jy = neighbors[j][1];
      int jz = neighbors[j][2];
      int iz_jz = iz + jz;
      int iy_jy = iy + jy;
      int ix_jx = ix + jx;
      if ((iz_jz < 0) || (image_size[2] <= iz_jz))
        continue;
      if ((iy_jy < 0) || (image_size[1] <= iy_jy))
        continue;
      if ((ix_jx < 0) || (image_size[0] <= ix_jx))
        continue;

      if (aaafMask && (aaafMask[iz_jz][iy_jy][ix_jx] == 0.0))
        continue;

      if (aaaiDest[iz_jz][iy_jy][ix_jx] == WATERSHED_BOUNDARY)
        continue;

      else if (aaaiDest[iz_jz][iy_jy][ix_jx] == QUEUED)
        continue;

      else if (aaaiDest[iz_jz][iy_jy][ix_jx] == UNDEFINED)
      {

        aaaiDest[iz_jz][iy_jy][ix_jx] = QUEUED;

        // and push this neighboring voxel onto the queue.
        // (...if they have not been assigned to a basin yet.  This
        //  insures that the same voxel is never pushed more than once.)
        array<Coordinate, 3> neighbor_crds;
        neighbor_crds[0] = ix_jx;
        neighbor_crds[1] = iy_jy;
        neighbor_crds[2] = iz_jz;
        Scalar neighbor_score = aaafSource[iz_jz][iy_jy][ix_jx]*SIGN_FACTOR;
        q.push(make_tuple(-neighbor_score,
                          i_which_basin,
                          neighbor_crds));
      }
      else {
        if (aaaiDest[iz_jz][iy_jy][ix_jx] != aaaiDest[iz][iy][ix])
        {
          assert((aaaiDest[iz][iy][ix] == i_which_basin + 1) ||
                 (aaaiDest[iz][iy][ix] == WATERSHED_BOUNDARY));

          if (show_boundaries)
          {
            // If these two neighboring voxels already belong to different 
            // basins, then we should assign the more SHALLOW voxel to 
            // WATERSHED_BOUNDARY. Which of the two voxels is "more shallow"?
            // This is the voxel with the higher intensity (*SIGN_FACTOR).
            assert(SIGN_FACTOR * aaafSource[iz_jz][iy_jy][ix_jx]  <=
                   SIGN_FACTOR * aaafSource[iz][iy][ix]);
            aaaiDest[iz][iy][ix] = WATERSHED_BOUNDARY;
          }
        } // if (aaaiDest[iz_jz][iy_jy][ix_jx] != aaaiDest[iz][iy][ix])
      } // else clause for "if (aaaiDest[iz_jz][iy_jy][ix_jx] == UNDEFINED)"
    } // for (int j = 0; j < num_neighbors; j++)
  } //while (q.size() != 0)

  #ifndef NDEBUG
  // DEBUGGING
  // All of the voxels should either be assigned to a basin,
  // OR assigned to "WATERSHED_BOUNDARY" or to "UNDEFINED"
  // (but not "QUEUED").  Check for that below:
  for (int iz=0; iz<image_size[2]; iz++) {
    for (int iy=0; iy<image_size[1]; iy++) {
      for (int ix=0; ix<image_size[0]; ix++) {
        if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
          continue;
        assert(aaaiDest[iz][iy][ix] != QUEUED);
      }
    }
  }
  #endif //#ifndef NDEBUG

  if (label_boundary != WATERSHED_BOUNDARY) {
    // If the caller explicitly specified the brightness that voxels 
    // on the boundaries between basins should be assigned to,
    // then take care of that now.
    for (int iz=0; iz<image_size[2]; iz++) {
      for (int iy=0; iy<image_size[1]; iy++) {
        for (int ix=0; ix<image_size[0]; ix++) {
          if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
            continue;
          if (aaaiDest[iz][iy][ix] == WATERSHED_BOUNDARY)
            aaaiDest[iz][iy][ix] = label_boundary;
        }
      }
    }
  }

  if (label_undefined != UNDEFINED) {
    // If the caller explicitly specified the brightness that "undefind" voxels
    // (ie. voxels whose brightnesses are disqualified for being too bright
    //      or not bright enough)
    // ...should have, then take care of that now.
    for (int iz=0; iz<image_size[2]; iz++) {
      for (int iy=0; iy<image_size[1]; iy++) {
        for (int ix=0; ix<image_size[0]; ix++) {
          if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
            continue;
          if (aaaiDest[iz][iy][ix] == UNDEFINED)
            aaaiDest[iz][iy][ix] = label_undefined;
        }
      }
    }
  }

  if (aaaiMarkers) {
    map<ptrdiff_t, ptrdiff_t> labels_old2new;
    for (int iz=0; iz<image_size[2]; iz++) {
      for (int iy=0; iy<image_size[1]; iy++) {
        for (int ix=0; ix<image_size[0]; ix++) {
          ptrdiff_t label_old = aaaiDest[iz][iy][ix];
          ptrdiff_t label_new = aaaiMarkers[iz][iy][ix];
          if ((label_new > 0) &&
              (label_old != label_boundary) &&
              (label_old != label_undefined))
            labels_old2new[label_old] = label_new;
        }
      }
    }
    for (int iz=0; iz<image_size[2]; iz++) {
      for (int iy=0; iy<image_size[1]; iy++) {
        for (int ix=0; ix<image_size[0]; ix++) {
          if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
            continue;
          if ((aaaiDest[iz][iy][ix] == label_boundary) ||
              (aaaiDest[iz][iy][ix] == label_undefined))
            continue;
          if (labels_old2new.find(aaaiDest[iz][iy][ix])!=labels_old2new.end()) {
            aaaiDest[iz][iy][ix] = labels_old2new[aaaiDest[iz][iy][ix]];
          }
          else
            aaaiDest[iz][iy][ix] = label_undefined;
        }
      }
    }
  } //if (aaaiMarkers)

  if (pReportProgress)
    *pReportProgress << "Number of basins found: "
                     << num_basins << endl;

  delete [] neighbors;

  return num_basins;

} //Watershed()






} //namespace visfd



#endif //#ifndef _SEGMENTATION_HPP
