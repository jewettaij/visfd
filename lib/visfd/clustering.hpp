///   @file segmentation.hpp
///   @brief a collection of functions useful for clustering
///          (grouping similar voxels together)
///   @author Andrew Jewett
///   @date 2019-4-15


#ifndef _CLUSTERING_HPP
#define _CLUSTERING_HPP


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
#include <lin3_utils.hpp>     // defines "Normalize3()"
#include <visfd_utils.hpp>    // defines invert_permutation(), AveArray(), ...
#include <alloc2d.hpp>    // defines Alloc2D() and Dealloc2D()
#include <alloc3d.hpp>    // defines Alloc3D() and Dealloc3D()
#include <filter1d.hpp>   // defines "Filter1D" (used in ApplySeparable())
#include <filter3d.hpp>   // defines common 3D image filters




namespace visfd {



/// @brief   The following enumerated constants are used as arguments
///          to the ClusterConnected() function.
///          You can use these options to control the order
///          the order that the clusters are returned to the caller.
///
///          SORT_BY_SIZE will sort clusters by their size.
///
///          SORT_BY_SALIENCY will sort clusters by their brightest voxel.
///          (The voxel with the highest saliency.)

typedef enum eClusterSortCriteria {
  SORT_BY_SALIENCY,
  SORT_BY_SIZE
} ClusterSortCriteria;



/// @brief  WARNING: EXPERIMENTAL CODE.
///        THIS FUNCTION'S ARGUMENT LIST MAY CHANGE SIGNIFICANTLY IN THE FUTURE.
///         This function is used to cluster voxels of high saliency
///         into islands which are connected together by adjacent voxels,
///         and are separated by regions of low saliency
///         (or, if aaaafVector or aaaafSymmetricTensor are not nullptr,
///          regions where the voxels point in incompatible directions
///          from their neighbors).
///         These different islands may correspond to different objects
///         in the source image.
///         Based on the watershed algorithm, this algorithm starts with
///         voxels which lie near a maxima (or a minima) of saliency, 
///         and groups voxels of similar saliency together (as well as voxels
///         whose vector and/or tensor directions are compatible, if applicable)
///
/// @return The function does not have a return value.
///         After the function is finished, the aaaiDest[][][] array will
///         store the cluster-ID associated with each voxel.
///         The clusters are sorted by their size (by default), and this is
///         relfected in their cluster-IDs.  (The largest cluster has ID 1)
///         (Cluster-ID numbers begin at 1, not 0.)
///         (The "pv_cluster_maxima" argument, if != nullptr, will store
///          the location of the saliency maxima (or minima) for each cluster.)
///
/// @note   Voxels below the saliency threshold are ignored, and will
///         assigned to the value of UNDEFINED, which (by default)
///         is equal to the number of clusters detected in the image plus 1.
///         (This value can be overridden using
///          the "label_undefined", and "undefined_label_is_max" arguments.)
///
/// @note   If a aaafMask array is supplied by the caller, then voxels located
///         in regions where aaafMask[iz][iy][ix]=0 will be ignored.
///
/// The remaining notes below describe the behavior of the optional
/// aaaafVector, aaaafSymmetricTensor, aaaafVectorStandardized arguments,
/// which are not generally useful in all situations
/// @note   If aaaafVector != nullptr, then
///         neighboring voxels whose direction (vector) changes discontinuously,
///         will not be added to an existing cluster.
/// @note   If aaaafSymmetricTensor != nullptr, then
///         neighboring voxels whose direction (tensor) changes discontinuously,
///         will not be added to an existing cluster.
/// @note   In addition, voxels whose Hessian (2nd derivative matrix) of the
///         aaafSaliency array is INCOMPATIBLE with the corresponding entry from
///         aaaafSymmetricTensor (if != nullptr) will not be added to any cluster.
///         (Tensors are compared using the normalized Trace-product.)
/// @note   In addition, voxels whose principle eigenvector from the Hessian
///         (2nd derivative matrix) of the aaafSaliency array is POINTING IN A
///         SIGNIFICANTLY DIFFERENT DIRECTION FROM the corresponding entry from
///         the aaaafVector array (if != nullptr) will not be added to any cluster.
/// @note   If (consider_dot_product_sign == false), the difference in vector
///         directions will be calculated by considering only the magnitude
///         of the dot product between them (not the sign).
///         (This is useful because some image processing algorithms generate
///         vectors which lack polarity.  For example, planar surface "ridge"
///         detectors calculate vectors normal to a surface, which are equally
///         likely to point inward or outward.)
/// @note   If (consider_dot_product_sign == false) AND
///         if aaaafVector and aaaafVectorStandardized are both non-null, THEN
///         aaaafVectorStandardized array will store a version of aaaafVector
///         array whose signs have been flipped to preserve consistency
///         of directionality as much as possible.  If any clusters contain any
///         closed loops loops that have directions which cannot be chosen
///         consistently, they will be cut at a location where their saliency
///         is weak.  (The goal is to avoid Möbius-strip-like defects.)

template<typename Scalar, typename Label, typename Coordinate, typename VectorContainer=Scalar*, typename TensorContainer=Scalar*>

void
ClusterConnected(int const image_size[3],                   //!< #voxels in xyz
                 Scalar const *const *const *aaafSaliency,  //!< intensity of each voxel
                 Label ***aaaiDest,                       //!< watershed segmentation results go here
                 Scalar const *const *const *aaafMask,    //!< optional: Ignore voxels whose mask value is 0
                 Scalar threshold_saliency=-std::numeric_limits<Scalar>::infinity(), //!< don't consider voxels with saliencies below this value
                 Label label_undefined = 0,               //!< voxels storing this value do not belong to any clusters
                 bool undefined_label_is_max = false,     //!< set label_undefined to number of clusters + 1 (overrides label_undefined)
                 VectorContainer const *const *const *aaaafVector=nullptr,
                 Scalar threshold_vector_saliency=M_SQRT1_2,    //!< voxels with incompatible saliency and vector are ignored (=-1.001 disables)
                 Scalar threshold_vector_neighbor=M_SQRT1_2,    //!< neighboring voxels with incompatible vectors are ignored (=-1.001 disables)
                 bool consider_dot_product_sign = true,        //!< does the sign of the dot product matter?  If not, compare abs(DotProduct()) with threshold_vector variables
                 TensorContainer const *const *const *aaaafSymmetricTensor=nullptr,
                 Scalar threshold_tensor_saliency=M_SQRT1_2,    //!< voxels with incompatible saliency and tensor are ignored (=-1.1 disables)
                 Scalar threshold_tensor_neighbor=M_SQRT1_2,    //!< neighboring voxels with incompatible tensors are ignored (=-1.1 disables)
                 bool tensor_is_positive_definite_near_ridge=true, //!< what is the sign of the principal tensor eigenvalue(s) near a ridge we care about?
                 int connectivity=1,                      //!< square root of the search radius around each voxel (1=nearest_neighbors, 2=2D_diagonal, 3=3D_diagonal)
                 vector<array<Coordinate, 3> > *pv_cluster_maxima=nullptr, //!< optional: the location of saliency minima or maxima which seeded each cluster
                 vector<Scalar> *pv_cluster_sizes=nullptr, //!< optional: what was the size of each cluster? (either the number of voxels, or the sum of voxel weights)
                 vector<Scalar> *pv_cluster_saliencies=nullptr, //!< optional: what was the saliency (brightness) of each cluster's brightest voxel?
                 ClusterSortCriteria sort_criteria = ClusterSortCriteria::SORT_BY_SIZE, //!< which clusters get reported first? (by default, the biggest ones)
                 Scalar const *const *const *aaafVoxelWeights=nullptr, //!< optional: weights of each voxel used when sort_criteria==SORT_BY_SIZE
                 #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
                 VectorContainer ***aaaafVectorStandardized=nullptr, //!< optional: place to store "standardized" vector directions
                 #endif
                 const vector<vector<array<Coordinate, 3> > > *pMustLinkConstraints=nullptr,  //!< Optional: a list of sets of voxel locations.  This insures that voxels in each set will belong to the same cluster.
                 bool start_from_saliency_maxima=true,             //!< start from local maxima? (if false, minima will be used)  WARNING: As of 2019-2-28, this function has not yet been tested with the non-default value (false)
                 ostream *pReportProgress=nullptr)  //!< print progress to the user?
{
  selfadjoint_eigen3::EigenOrderType eival_order;
  if (start_from_saliency_maxima)
    eival_order = selfadjoint_eigen3::DECREASING_EIVALS; //<--first eigenvalue will be the largest eigenvalue
  else
    eival_order = selfadjoint_eigen3::INCREASING_EIVALS; //<--first eigenvalue will be the smallest (most negative) eigenvalue

  assert(image_size);
  assert(aaafSaliency);
  assert(aaaiDest);

  if (! consider_dot_product_sign) {
    // Weird detail:
    // A negative threshold_vector should mean that the caller wants us to
    // ignore the vector argument and not attempt to compare vectors
    // associated with each voxel.
    // (At least that's how I invoke it. -AJ 2019-2-11)
    // However if (consider_dot_product_sign==false), then the absolute
    // value of the threshold is used, which will be a positive number.
    // This will cause voxels whose dot product is less than the absolute
    // value of this threshold to be ignored.
    // This is not what we want.
    // In that case we want to choose a threshold which is as low as possible
    // for the absolute value of a dot product, which is 0.0;
    if (threshold_vector_saliency < 0)
      threshold_vector_saliency = 0.0;
    if (threshold_vector_neighbor < 0)
      threshold_vector_neighbor = 0.0;
  }


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


  Scalar SIGN_FACTOR = 1.0;
  if (start_from_saliency_maxima) {
    SIGN_FACTOR = -1.0;
  }

  vector<array<Coordinate, 3> > extrema_locations; //where is the minima/maxima?
  vector<Scalar> extrema_scores;  //how bright is this minima or maxima?
  vector<size_t> extrema_nvoxels; //(needed to pacify syntax of FindExtrema3d())

  // Find all the local minima (or maxima?) in the image.

  FindExtrema(image_size,
              aaafSaliency,
              aaafMask,
              extrema_locations,
              extrema_scores,
              extrema_nvoxels,
              (! start_from_saliency_maxima), //<-- minima or maxima?
              threshold_saliency,
              connectivity,
              true,  // maxima are allowed to be located on the image border
              static_cast<Label***>(nullptr),
              pReportProgress);

  ptrdiff_t UNDEFINED = extrema_locations.size() + 1; //an impossible value
  ptrdiff_t QUEUED = extrema_locations.size() + 2; //an impossible value

  //initialize aaaiDest[][][]
  for (int iz=0; iz<image_size[2]; iz++) {
    for (int iy=0; iy<image_size[1]; iy++) {
      for (int ix=0; ix<image_size[0]; ix++) {
        if (aaafMask && aaafMask[iz][iy][ix] == 0.0)
          continue;
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
      "-- Clustering voxels belonging to different objects --\n"
      "starting from " << extrema_locations.size() << " different local "
                     << (start_from_saliency_maxima ? "maxima" : "minima")
                     << endl;

  // Initialize the queue with the voxels at these minima locations

  for (size_t i=0; i < extrema_locations.size(); i++) {
    // Create an entry in q for each of the local minima

    // Assign a different integer to each of these minima, starting at 1
    ptrdiff_t which_basin = i;

    int ix = extrema_locations[i][0];
    int iy = extrema_locations[i][1];
    int iz = extrema_locations[i][2];

    // These entries in the priority queue will be sorted by "score"
    // which is the intensity of the image at this voxel location.
    Scalar score = extrema_scores[i];
    assert(score == aaafSaliency[iz][iy][ix]);
    score *= SIGN_FACTOR; //(enable search for either local minima OR maxima)

    // Note:FindExtrema() should avoid minima above threshold_saliency,
    //     or maxima below threshold_saliency. We check for that with an assert:
    assert(score <= threshold_saliency * SIGN_FACTOR);

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

  } // for (size_t i=0; i < extrema_locations.size(); i++)


  // Each cluster contains one or more "basins"
  // A "basin" is a group of connected voxels which are all part of the
  // same local minima (or maxima).
  // Initially, each basin is it's own cluster.
  // (Later we may merge multiple basins into a single cluster.)
  vector<ptrdiff_t> basin2cluster(extrema_locations.size());
  for (size_t i=0; i < basin2cluster.size(); i++)
    basin2cluster[i] = i;

  // Inverse lookup table
  vector<set<ptrdiff_t> > cluster2basins(extrema_locations.size());
  for (size_t i=0; i < basin2cluster.size(); i++) {
    cluster2basins[i] = set<ptrdiff_t>();
    cluster2basins[i].insert(i);
  }


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
        if (SIGN_FACTOR * aaafSaliency[iz][iy][ix]
            >
            SIGN_FACTOR * threshold_saliency)
          continue;
        n_voxels_image++;
      }
    }
  }

  #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
  // initialize aaaafVectorStandardized[][][]
  if (aaaafVector && aaaafVectorStandardized) {
    for (int iz=0; iz<image_size[2]; iz++) {
      for (int iy=0; iy<image_size[1]; iy++) {
        for (int ix=0; ix<image_size[0]; ix++) {
          aaaafVectorStandardized[iz][iy][ix] = aaaafVector[iz][iy][ix];
        }
      }
    }
  }
  vector<char> basin2polarity(extrema_locations.size(), 1);
  bool voxels_discarded_due_to_polarity = false;
  Scalar voxel_discarded_due_to_polarity_saliency = -1.0;
  int voxel_discarded_due_to_polarity_ix = -1; // an impossible value
  int voxel_discarded_due_to_polarity_iy = -1; // an impossible value
  int voxel_discarded_due_to_polarity_iz = -1; // an impossible value
  #endif // #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
 
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
    Scalar i_score = -std::get<0>(p); //abs(i_score) = voxel intensity = aaafSaliency[iz][iy][ix]
    Scalar i_which_basin = std::get<1>(p); // the basin to which that voxel belongs (tentatively)
    int ix = std::get<2>(p)[0]; // voxel location
    int iy = std::get<2>(p)[1]; //   "      "
    int iz = std::get<2>(p)[2]; //   "      "

    // Should we ignore this voxel?

    if (i_score > threshold_saliency * SIGN_FACTOR) {
      // stop when the voxel brightness(*SIGN_FACTOR) exceeds the threshold_saliency
      aaaiDest[iz][iy][ix] = UNDEFINED;
      continue;
    }

    if (aaafMask && aaafMask[iz][iy][ix] == 0.0) {
      // ignore voxel if the user specified a mask and the voxel does not belong
      aaaiDest[iz][iy][ix] = UNDEFINED;
      continue;
    }


    { // Use inconsistencies in aaafSaliency to discard voxel ix,iy,iz?
      
      // -----------------------------------------------
      // First, condsider discarding voxel ix,iy,iz due to
      // inconsistencies between aaafSaliency and aaaafSymTensor here
      // -----------------------------------------------

      Scalar saliency_hessian3x3[3][3];

      CalcHessianFiniteDifferences(aaafSaliency, //!< source image
                                   ix, iy, iz,
                                   saliency_hessian3x3,
                                   image_size);

      // Confusing sign compatibility issue:
      // We assume that the tensor passed
      // by the caller will be positive-definite near a ridge.
      // In other words, it's largest eigenvalue will be positive there,
      // and its corresponding eigenvector is assumed to point in a direction
      // where the 2nd-derivative is largest.
      // In that case we must invert the sign of the hessian matrix before we
      // compare it with the tensor, because we also assume the saliency is
      // bright on a dark background.  Hence it's second derivative (ie. its 
      // Hessian) along that direction will be negative instead of positive
      // at that location.  Hence if both of these assumptions are true, 
      // we must multiply the entries in saliency_hessian by -1 before 
      // comparing them with the tensor.

      if (tensor_is_positive_definite_near_ridge ==
          start_from_saliency_maxima) {
        for (int di = 0; di < 3; di++)
          for (int dj = di; dj < 3; dj++)
            saliency_hessian3x3[di][dj] *= -1.0;
      }

      // To reduce memory consumption,
      // save the resulting 3x3 matrix in a smaller 1-D array (with 6 entries)
      // whose index is given by MapIndices_3x3_to_linear[][]

      Scalar saliency_hessian[6];

      for (int di = 0; di < 3; di++)
        for (int dj = di; dj < 3; dj++)
          saliency_hessian[ MapIndices_3x3_to_linear[di][dj] ]
            = saliency_hessian3x3[di][dj];

      bool discard_this_voxel = false;

      if (aaaafSymmetricTensor) {

        Scalar tp = TraceProductSym3(saliency_hessian, 
                                     aaaafSymmetricTensor[iz][iy][ix]);
        Scalar fs = FrobeniusNormSym3(saliency_hessian);
        Scalar ft = FrobeniusNormSym3(aaaafSymmetricTensor[iz][iy][ix]);

        if (tp < threshold_tensor_saliency * fs * ft) {
          discard_this_voxel = true;
        }
      }

      // -----------------------------------------------
      // Then, condsider discarding voxel ix,iy,iz due to
      // inconsistencies between aaafSaliency and aaaafVector at this location
      // -----------------------------------------------

      Scalar s_eivals[3];
      Scalar s_eivects[3][3];

      // the eigenvector of the saliency_hessian that we care about is the
      // one with the largest eigenvalue (which is assumed to be positive).

      ConvertFlatSym2Evects3(saliency_hessian,
                             s_eivals,
                             s_eivects,
                             eival_order);

      if (aaaafVector) {
        bool vect_threshold_exceeded = false;

        if (consider_dot_product_sign) {
          if (DotProduct3(s_eivects[0], //principal (first) eivenvector
                          aaaafVector[iz][iy][ix])
              <
              (threshold_vector_saliency *
               Length3(s_eivects[0]) *
               Length3(aaaafVector[iz][iy][ix])))
            vect_threshold_exceeded = true;
        }
        else {
          if (SQR(DotProduct3(s_eivects[0], //principal (first) eivenvector
                              aaaafVector[iz][iy][ix]))
              <
              (SQR(threshold_vector_saliency) *
               SquaredNorm3(s_eivects[0]) *
               SquaredNorm3(aaaafVector[iz][iy][ix])))
            vect_threshold_exceeded = true;
        }
        if (vect_threshold_exceeded) {
          discard_this_voxel = true;
        }
      }

      if (discard_this_voxel) {
        aaaiDest[iz][iy][ix] = UNDEFINED;
        // Are we deleting a voxel which is also a basin local minima/maxima?
        if ((ix == extrema_locations[i_which_basin][0]) &&
            (iy == extrema_locations[i_which_basin][1]) &&
            (iz == extrema_locations[i_which_basin][2]))
          // If so, then this entire basin should be deleted from consieration
          // Next line effectively deletes i_which_basin from basin2cluster[]
          basin2cluster[i_which_basin] = -1;
        continue;
      }

    } // Use inconsistencies in aaafSaliency to discard voxel ix,iy,iz?


    assert(aaaiDest[iz][iy][ix] == QUEUED);
    // Now we assign this voxel to the basin
    aaaiDest[iz][iy][ix] = i_which_basin;
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
    // Check the voxels that surround voxel (ix,iy,iz)

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


      { // Difference between voxel ix,iy,iz and ix_jx,iy_jy,iz_jz too large?

        // -----------------------------------------------
        // Then, condsider discarding voxel ix_jx,iy_jy,iz_jz due to
        // inconsistencies between aaaafVector[iz][iy][ix]
        //                     and aaaafVector[iz_jz][iy_jy][ix_jx]
        // -----------------------------------------------
        if (aaaafSymmetricTensor) {
          if (TraceProductSym3(aaaafSymmetricTensor[iz][iy][ix],
                               aaaafSymmetricTensor[iz_jz][iy_jy][ix_jx])
              <
              (threshold_tensor_neighbor *
               FrobeniusNormSym3(aaaafSymmetricTensor[iz][iy][ix]) *
               FrobeniusNormSym3(aaaafSymmetricTensor[iz_jz][iy_jy][ix_jx])))
          {
            continue;
          }
        }

        // -----------------------------------------------
        // Then, condsider discarding voxel ix_jx,iy_jy,iz_jz due to
        // inconsistencies between aaaafSymTensor[iz][iy][ix]
        //                     and aaaafSymTensor[iz_jz][iy_jy][ix_jx]
        // -----------------------------------------------
        if (aaaafSymmetricTensor) {
          if (consider_dot_product_sign) {
            if (DotProduct3(aaaafVector[iz][iy][ix],
                            aaaafVector[iz_jz][iy_jy][ix_jx])
                <
                (threshold_tensor_neighbor *
                 Length3(aaaafVector[iz][iy][ix])*
                 Length3(aaaafVector[iz_jz][iy_jy][ix_jx])))
            {
              continue;
            }
          }
          else {
            if (SQR(DotProduct3(aaaafVector[iz][iy][ix],
                                aaaafVector[iz_jz][iy_jy][ix_jx]))
                <
                (SQR(threshold_vector_neighbor) *
                 SquaredNorm3(aaaafVector[iz][iy][ix])*
                 SquaredNorm3(aaaafVector[iz_jz][iy_jy][ix_jx])))
            {
              continue;
            }
          }
        }
      } // Difference between voxel ix,iy,iz and ix_jx,iy_jy,iz_jz too large?


      if (aaaiDest[iz_jz][iy_jy][ix_jx] == QUEUED) {
        continue;
      }
      else if (aaaiDest[iz_jz][iy_jy][ix_jx] == UNDEFINED)
      {

        aaaiDest[iz_jz][iy_jy][ix_jx] = QUEUED;

        // and push this neighboring voxels onto the queue.
        // (...if they have not been assigned to a basin yet.  This
        //  insures that the same voxel is never pushed more than once.)
        array<Coordinate, 3> neighbor_crds;
        neighbor_crds[0] = ix_jx;
        neighbor_crds[1] = iy_jy;
        neighbor_crds[2] = iz_jz;
        Scalar neighbor_score = aaafSaliency[iz_jz][iy_jy][ix_jx]*SIGN_FACTOR;
        q.push(make_tuple(-neighbor_score,
                          i_which_basin,
                          neighbor_crds));


        #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
        if (aaaafVector && aaaafVectorStandardized &&
            (! consider_dot_product_sign))
        {
          if (DotProduct3(aaaafVectorStandardized[iz][iy][ix],
                          aaaafVectorStandardized[iz_jz][iy_jy][ix_jx])
              <
              0.0)
          {
            // Voxels added as the basin is growing should always have
            // directions which are compatible with each other.
            // (IE the dot product between them should not be negative.)
            // (Why?  See "Discussion" below)
            // If this is not the case, then invert the sign of the newcommers.
            // (Note: The voxels within different basins are made compatible
            //        by flipping the sign of their entry in "basin2polarity[]")
            aaaafVectorStandardized[iz_jz][iy_jy][ix_jx][0] *= -1.0;
            aaaafVectorStandardized[iz_jz][iy_jy][ix_jx][1] *= -1.0;
            aaaafVectorStandardized[iz_jz][iy_jy][ix_jx][2] *= -1.0;
            // Discussion:
            // It is possible to assume this because the topology of the
            // basins that we are growing at this point cannot have loops.
            // This is because each neighboring voxel we add in this step
            // is of type UNDEFINED.  Later we may encounter loops when we
            // run into neighboring voxels which have already beem processed.
            // We will worry about Möbius loops then.
          }
        } // if (aaaafVectorStandardized && (! consider_dot_product_sign))
        #endif


      } // else if (aaaiDest[iz_jz][iy_jy][ix_jx] == UNDEFINED)
      else {

        assert(aaaiDest[iz][iy][ix] == i_which_basin);

        ptrdiff_t basin_i = aaaiDest[iz][iy][ix];
        ptrdiff_t basin_j = aaaiDest[iz_jz][iy_jy][ix_jx];

        ptrdiff_t cluster_i = basin2cluster[basin_i];
        ptrdiff_t cluster_j = basin2cluster[basin_j];

        #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
        bool polarity_match = true;
        if (aaaafVectorStandardized && (! consider_dot_product_sign))
        {
          if (DotProduct3(aaaafVectorStandardized[iz][iy][ix],
                          aaaafVectorStandardized[iz_jz][iy_jy][ix_jx])
              *
              basin2polarity[basin_i]
              *
              basin2polarity[basin_j]
              <
              0.0)
            polarity_match = false;
        }
        #endif

        if (cluster_i == cluster_j) {
          #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
          if (! polarity_match) {
            voxels_discarded_due_to_polarity = true;
            voxel_discarded_due_to_polarity_ix = ix_jx;
            voxel_discarded_due_to_polarity_iy = iy_jy;
            voxel_discarded_due_to_polarity_iz = iz_jz;
            voxel_discarded_due_to_polarity_saliency = aaafSaliency[iz_jz][iy_jy][ix_jx];

            // In that case throw away this voxel.
            // Hopefully this will prevent linking of voxels containing
            // vector directors which cannot be reconciled.  Example:
            // This should prevent surfaces which resemble a Möbius strip
            // from forming a complete closed loop.  Hopefully it will
            // cut the loop at the voxels with the most tenuous connections.
            continue;
          }
          #endif
        }
        else //if (cluster_i != cluster_j)
        {
          // -- merge the two clusters to which basin_i and basin_j belong --

          // (arbitrarily) choose the cluster with the smaller ID number
          // to absorb the other cluster, and delete the other cluster.
          ptrdiff_t merged_cluster_id  = std::min(cluster_i, cluster_j);
          ptrdiff_t deleted_cluster_id = std::max(cluster_i, cluster_j);
          assert(cluster_i != cluster_j);

          // copy the basins from the deleted cluster into the merged cluster
          for (auto p = cluster2basins[deleted_cluster_id].begin();
               p != cluster2basins[deleted_cluster_id].end();
               p++)
          {
            ptrdiff_t basin_id = *p;
            cluster2basins[merged_cluster_id].insert(basin_id);
            basin2cluster[basin_id] = merged_cluster_id;

            #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
            if (aaaafVector && aaaafVectorStandardized &&
                (! consider_dot_product_sign))
            {
              if (! polarity_match)
                basin2polarity[basin_id] *= -1.0;
            }
            #endif
          }

          // -- delete all the basins from the deleted cluster --
          cluster2basins[deleted_cluster_id].clear();

        } // if (aaaiDest[iz_jz][iy_jy][ix_jx] != aaaiDest[iz][iy][ix])
      } // else clause for "if (aaaiDest[iz_jz][iy_jy][ix_jx] == UNDEFINED)"
    } // for (int j = 0; j < num_neighbors; j++)
  } //while (q.size() != 0)


  #ifndef NDEBUG
  // DEBUGGING
  // All of the voxels should either be assigned to a basin or to "UNDEFINED"
  // (but not "QUEUED").  Check for that below:
  for (int iz=0; iz<image_size[2]; iz++)
    for (int iy=0; iy<image_size[1]; iy++)
      for (int ix=0; ix<image_size[0]; ix++)
        assert(aaaiDest[iz][iy][ix] != QUEUED);
  #endif //#ifndef NDEBUG



  if (pMustLinkConstraints)
  {

    // loop over the different groups of voxels
    //    (we will force each voxel in a group to belong to the same cluster)

    for (size_t i_group = 0; i_group < pMustLinkConstraints->size(); i_group++)
    {

      // some variables we will need

      ptrdiff_t FIRST_ITER = -1; // an impossible value used below

      //r_i_init = coordinates of most recent voxel selected by the user
      array<int,3> r_i_init={-1,-1,-1};
      //r_i = coordinates of the voxel in aaaiDest[][][] nearest to r_i_init
      array<int,3> r_i;
      //basin_i = the ID number for the basin to which voxel r_i belongs
      ptrdiff_t basin_i = FIRST_ITER;

      //r_j_init, r_j, basin_j correspond to the PREVIOUSLY processed voxel
      ptrdiff_t basin_j = FIRST_ITER;
      array<int,3> r_j_init={-1,-1,-1}; //an impossible value
      array<int,3> r_j={-1,-1,-1};      //an impossible value


      // Loop over the voxels in each group 
      // and force them to belong to the same cluster

      for (auto pLocation = (*pMustLinkConstraints)[i_group].begin();
           pLocation != (*pMustLinkConstraints)[i_group].end();
           pLocation++)
      {
        r_i_init[0] = int(floor((*pLocation)[0]+0.5));
        r_i_init[1] = int(floor((*pLocation)[1]+0.5));
        r_i_init[2] = int(floor((*pLocation)[2]+0.5));
        if (*pReportProgress) {
          *pReportProgress << "  finding voxel nearest to (";
          for (int d = 0; d < 3; d++) {
            *pReportProgress << r_i_init[d];
            if (d+1 != 3) *pReportProgress << ",";
          }
          *pReportProgress << ") -> (";
        }

        set<Label> ignore_these_voxel_types;
        ignore_these_voxel_types.insert(UNDEFINED);

        FindNearestVoxel(image_size,
                         r_i_init, //find voxel in aaaiDest closest to this
                         r_i, //store the location here
                         aaaiDest, //inverse lookup voxel location->basinID
                         aaafMask,
                         ignore_these_voxel_types, // skip over these voxels
                         true); //invert selection (skip over them)

        if (*pReportProgress) {
          for (int d = 0; d < 3; d++) {
            *pReportProgress << r_i[d];
            if (d+1 != 3) *pReportProgress << ",";
          }
          *pReportProgress << ")";
          if (aaaafVectorStandardized) {
            *pReportProgress << ", normal=(";
            for (int d = 0; d < 3; d++) {
              aaaafVectorStandardized[r_i[2]][r_i[1]][r_i[0]][d];
              *pReportProgress << aaaafVectorStandardized[r_i[2]][r_i[1]][r_i[0]][d];
              if (d+1 != 3) *pReportProgress << ",";
            }
            *pReportProgress << ")";
          }
          *pReportProgress << "\n";
        }

        if ((r_i[0] == -1) && (r_i[1] == -1) && (r_i[2] == -1))
          throw VisfdErr("Error: No voxels clustered. Empty image. Your cluster criteria are too strict.\n"
                         "       (Attempting the find the nearest voxel from an empty set.)\n");

        basin_i = aaaiDest[r_i[2]][r_i[1]][r_i[0]];

        assert((basin_i != UNDEFINED) && (basin_i != QUEUED));

        if ((basin_j != FIRST_ITER) && (basin_i != basin_j))
        {
          // Then basin_i should be linked with basin_j 
          // by merging the clusters to which they belong.

          ptrdiff_t cluster_i = basin2cluster[basin_i];
          ptrdiff_t cluster_j = basin2cluster[basin_j];

          if (cluster_i != cluster_j)
          {
            ptrdiff_t merged_cluster_id  = std::min(cluster_i, cluster_j);
            ptrdiff_t deleted_cluster_id = std::max(cluster_i, cluster_j);

            #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
            // Try to infer whether or not the voxels located at r_i and r_j
            // are pointing in a compatible direction.
            // If not, set "polarity_match" to false.
            // Later we will flip polarity of the voxels in one of the 
            // joined clusters in order to make sure the resulting cluster.
            // Note:
            // Since r_i and r_j could be separated by a great distance
            // the surface direction (or curve direction) is likely to be
            // not pointing in the same direction, or even within 90 degrees
            // of that direction. It makes no sense to check if they are
            // pointing in the same direction.  Instead, we try to infer
            // if they are pointing in a "compatible" direction.
            bool polarity_match;
            array<Scalar, 3> r_ij, n_i, n_j, ni_cross_nj;
            for (int d = 0; d < 3; d++) {
              n_i[d] = aaaafVectorStandardized[r_i[2]][r_i[1]][r_i[0]][d];
              n_j[d] = aaaafVectorStandardized[r_j[2]][r_j[1]][r_j[0]][d];
              r_ij[d] = r_i[d] - r_j[d];
            }
            Normalize3(r_ij);
            Scalar ni_dot_rij = DotProduct3(n_i, r_ij);
            Scalar nj_dot_rij = DotProduct3(n_j, r_ij);
            // If each plane is at an angle less than 45 degrees from r_ij
            // then consider the two planes to be facing the same direction
            // and part of the same portion of the curve.
            Scalar theta0 = M_PI/4;  // (45 degrees)
            Scalar theta_ni_rij = asin(abs(ni_dot_rij));
            Scalar theta_nj_rij = asin(abs(nj_dot_rij));
            if (*pReportProgress) {
              *pReportProgress << "    (theta_ni_rij == "
                               << theta_ni_rij*180.0/M_PI
                               << ", theta_nj_rij == "
                               << theta_nj_rij*180.0/M_PI << ")\n";
            }
            if ((theta_ni_rij < theta0) &&
                (theta_nj_rij < theta0))
            {
              polarity_match = (DotProduct3(n_i, n_j) > 0);
            }
            else {
              // otherwise, the two planes are probably facing the opposite
              // direction.  In that case, use the following criteria:
              polarity_match = (ni_dot_rij * nj_dot_rij <= 0);
            }

            // Confusing code:
            // We only flip the polarity of the basins from the deleted cluster
            // if the basin2polarity[basin_i] & basin2polarity[basin_j] entries
            // are NOT consistent with "polarity_match", which we observed
            // directly (for example by calculating the dot products of the
            // surface normals at these two locations r_i and r_j).
            
            bool flip_polarity = (polarity_match !=
                                  (basin2polarity[basin_i] ==
                                   basin2polarity[basin_j]));

            #endif //#ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION

            // copy the basins from the deleted cluster into the merged cluster
            for (auto p = cluster2basins[deleted_cluster_id].begin();
                 p != cluster2basins[deleted_cluster_id].end();
                 p++)
            {
              ptrdiff_t basin_id = *p;
              cluster2basins[merged_cluster_id].insert(basin_id);
              basin2cluster[basin_id] = merged_cluster_id;

              #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
              if (aaaafVector && aaaafVectorStandardized &&
                  (! consider_dot_product_sign))
              {
                if (flip_polarity)
                  basin2polarity[basin_id] *= -1.0;
              }
              #endif
            }

            // -- delete all the basins from the deleted cluster --
            cluster2basins[deleted_cluster_id].clear();

          } //if (cluster_i != cluster_j) {


        } //if ((basin_j != FIRST_ITER) && (basin_i != basin_j))

        basin_j = basin_i;
        r_j_init = r_i_init;
        r_j = r_i;

      } //for (auto pLocation = (*pMustLinkConstraints)[i].begin(); ...)
      
    } //for (size_t i = 0; i < pMustLinkConstraints->size(); i++)
    
  } // if (pMustLinkConstraints)


  // Optional: We don't need "cluster2basins" any more.  Free up its memory.
  cluster2basins.clear();


  // Count the number of clusters
  size_t n_clusters = 0;
  size_t n_basins = basin2cluster.size();
  vector<size_t> clusterold2clusternew(n_basins);
  vector<size_t> cluster2deepestbasin; //inverse lookup table
  for (size_t i=0; i < n_basins; i++) {
    // (Originally, clusters were assigned to the ID number of the basin they
    //  started from.  But after merging multiple basins into the same cluster
    //  there are more basins than clusters.)
    clusterold2clusternew[i] = n_clusters;
    if (basin2cluster[i] == i) {
      cluster2deepestbasin.push_back(i);
      assert(cluster2deepestbasin[n_clusters] == i);
      n_clusters++;
    }
  }

  if (pReportProgress)
    *pReportProgress <<
      "Number of clusters found: " << n_clusters << "\n";

  // Renumber the clusters from 0 ... n_clusters-1
  for (size_t i=0; i < basin2cluster.size(); i++) {
    basin2cluster[i] = clusterold2clusternew[basin2cluster[i]];
  }



  #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
  if (aaaafVector && aaaafVectorStandardized &&
    (! consider_dot_product_sign))
  {
    // Now, finally incorporate the data we stored in basin2polarity[] into
    // aaaafVectorStandardized.
    // (We didn't do this earlier because entire basins flip polarity often)
    for (int iz=0; iz<image_size[2]; iz++) {
      for (int iy=0; iy<image_size[1]; iy++) {
        for (int ix=0; ix<image_size[0]; ix++) {
          if ((aaafMask && aaafMask[iz][iy][ix] == 0.0) ||
              (aaaiDest[iz][iy][ix] == UNDEFINED)) {
            aaaafVectorStandardized[iz][iy][ix][0] = 0.0;
            aaaafVectorStandardized[iz][iy][ix][1] = 0.0;
            aaaafVectorStandardized[iz][iy][ix][2] = 0.0;
          }
          else {
            ptrdiff_t basin_id = aaaiDest[iz][iy][ix];
            assert((0 <= basin_id) && (basin_id < n_basins));
            Scalar polarity = basin2polarity[basin_id];
            aaaafVectorStandardized[iz][iy][ix][0] *= polarity;
            aaaafVectorStandardized[iz][iy][ix][1] *= polarity;
            aaaafVectorStandardized[iz][iy][ix][2] *= polarity;
          }
        }
      }
    }

    if (pReportProgress && voxels_discarded_due_to_polarity) {
      *pReportProgress
        << "WARNING: During clustering, some voxels were discarded to avoid\n"
        << "         directional singularities.  More specifically:\n"
        << "         some voxels were discarded in order to \"cut\" the conncted clusters\n"
        << "         so that each remaining object (cluster) consist of voxels which\n"
        << "         are consistently orientable.\n"
        << "         (For example, a Möbius strip forms a closed loop whose surface\n"
        << "          normals are not consistently orientable.  Such a loop would be cut\n"
        << "          somewhere along its length.)\n"
        << "         This cutting was done where the saliency (brightness/darkness) of the\n"
        << "         object was weakest. The first such voxel discarded (cut) had saliency\n"
        << "         " << voxel_discarded_due_to_polarity_saliency << "\n"
        << "         ...and was located at position: ("
        << voxel_discarded_due_to_polarity_ix << ", "
        << voxel_discarded_due_to_polarity_iy << ", "
        << voxel_discarded_due_to_polarity_iz << ")\n"
        << "         (Note: Zero-indexing.  The first voxel has position 0,0,0 not 1,1,1)\n"
        << endl;
    } // if (pReportProgress && voxels_discarded_due_to_polarity)
  } // if (aaaafVectorStandardized && (! consider_dot_product_sign))
  #endif // #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION


  // Now, assign all the voxels in aaaiDest[][][]
  // to their clusters instead of their basins.
  for (int iz=0; iz<image_size[2]; iz++) {
    for (int iy=0; iy<image_size[1]; iy++) {
      for (int ix=0; ix<image_size[0]; ix++) {
        if (aaaiDest[iz][iy][ix] == UNDEFINED)
          continue;
        if (aaafMask && aaafMask[iz][iy][ix] == 0.0)
          continue;
        ptrdiff_t basin_id = aaaiDest[iz][iy][ix];
        assert((0 <= basin_id) && (basin_id < n_basins));
        aaaiDest[iz][iy][ix] = basin2cluster[ basin_id ];
      }
    }
  }

  // Keep track of the number of voxels in each cluster.
  // I call this the "size" of each cluster.
  //
  // (If the caller supplied a non-null "aaafVoxelWeights" array, then instead
  //  I define the "size" of a cluster as the sum of the entries in the
  //  "aaafVoxelWeights[][][]" array for the voxels within that cluster.)
  //
  // Later we will sort the clusters by their sizes.

  vector<long double> cluster_sizes(n_clusters, 0.0);

  // I use floating point numbers because "cluster_sizes" need not be integers.
  // (Note: "vector<Label>" or "vector<Scalar>" won't work here since, in
  //        typical usage, both Label and Scalar are of type "float".
  //        Typical MRC/REC files have 10^9 voxels, and adding 10^9 32-bit
  //        floats together will definitely result in numeric underflow.
  //        So we use "vector<double>" or "vector<long double>" precision 
  //        instead to make sure we get an accurate total when adding.

  for (int iz=0; iz<image_size[2]; iz++) {
    for (int iy=0; iy<image_size[1]; iy++) {
      for (int ix=0; ix<image_size[0]; ix++) {
        if (aaaiDest[iz][iy][ix] == UNDEFINED)
          continue;
        ptrdiff_t cluster_id = aaaiDest[iz][iy][ix];
        if (aaafVoxelWeights)
          cluster_sizes[cluster_id] += aaafVoxelWeights[iz][iy][ix];
        else
          cluster_sizes[cluster_id] += 1.0;
      }
    }
  }


  #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION
  if (aaaafVectorStandardized && (! consider_dot_product_sign))
  {
    // At this point we have got the voxels pointing in directions which are
    // consistent with each other at least, but we have not determined
    // whether they are currently pointing outside or inside.
    // Assuming that this cluster represents a closed surface,
    // we would like them to be pointing outward.
    //
    // We use an extremely crude method to decide this:
    // If the majority of voxels in a cluster are pointing away from 
    // the cluster's center of mass, then we say that they are pointing outward.
    // (More precisely, if the sum of the dot-products is positive...)
    // If they are not, then we should flip all of the voxel directions...

    // First, calculate the center of mass of each voxel
    vector<array<long double, 3> > cluster_centers_of_mass(n_clusters);
    for (int iz=0; iz<image_size[2]; iz++) {
      for (int iy=0; iy<image_size[1]; iy++) {
        for (int ix=0; ix<image_size[0]; ix++) {
          if (aaafMask && aaafMask[iz][iy][ix] == 0.0)
            continue;
          if (aaaiDest[iz][iy][ix] == UNDEFINED)
            continue;
          ptrdiff_t cluster_id = aaaiDest[iz][iy][ix];
          if (aaafVoxelWeights) {
            // Example:
            // Suppose the voxel directions correspond to surface normals.
            // One reason the caller might pass a non-null aaafVoxelWeights
            // array would be to store the surface area corresponding to
            // each voxel.  If that's the case, then we are performing
            // a surface area weighted center-of-mass for each cluster.
            cluster_centers_of_mass[cluster_id][0] +=
              ix * aaafVoxelWeights[iz][iy][ix];
            cluster_centers_of_mass[cluster_id][1] +=
              iy * aaafVoxelWeights[iz][iy][ix];
            cluster_centers_of_mass[cluster_id][2] +=
              iz * aaafVoxelWeights[iz][iy][ix];
          }
          else {
            // Otherwise just calculate the average position of each cluster
            cluster_centers_of_mass[cluster_id][0] += ix;
            cluster_centers_of_mass[cluster_id][1] += iy;
            cluster_centers_of_mass[cluster_id][2] += iz;
          }
        }
      }
    }
    for (int i = 0; i < n_clusters; i++)
      for (int d = 0; d < 3; d++)
        // cluster_size[i] = sum of the weights of the i'th cluster
        cluster_centers_of_mass[i][d] /= cluster_sizes[i];

    // Do the majority of voxels point toward or away from the center of mass?
    vector<long double> sum_dot_products(n_clusters, 0.0);
    for (int iz=0; iz<image_size[2]; iz++) {
      for (int iy=0; iy<image_size[1]; iy++) {
        for (int ix=0; ix<image_size[0]; ix++) {
          if (aaafMask && aaafMask[iz][iy][ix] == 0.0)
            continue;
          if (aaaiDest[iz][iy][ix] == UNDEFINED)
            continue;
          ptrdiff_t cluster_id = aaaiDest[iz][iy][ix];
          Scalar r_minus_rcom[3];
          r_minus_rcom[0] = ix - cluster_centers_of_mass[cluster_id][0];
          r_minus_rcom[1] = iy - cluster_centers_of_mass[cluster_id][1];
          r_minus_rcom[2] = iz - cluster_centers_of_mass[cluster_id][2];
          long double delta_sum =
            DotProduct3(r_minus_rcom, aaaafVectorStandardized[iz][iy][ix]);
          if (aaafVoxelWeights) {
            // Example:
            // Suppose the voxel directions correspond to surface normals.
            // One reason the caller might pass a non-null aaafVoxelWeights
            // array would be to store the surface area corresponding to
            // each voxel.  If that's the case, then the
            // sum we are calculating here is the surface integral of the
            // dot product between "r_minus_rcom" and the voxel direction.
            delta_sum *= aaafVoxelWeights[iz][iy][ix];
          }
          sum_dot_products[cluster_id] += delta_sum;
        } // for (int ix=0; ix<image_size[0]; ix++)
      } // for (int iy=0; iy<image_size[1]; iy++)
    } // for (int iz=0; iz<image_size[2]; iz++)

    // note:the next 3 lines are not necessary since we only care about the sign
    //for (int i = 0; i < n_clusters; i++)
    //  // cluster_size[i] = sum of the weights of the i'th cluster
    //  sum_dot_products[i] /= cluster_size[i];

    for (int iz=0; iz<image_size[2]; iz++) {
      for (int iy=0; iy<image_size[1]; iy++) {
        for (int ix=0; ix<image_size[0]; ix++) {
          if (aaafMask && aaafMask[iz][iy][ix] == 0.0)
            continue;
          if (aaaiDest[iz][iy][ix] == UNDEFINED)
            continue;
          ptrdiff_t cluster_id = aaaiDest[iz][iy][ix];
          if (sum_dot_products[cluster_id] < 0.0)
            for (int d = 0; d < 3; d++)
              aaaafVectorStandardized[iz][iy][ix][d] *= -1.0;
        }
      }
    }
  } // if (aaaafVectorStandardized && (! consider_dot_product_sign))
  #endif // #ifndef DISABLE_STANDARDIZE_VECTOR_DIRECTION



  if (pv_cluster_maxima != nullptr) {
    pv_cluster_maxima->resize(n_clusters);
    for (size_t i = 0; i < n_clusters; i++)
      (*pv_cluster_maxima)[i] = extrema_locations[ cluster2deepestbasin[i] ];
  }

  if (pv_cluster_saliencies) {
    pv_cluster_saliencies->resize(n_clusters);
    for (size_t i = 0; i < n_clusters; i++) {
      int ix = extrema_locations[ cluster2deepestbasin[i] ][0];
      int iy = extrema_locations[ cluster2deepestbasin[i] ][1];
      int iz = extrema_locations[ cluster2deepestbasin[i] ][2];
      (*pv_cluster_saliencies)[i] = aaafSaliency[iz][iy][ix];
    }
  }
                   
  if (pv_cluster_sizes) {
    pv_cluster_sizes->resize(n_clusters);
    for (size_t i = 0; i < n_clusters; i++)
      (*pv_cluster_sizes)[i] = cluster_sizes[i];
  }
    
  // In what order should we present the voxels?
  if ((sort_criteria == SORT_BY_SALIENCY) && (n_clusters > 0)) {
    // Do nothing.
    // (The blobs are already sorted by the height of their saliency maxima.)
  }
  else if ((sort_criteria == SORT_BY_SIZE) && (n_clusters > 0))
  {
    vector<tuple<Scalar, size_t> > size_index(n_clusters);

    for (size_t i = 0; i < n_clusters; i++)
      size_index[i] = make_tuple(cluster_sizes[i], i);

    if (pReportProgress)
      *pReportProgress << "  sorting these "<<n_clusters<<" clusters according to their sizes... ";

    // Sort the clusters in descending order,
    // with the largest clusters appearing first.
    sort(size_index.rbegin(),
         size_index.rend());

    vector<size_t> permutation(n_clusters);
    for (size_t i = 0; i < size_index.size(); i++)
      permutation[i] = get<1>(size_index[i]);
    size_index.clear();

    vector<size_t> permutation_inv;
    invert_permutation(permutation, permutation_inv);

    if (pReportProgress)
      *pReportProgress << "done\n";

    if (pv_cluster_maxima)
      apply_permutation(permutation, *pv_cluster_maxima);

    for (int iz=0; iz<image_size[2]; iz++) {
      for (int iy=0; iy<image_size[1]; iy++) {
        for (int ix=0; ix<image_size[0]; ix++) {
          if (aaaiDest[iz][iy][ix] == UNDEFINED)
            continue;
          ptrdiff_t cluster_id = aaaiDest[iz][iy][ix];
          aaaiDest[iz][iy][ix] = permutation_inv[cluster_id];
        }
      }
    }

  } // if (aaafVoxelWeights != nullptr)



  // Now, deal with voxels which are "undefined"
  // voxels to have a high instead of a low value.
  if (undefined_label_is_max)
    label_undefined = n_clusters+1;

  // Final procesing: Replace "UNDEFINED" with "label_undefined", or add 1:
  for (int iz=0; iz<image_size[2]; iz++) {
    for (int iy=0; iy<image_size[1]; iy++) {
      for (int ix=0; ix<image_size[0]; ix++) {

        if (aaafMask && aaafMask[iz][iy][ix] == 0.0)
          // as before: ignore voxels excluded by the mask
          continue;

        if (aaaiDest[iz][iy][ix] == UNDEFINED) {
          aaaiDest[iz][iy][ix] = label_undefined;
          continue;
        }

        // According to the function specifications, the aaaiDest[][][] array
        // should store the cluster to which the voxel belongs.  This should be
        // a number between 1 and n_clusters.  Currently, the number in the
        // aaaiDest[][][] array lies between 0 and n_clusters-1 (0-indexing).
        assert((0 <= aaaiDest[iz][iy][ix]) && (aaaiDest[iz][iy][ix]<n_basins));

        // Hence we should add 1 to it:
        aaaiDest[iz][iy][ix] += 1;
      }
    }
  }

  #ifndef NDEBUG
  {
    vector<bool> clusters_visited(n_clusters, false);

    for (int iz=0; iz<image_size[2]; iz++) {
      for (int iy=0; iy<image_size[1]; iy++) {
        for (int ix=0; ix<image_size[0]; ix++) {
          if (aaafMask && aaafMask[iz][iy][ix] == 0.0)
            // as before: ignore voxels excluded by the mask
            continue;
          if (aaaiDest[iz][iy][ix] == UNDEFINED)
            continue;
          assert((0<aaaiDest[iz][iy][ix]) && (aaaiDest[iz][iy][ix]<=n_basins));
          long cluster_id = aaaiDest[iz][iy][ix] - 1;
          clusters_visited[cluster_id] = true;
        }
      }
    }
    // Check to make sure that every cluster in the list
    // has at least one voxel associated with it.
    for (int i=0; i < n_clusters; i++)
      assert(clusters_visited[i] == true);
    // NOTE: THIS ASSERTION COULD FAIL even if the program is working when 
    //       the image is large and contains millions of clusters.
    //       This is because aaaiDest[iz][iy][ix] is of type "Label" which
    //       is usually "float", and single-precision floats cannot represent
    //       integers above 1e+06 (I forget the exact number.)
  }
  #endif // #ifndef NDEBUG

  delete [] neighbors;

} // ClusterConnected()


} //namespace visfd



#endif //#ifndef _CLUSTERING_HPP
