#ifndef _UNSUPPORTED_HPP
#define _UNSUPPORTED_HPP

/// UNSUPPORTED CODE
/// DEPRECIATION WARNING:
/// THIS FILE IMPLEMENTS FEATURES WHICH ARE NO LONGER MAINTAINED AND TESTED.
/// THIS CODE WILL LIKELY BE DELETED IN THE FUTURE.


#include <filter3d.hpp>
#include <mrc_simple.hpp>
#include <random_gen.h>
#include "settings.hpp"



/// @brief  This function is used to cluster voxels which are connected 
///         together and belong to the same object (for example, a surface).
///         Based on the watershed algorithm, it starts with voxels
///         which lie near a minima (or maxima) of saliency, 
///         and groups voxels of similar saliency,
///         (and compatible tensor direction) together.
/// @note   Neighboring voxels whose direction (vector) changes discontinuously,
///         will not be added to an existing cluster.
/// @note   Neighboring voxels whose direction (tensor) changes discontinuously,
///         will not be added to an existing cluster.
/// @note   In addition,
///         voxels whose saliency Hessian and tensor are incompatible
///         will not be added to any cluster.
/// @note   Voxels whose saliency Hessian's principal eigenvector
///         and vector are incompatible will not be added to any cluster.
///
/// @note DELETE THE NEXT COMMENT (not clear if relevant).
///
/// @code
/// Suppose tensor A is approximated by its principal eigenvector A~=l1|e1><e1|
/// Suppose tensor B is approximated by its principal eigenvector B~=L1|E1><E1|
/// Then  trace(A B) = l1 * L1 * dot_product(e1, E1)^2  =  l1 * L1 * <e1|E1>^2
/// @endcode
/// Proof:
/// @code
///   <e1|E1>^2 =       <E1|e1>  <e1|E1>
///             = <E1| (sum_i|i><i|) |e1><e1|E1>  (resolution of the identity)
///             = sum_i <i| |e1><e1| |E1><E1| |i>
///             = trace( |e1><e1||E1><E1| )
/// @endcode
/// Hence the trace-product (trace(A B)) might be a cheap way to decide whether
/// two tensors' (A,B) eigenvectors are pointing a compatible directions
/// (without needing to diagonalize them to find their eigenvectors).

template<class Scalar, class Label, class Coordinate, class VectorContainer, class TensorContainer>
void
ConnectedClusters(int const image_size[3],                   //!< #voxels in xyz
                  Scalar const *const *const *aaafSaliency,  //!< intensity of each voxel
                  Label ***aaaiDest,                       //!< watershed segmentation results go here
                  Scalar const *const *const *aaafMask,    //!< optional: Ignore voxels whose mask value is 0
                  Scalar threshold_saliency=std::numeric_limits<Scalar>::infinity(), //!< don't segment voxels exceeding this height
                  Label label_undefined = 0,               //!< voxels storing this value do not belong to any clusters
                  bool undefined_label_is_max = false,     //!< set label_undefined to number of clusters + 1 (overrides label_undefined)
                  VectorContainer const *const *const *aaaafVector=NULL,
                  Scalar threshold_vector_saliency=-1.0,    //!< voxels with incompatible saliency and vector are ignored (=1.0 disables)
                  Scalar threshold_vector_neighbor=-1.0,    //!< neighboring voxels with incompatible vectors are ignored (=1.0 disables)
                  TensorContainer const *const *const *aaaafSymmetricTensor=NULL,
                  Scalar threshold_tensor_saliency=-1.0,    //!< voxels with incompatible saliency and tensor are ignored (=1.0 disables)
                  Scalar threshold_tensor_neighbor=-1.0,    //!< neighboring voxels with incompatible tensors are ignored (=1.0 disables)
                  bool start_from_minima=true,             //!< start from local minima? (if false, maxima will be used)
                  int connectivity=3,                      //!< square root of the search radius around each voxel (1=nearest_neighbors, 2=2D_diagonal, 3=3D_diagonal)
                  EigenOrderType eival_order = selfadjoint_eigen3::INCREASING_EIVALS,
                  vector<array<Coordinate, 3> > *pv_cluster_centers=NULL, //!< optional: the location of saliency minima or maxima which seeded each cluster
                  ostream *pReportProgress=NULL)  //!< print progress to the user?
{

  assert(image_size);
  assert(aaafSaliency);
  assert(aaaiDest);

  // Figure out which neighbors to consider when searching neighboring voxels
  int (*neighbors)[3] = NULL; //a pointer to a fixed-length array of 3 ints
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
  if (! start_from_minima) {
    SIGN_FACTOR = -1.0;
    if (threshold_saliency == std::numeric_limits<Scalar>::infinity())
      // Then the caller did not bother to specify a threshold.  Note: the
      // default threshold assumes we are searching for minima, not maxima.
      // If the caller did not specify a halting threshold, AND if we ARE
      // looking for maxima (not minima), then we have to reverse the sign of
      // the default threshold as well (this is a little confusing to me too):
      threshold_saliency = -std::numeric_limits<Scalar>::infinity();
  }

  vector<array<Coordinate, 3> > extrema_locations; //where is the minima/maxima?
  vector<Scalar> extrema_scores;  //how bright is this minima or maxima?
  vector<size_t> extrema_nvoxels; //(needed to pacify syntax of FindExtrema3d())

  // Find all the local minima (or maxima?) in the image.

  FindExtrema3D(image_size,
                aaafSaliency,
                aaafMask,
                extrema_locations,
                extrema_scores,
                extrema_nvoxels,
                start_from_minima, //<-- minima or maxima?
                threshold_saliency,
                connectivity,
                true,
                static_cast<Scalar***>(NULL),
                pReportProgress);

  Label UNDEFINED = extrema_locations.size() + 1; //an impossible value
  Label QUEUED = extrema_locations.size() + 2; //an impossible value

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
                   Label,     // the basin-ID to which the voxel belongs
                   array<Coordinate, 3>  // location of the voxel
                  >
                > q;
               

  if (pReportProgress)
    *pReportProgress <<
      " ---- Clustering voxels belonging to different objects ----\n"
      "starting from " << extrema_locations.size() << " different local "
                       << (start_from_minima ? "minima" : "maxima") << endl;

  // Initialize the queue with the voxels at these minima locations

  for (size_t i=0; i < extrema_locations.size(); i++) {
    // Create an entry in q for each of the local minima

    // Assign a different integer to each of these minima, starting at 1
    Label which_basin = i;

    int ix = extrema_locations[i][0];
    int iy = extrema_locations[i][1];
    int iz = extrema_locations[i][2];

    // These entries in the priority queue will be sorted by "score"
    // which is the intensity of the image at this voxel location.
    Scalar score = extrema_scores[i];
    assert(score == aaafSaliency[iz][iy][ix]);
    score *= SIGN_FACTOR; //(enable search for either local minima OR maxima)

    // Note:FindExtrema3D() should avoid minima above threshold_saliency,
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

  // Initially each basin is it's own cluster.
  // (Later we may merge multiple basins into a single cluster.)
  vector<Label> basin2cluster(extrema_locations.size());
  for (size_t i=0; i < basin2cluster.size(); i++)
    basin2cluster[i] = i;

  // Loop over all the voxels on the periphery
  // of the voxels with lower intensity (brightness).
  // This is a priority-queue of voxels which lie adjacent to 
  // the voxels which have been already considered.
  // As we consider new voxels, add their neighbors to this queue as well
  // until we process all voxels in the image.

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
        // count the number of voxels in the image we will be looping over later
        n_voxels_image++;
      }
    }
  }

  //#ifndef NDEBUG
  //set<array<Coordinate, 3> > visited;
  //#endif

  while (! q.empty())
  {
    tuple<Scalar, Label, array<Coordinate, 3> > p = q.top();
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
      CalcHessianFiniteDifferences3D(aaafSaliency, //!< source image
                                     ix, iy, iz,
                                     saliency_hessian3x3);
      Scalar saliency_hessian[6];

      // To reduce memory consumption,
      // save the resulting 3x3 matrix in a smaller 1-D array (with 6 entries)
      // whose index is given by MapIndices_3x3_to_linear[][]
      for (int di = 0; di < 3; di++)
        for (int dj = di; dj < 3; dj++)
          saliency_hessian[ MapIndices_3x3_to_linear[di][dj] ]
            = saliency_hessian3x3[di][dj];

      if (aaaafSymmetricTensor) {

        float tp = TraceProductSym3(saliency_hessian, 
                                    aaaafSymmetricTensor[iz][iy][ix]);
        float fs = FrobeniusNormSym3(saliency_hessian);
        float ft = FrobeniusNormSym3(aaaafSymmetricTensor[iz][iy][ix]);

        if (TraceProductSym3(saliency_hessian, 
                             aaaafSymmetricTensor[iz][iy][ix])
            <
            (threshold_tensor_saliency *
             FrobeniusNormSym3(saliency_hessian) *
             FrobeniusNormSym3(aaaafSymmetricTensor[iz][iy][ix])))
        {
          aaaiDest[iz][iy][ix] = UNDEFINED;
          continue;
        }
      }

      // -----------------------------------------------
      // Then, condsider discarding voxel ix,iy,iz due to
      // inconsistencies between aaafSaliency and aaaafVector at this location
      // -----------------------------------------------

      Scalar s_eivals[3];
      Scalar s_eivects[3][3];
      ConvertFlatSym2Evects3(saliency_hessian,
                             s_eivals,
                             s_eivects,
                             eival_order);

      if (aaaafVector) {
        if (SQR(DotProduct3(s_eivects[0], //principal eivenvector
                            aaaafVector[iz][iy][ix]))
            <
            (SQR(threshold_vector_saliency) *
             SquaredNorm3(s_eivects[0]) *
             SquaredNorm3(aaaafVector[iz][iy][ix])))
        {
          aaaiDest[iz][iy][ix] = UNDEFINED;
          continue;
        }
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
            continue;
        }

        // -----------------------------------------------
        // Then, condsider discarding voxel ix_jx,iy_jy,iz_jz due to
        // inconsistencies between aaaafSymTensor[iz][iy][ix]
        //                     and aaaafSymTensor[iz_jz][iy_jy][ix_jx]
        // -----------------------------------------------
        if (aaaafSymmetricTensor) {
          if (SQR(DotProduct3(aaaafVector[iz][iy][ix],
                              aaaafVector[iz_jz][iy_jy][ix_jx]))
              <
              (SQR(threshold_tensor_neighbor) *
               SquaredNorm3(aaaafVector[iz][iy][ix])*
               SquaredNorm3(aaaafVector[iz_jz][iy_jy][ix_jx])))
            continue;
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
      }
      else {
        if (aaaiDest[iz_jz][iy_jy][ix_jx] != aaaiDest[iz][iy][ix])
        {
          assert(aaaiDest[iz][iy][ix] == i_which_basin);

          // merge the two basins into the same cluster
          Label basin_i = aaaiDest[iz][iy][ix];
          Label basin_j = aaaiDest[iz_jz][iy_jy][ix_jx];
          Label merged_cluster_id = MIN(basin2cluster[basin_i],
                                        basin2cluster[basin_j]);
          basin2cluster[basin_i] = merged_cluster_id;
          basin2cluster[basin_j] = merged_cluster_id;

        } // if (aaaiDest[iz_jz][iy_jy][ix_jx] != aaaiDest[iz][iy][ix])
      } // else clause for "if (aaaiDest[iz_jz][iy_jy][ix_jx] == UNDEFINED)"
    } // for (int j = 0; j < num_neighbors; j++)
  } //while (q.size() != 0)


  #ifndef NDEBUG
  // DEBUGGING
  // All of the voxels should either be assigned to a basin,
  // OR assigned to "WATERSHED_BOUNDARY" or to "UNDEFINED"
  // (but not "QUEUED").  Check for that below:
  for (int iz=0; iz<image_size[2]; iz++)
    for (int iy=0; iy<image_size[1]; iy++)
      for (int ix=0; ix<image_size[0]; ix++)
        assert(aaaiDest[iz][iy][ix] != QUEUED);
  #endif //#ifndef NDEBUG


  // Count the number of clusters
  size_t n_clusters = 0;
  size_t n_basins = basin2cluster.size();
  vector<size_t> clusterold2clusternew(n_basins);
  vector<size_t> cluster2basin; //inverse lookup table
  for (size_t i=0; i < n_basins; i++) {
    // (Originally, clusters were assigned to the ID number of the basin they
    //  started from.  But after merging multiple basins into the same cluster
    //  there are more basins than clusters.)
    clusterold2clusternew[i] = n_clusters;
    if (basin2cluster[i] == i) {
      cluster2basin.push_back(i);
      assert(cluster2basin[n_clusters] == i);
      n_clusters++;
    }
  }

  // Renumber the clusters from 0 ... n_clusters-1
  for (size_t i=0; i < basin2cluster.size(); i++) {
    basin2cluster[i] = clusterold2clusternew[basin2cluster[i]];
  }
    
  // Now, assign all the voxels in aaaiDest[][][]
  // to their clusters instead of their basins.
  for (int iz=0; iz<image_size[2]; iz++)
    for (int iy=0; iy<image_size[1]; iy++)
      for (int ix=0; ix<image_size[0]; ix++)
        aaaiDest[iz][iy][ix] = basin2cluster[ aaaiDest[iz][iy][ix] ];

  if (pv_cluster_centers != NULL) {
    pv_cluster_centers->resize(n_clusters);
    for (size_t i = 0; i < n_clusters; i++)
      (*pv_cluster_centers)[i] = extrema_locations[ cluster2basin[i] ];
  }

  // Now, deal with the clusters which the caller wants the "undefined"
  // voxels to have a high instead of a low value.
  if (undefined_label_is_max)
    label_undefined = n_clusters+1;

  for (int iz=0; iz<image_size[2]; iz++) {
    for (int iy=0; iy<image_size[1]; iy++) {
      for (int ix=0; ix<image_size[0]; ix++) {
        if (aaaiDest[iz][iy][ix] == UNDEFINED)
          aaaiDest[iz][iy][ix] == label_undefined;
      }
    }
  }

} // ConnectedClusters()






/// @brief  IGNORE THIS CLASS.  NOT CURRENTLY AVAILABLE
template<class Scalar, class Integer, class VectorContainer, class TensorContainer>
class TV3D_ACO : public TV3D<Scalar,Integer,VectorContainer,TensorContainer>
{
  /// @brief  IGNORE THIS FUNCTION.
  ///         This is a version of "Ant Colony Optimization"(ACO) which attempts
  ///         to exploit the directionality stored within voxels in an image.
  ///         (Such directionality can be generated by calculating the gradient,
  ///          or Hessian, or applying tensor-voting to the original image.)
  /// @note   This algorithm is unfinished as of 2019-1-23
  ///         (and will honestly probably never by finished).
  /// @note   As of 2019-1-23, there is no particular reason yet that this 
  ///         function needs to be a member function of "TV3D".  
  ///         For now at least, I could make it a stand-alone function.
  void
  ACO(Integer const image_size[3],  //!< source image size
      Scalar const *const *const *aaafSaliency,  //!< saliency (score) of each voxel (usually based on Hessian eigenvalues)
      VectorContainer const *const *const *aaaafV,  //!< vector associated with each voxel
      Scalar ***aaafPheremones,  //!< the pheremone density of every voxel is stored here
      Scalar const *const *const *aaafMask, //!< optional: if not NULL then we ignore voxel ix,iy,iz if aaafMask[iz][iy][ix]==0
      Integer num_ants, //!< number of independent "ant" agents in the simulation
      Scalar evaporation_rate, //!< the fraction of pheremones which evaporate at every simulation timestep (between 0 and 1.  Usually << 1)
      Scalar exponent_saliency = 1, //!< parameter used in ACO
      Scalar exponent_pheremone = 1, //!< parameter used in ACO
      Integer connectivity=3,        //!< square root of the search radius around each voxel (1=nearest_neighbors, 2=2D_diagonal, 3=3D_diagonal)
      //Scalar threshold_saliency = 0.0,
      //bool normalize=true,           //!< REMOVE THIS?  not clear if this parameter is still relevant
      Integer sim_duration=0,        //!< run the ACO simulation for this may iterations
      //Integer backup_interval=0,    //!< periodically back up the simulation
      //ostream *pReportBackups=NULL, //!< periodically back up the simulation
      ostream *pReportProgress=NULL   //!< print progress to the user?
      )
  {
    assert(image_size);
    assert(aaafSaliency);
    assert(aaaafV);
    assert(aaafPheremones);

    size_t n_voxels = 0;

    // Initialize the aaafPheremones[][][] array to be a constant everywhere
    for (int iz = 1; iz < image_size[2]-1; iz++) {
      for (int iy = 1; iy < image_size[1]-1; iy++) {
        for (int ix = 1; ix < image_size[0]-1; ix++) {
          aaafPheremones[iz][iy][ix] = 1.0;
          if ((! aaafMask) && (aaafMask[iz][iy][ix] == 0.0)) {
            // ants can never visit outside the mask
            // so their pheremone density there should be zero
            aaafPheremones[iz][iy][ix] = 0.0;
            continue;
          }
          n_voxels++;
        }
      }
    }

    // "delta_pheremones" is the amount of "pheremone" deposited on a voxel
    // each time an ant visits it.  What is a good value for this parameter?

    Scalar delta_pheremones = (evaporation_rate * n_voxels) / num_ants;

    // Hopefully setting "delta_pheremones" to this value insures the total
    // amount of pheremones in the image remains roughly constant over time.
    // (Why?  The total amount of pheremones that evaporate per iteration is
    //  n_voxels * evaporation rate.  In every iteration of the simulation, 
    //  each ant deponsits "delta_pheremones".  There are "num_ants" of them.)
    // Note: The total amount of pheremones might still drift over long
    //       simulation times due to numerical underflow.

    // Randomize the initial position of the ants
    Integer (*aaiAntCrds)[3] = NULL;
    aaiAntCrds = new int[num_ants][3];
    for (Integer i = 0; i < num_ants; i++) {
      bool inside = false;
      while (! inside) {
        int ix = RANDOM_INT(image_size[0]);
        int iy = RANDOM_INT(image_size[1]);
        int iz = RANDOM_INT(image_size[2]);
        if ((! aaafMask) || (aaafMask[iz][iy][ix] != 0.0))
          inside = true;
      }
    } //for (Integer i = 0; i < num_ants; i++)


    // In the next step, simulated "ants" will wander from one voxel
    // to its neighbors.  Before we continue, we should
    // figure out which neighbors to consider when searching neighboring voxels
    int (*neighbors)[3] = NULL; //a pointer to a fixed-length array of 3 ints
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
            //if ((jx == 0) && (jy == 0) && (jz == 0))
            //  continue;
            if (jx*jx+jy*jy+jz*jz > r_neigh_sqd)
              continue;
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

    Scalar *prob_neigh = new Scalar [num_neighbors];



    // Ant colony optimization (ACO) simulation method
    for (Integer t=0; t < sim_duration; t++) {
      // In each iteration, move each ant to a new location
      #pragma omp parallel
      for (Integer i = 0; i < num_ants; i++) {
        int ix = aaiAntCrds[i][0];
        int iy = aaiAntCrds[i][1];
        int iz = aaiAntCrds[i][2];

        // loop over all of the neighbors
        for (int j = 0; j < num_neighbors; j++) {
          int jx = neighbors[j][0];
          int jy = neighbors[j][1];
          int jz = neighbors[j][2];
          int iz_jz = iz + jz;
          int iy_jy = iy + jy;
          int ix_jx = ix + jx;

          // ignore voxels which lie outside the accessible portion of the image
          if (((iz_jz < 0) || (image_size[2] <= iz_jz))
              ||
              ((iy_jy < 0) || (image_size[1] <= iy_jy))
              ||
              ((ix_jx < 0) || (image_size[0] <= ix_jx))
              ||
              (aaafMask && (aaafMask[iz_jz][iy_jy][ix_jx] == 0.0)))
          {
            prob_neigh[j] = 0.0;
            continue;
          }

          // Calculate the probability that the ant visits this voxel:
          // FILL THIS IN LATER !!!
          //prob_neigh[j] = SOME_FUNCTION_OF_PHEREMONES_AND_SALIENCY_AND_V

        } //for (int j = 0; j < num_neighbors; j++)


        // normalize the probabilities (stored in "prob_neigh[]")
        Scalar tot_p = 0.0;
        for (int j = 0; j < num_neighbors; j++)
          tot_p += prob_neigh[j];
        for (int j = 0; j < num_neighbors; j++)
          prob_neigh[j] /= tot_p;

        // choose one of the neighboring voxels at random
        Scalar dice_roll = RANDOM_REAL_0_1();
        Scalar prob_so_far = 0.0;
        int j;
        for (j = 0; j < num_neighbors; j++) {
          prob_neigh[j] /= tot_p;
          prob_so_far += prob_neigh[j];
          if (prob_so_far > dice_roll)
            break;
        } //for (j = 0; j < num_neighbors; j++)
        int jx = neighbors[j][0];
        int jy = neighbors[j][1];
        int jz = neighbors[j][2];
        aaiAntCrds[i][0] = ix+jx;
        aaiAntCrds[i][1] = iy+jy;
        aaiAntCrds[i][2] = iz+jz;
      } //for (Integer i = 0; i < num_ants; i++)

      #pragma omp parallel for collapse(3)
      for (int iz = 1; iz < image_size[2]-1; iz++) {
        for (int iy = 1; iy < image_size[1]-1; iy++) {
          for (int ix = 1; ix < image_size[0]-1; ix++) {
            if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
              continue;
            aaafPheremones[iz][iy][ix] *= (1.0 - evaporation_rate);
          }
        }
      }
      #pragma omp parallel
      for (Integer i = 0; i < num_ants; i++) {
        int ix = aaiAntCrds[i][0];
        int iy = aaiAntCrds[i][1];
        int iz = aaiAntCrds[i][2];
        aaafPheremones[iz][iy][ix] += delta_pheremones;
      } //for (Integer i = 0; i < num_ants; i++)

    } //for (Integer t=0; t < sim_duration; t++)

    delete [] aaiAntCrds;//(note:this 2D array is actually a contiguous 1Darray)
    delete [] prob_neigh;

  } //ACO()


}; //class TV3D_ACO





#ifndef DISABLE_TEMPLATE_MATCHING


template<class Scalar, class Integer>
class TemplateMatcher3D : public Filter3D<Scalar, Integer>
{
public:

  TemplateMatcher3D():
    Filter3D<Scalar, Integer>()
  { }

  TemplateMatcher3D(Filter3D<Scalar, Integer> copy_from):
    Filter3D<Scalar, Integer>(copy_from)
  { }
  
  /// @brief  Compute the sum weighted difference between the local entries
  ///         in an array centered around position ix,iy,iz, 
  ///         and another array we call the "template" (which is assumed to
  ///         be stored in the "aaafH" filter array member for this class).
  ///         It is assumed that the weighted average of both the
  ///          source and template arrays is 0 before calling this function.
  /// @return the weighted sum of the error (raised to the power of
  ///         err_exponent, which is 2 by default) summed over all entries
  ///         in the template.  Weights are typically chosen so that they
  ///         decay to 0 near the boundaries of the aaafH array.
  ///         This function is not yet intended for public use.
  Scalar
  _TemplateError(Integer ix, //!< voxel's position in the x,y,z directions
                 Integer iy, //!< voxel's position in the x,y,z directions
                 Integer iz, //!< voxel's position in the x,y,z directions
                 Integer const size_source[3], //!< size of the source array
                 Scalar const *const *const *aaafSource,   //!< the array 
                 Scalar scale, //!< how much to scale template intensities
                 Scalar const *const *const *aaafW, //!< weights used in all summations
                 Scalar err_exponent=2, //!< exponent used for calculating error
                 Scalar const *const *const *aaafMask = NULL, //!< optional array storing voxels we should ignore (if 0)
                 ostream *pReportProgress = NULL //!< optional ostream for printing out progress of the calculation
                 ) const 
  {

    Scalar g = 0.0;
    Scalar denominator = 0.0;

    for (Integer
         jz=-Filter3D<Scalar, Integer>::halfwidth[2];
         jz<=Filter3D<Scalar, Integer>::halfwidth[2]; jz++)
    {

      Integer iz_jz = iz-jz;
      if ((iz_jz < 0) || (size_source[2] <= iz_jz))
        continue;

      for (Integer
           jy=-Filter3D<Scalar, Integer>::halfwidth[1];
           jy<=Filter3D<Scalar, Integer>::halfwidth[1]; jy++)
      {

        Integer iy_jy = iy-jy;
        if ((iy_jy < 0) || (size_source[1] <= iy_jy))
          continue;

        for (Integer
             jx=-Filter3D<Scalar, Integer>::halfwidth[0];
             jx<=Filter3D<Scalar, Integer>::halfwidth[0]; jx++)
        {

          Integer ix_jx = ix-jx;
          if ((ix_jx < 0) || (size_source[0] <= ix_jx))
            continue;

          Scalar delta_g = 
            (scale * Filter3D<Scalar, Integer>::aaafH[jz][jy][jx]
             -
             aaafSource[iz_jz][iy_jy][ix_jx]);

          if (err_exponent == 2.0)
            delta_g *= delta_g;
          else
            delta_g = pow(delta_g, err_exponent);

          //if (! precompute_mask_times_source)
          delta_g *= aaafW[jz][jy][jx];


          if (aaafMask)     //<--rarely used. I may delete "aaafMask" later
            delta_g *= aaafMask[iz_jz][iy_jy][ix_jx];

          g += delta_g;

          //if (normalize) {
          if (aaafMask)     //<--rarely used. I may delete "aaafMask" later
            denominator +=
              aaafMask[iz_jz][iy_jy][ix_jx]
              *
              aaafW[jz][jy][jx];

          else
            //denominator += 1.0;
            denominator += aaafW[jz][jy][jx];

          //} //if (normalize)

        }
      }
    }

    //if (normalize) {
    if (denominator > 0.0)
      g /= denominator;
    else
      g = 1.0;

    return g;
  } // _TemplateError()



  /// @brief Calculate the TemplateError() everywhere in the source array
  ///        (aaafSource).  Store the result in aaafDest.
  ///        This function was not intended for public use.

  void
  _ScanTemplateError(Integer const size_source[3], //!< size of the source array
                     Scalar const *const *const *aaafSource, //!< source array
                     Scalar ***aaafDest, //!< store the template error results here
                     Scalar ***aaafC, //!< how much to scale template intensities
                     Scalar const *const *const *aaafW, //!< weights used in all summations
                     Scalar err_exponent=2, //!< exponent used for calculating error
                     Scalar const *const *const *aaafMask = NULL, //!< optional: indicate which entries should be ignored
                     ostream *pReportProgress = NULL //!< optional: print out the progress of the calculation
                     ) const
  {
    if (pReportProgress)
      *pReportProgress << "  progress: processing plane#" << endl;

    for (Integer iz=0; iz<size_source[2]; iz++) {

      if (pReportProgress)
        *pReportProgress << "  " << iz+1 << " / " << size_source[2] << "\n";

      for (Integer iy=0; iy<size_source[1]; iy++) {

        for (Integer ix=0; ix<size_source[0]; ix++) {

          // Calculate the effect of the filter on
          // the voxel located at position ix,iy,iz

          if ((aaafMask) && (aaafMask[iz][iy][ix] == 0.0)) {
            aaafDest[iz][iy][ix] = 1.0;
            continue;
          }
          
          aaafDest[iz][iy][ix] =
            TemplateError(ix, iy, iz,
                          size_source,
                          aaafSource,
                          aaafC[iz][iy][ix],
                          aaafW,
                          err_exponent,
                          NULL); //aaafMask);
                          //normalize);
        }
      }
    }
  } //_ScanTemplateError()

}; // class TemplateMatcher::Filter3D





/// @brief
///  Perform template matching with a Gaussian.
///  DEPRECIATION WARNING: This function may be removed in the future.
///  A comment describing of this function is provided in "unsupported.cpp"

void
HandleTemplateGauss(Settings settings,
                    MrcSimple &tomo_in,
                    MrcSimple &tomo_out,
                    MrcSimple &mask,
                    float voxel_width[3]);

/// @brief
///  Perform template matching with a (spherically symmetric)
///  generalized Gaussian.
///  DEPRECIATION WARNING: This function may be removed in the future.
///  A comment describing of this function is provided in "unsupported.cpp"

void
HandleTemplateGGauss(Settings settings,
                     MrcSimple &tomo_in,
                     MrcSimple &tomo_out,
                     MrcSimple &mask,
                     float voxel_width[3]);

#endif //#ifndef DISABLE_TEMPLATE_MATCHING






#ifndef DISABLE_BOOTSTRAPPING

/// @brief   Scramble the contents of an image by replacing the current voxel
///          with a randmly chosen voxel nearby (which lies within an
///          ellipsoid specified by "scramble_radius")
template<class Scalar>
void
ScrambleImage3D(int image_size[3],
                Scalar const *const *const *aaafSource,
                  //(ellipsoidal) scramble radius in x,y,z directions:
                int const scramble_radius[3], 
                Scalar ***aaafDest)
{
  for (int iz=0; iz < image_size[2]; iz++) {
    for (int iy=0; iy < image_size[1]; iy++) {
      for (int ix=0; ix < image_size[0]; ix++) {
        Scalar r = 2.0;
        int dx, dy, dz;
        while (r > 1.0) {
          dx = RANDOM_INT(1+2*scramble_radius[0]) - scramble_radius[0];
          dy = RANDOM_INT(1+2*scramble_radius[1]) - scramble_radius[1];
          dz = RANDOM_INT(1+2*scramble_radius[2]) - scramble_radius[2];
          r = (SQR(static_cast<float>(dx)/scramble_radius[0]) +
               SQR(static_cast<float>(dy)/scramble_radius[1]) +
               SQR(static_cast<float>(dz)/scramble_radius[2]));
        }
        aaafDest[iz][iy][ix] =
          aaafSource[iz+dz][iy+dy][iz+dz];
      }
    }
  }
} //ScrambleImage3D()


/// @brief
/// DEPRECIATION WARNING: This function will probably be removed in the future.
void
HandleBootstrapDogg(Settings settings,
                    MrcSimple &tomo_in,
                    MrcSimple &tomo_out,
                    MrcSimple &mask);

#endif //#ifndef DISABLE_BOOTSTRAPPING




#ifndef DISABLE_DOGGXY
void
HandleDoggXY(Settings settings,
             MrcSimple &tomo_in,
             MrcSimple &tomo_out,
             MrcSimple &mask,
             float voxel_width[3]);
#endif //#ifndef DISABLE_DOGGXY



#endif //#ifndef _UNSUPPORTED_HPP
