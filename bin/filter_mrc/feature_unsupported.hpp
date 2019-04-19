#ifndef _FEATURE_UNSUPPORTED_HPP
#define _FEATURE_UNSUPPORTED_HPP

#include <cmath>
#include <visfd.hpp>
#include <random_gen.h>





namespace visfd {





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





/// @brief  IGNORE THIS CLASS.  NOT CURRENTLY AVAILABLE
template<typename Scalar, typename Integer, typename VectorContainer, typename TensorContainer>
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
      Scalar const *const *const *aaafMask, //!< optional: if not nullptr then we ignore voxel ix,iy,iz if aaafMask[iz][iy][ix]==0
      Integer num_ants, //!< number of independent "ant" agents in the simulation
      Scalar evaporation_rate, //!< the fraction of pheremones which evaporate at every simulation timestep (between 0 and 1.  Usually << 1)
      Scalar exponent_saliency = 1, //!< parameter used in ACO
      Scalar exponent_pheremone = 1, //!< parameter used in ACO
      Integer connectivity=3,        //!< square root of the search radius around each voxel (1=nearest_neighbors, 2=2D_diagonal, 3=3D_diagonal)
      //Scalar threshold_saliency = 0.0,
      //bool normalize=true,           //!< REMOVE THIS?  not clear if this parameter is still relevant
      Integer sim_duration=0,        //!< run the ACO simulation for this may iterations
      //Integer backup_interval=0,    //!< periodically back up the simulation
      //ostream *pReportBackups=nullptr, //!< periodically back up the simulation
      ostream *pReportProgress=nullptr   //!< print progress to the user?
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
    Integer (*aaiAntCrds)[3] = nullptr;
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


template<typename Scalar, typename Integer>
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
                 Scalar const *const *const *aaafMask = nullptr, //!< optional array storing voxels we should ignore (if 0)
                 ostream *pReportProgress = nullptr //!< optional ostream for printing out progress of the calculation
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
                     Scalar const *const *const *aaafMask = nullptr, //!< optional: indicate which entries should be ignored
                     ostream *pReportProgress = nullptr //!< optional: print out the progress of the calculation
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
                          nullptr); //aaafMask);
                          //normalize);
        }
      }
    }
  } //_ScanTemplateError()

}; // class TemplateMatcher::Filter3D

#endif //#ifndef DISABLE_TEMPLATE_MATCHING




#ifndef DISABLE_BOOTSTRAPPING

/// @brief   Scramble the contents of an image by replacing the current voxel
///          with a randmly chosen voxel nearby (which lies within an
///          ellipsoid specified by "scramble_radius")
template<typename Scalar>
void
ScrambleImage(int image_size[3],
              Scalar const *const *const *aaafSource,
                // (ellipsoidal) scramble radius in x,y,z directions:
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
} //ScrambleImage()

#endif //#ifndef DISABLE_BOOTSTRAPPING




#ifndef DISABLE_INTENSITY_PROFILES

typedef enum eBlobCenterCriteria {
  MAXIMA,
  MINIMA,
  CENTER
} BlobCenterCriteria;



/// @brief  Create a plot of intensity vs. radius for each blob.
///         (It's not clear this is of general use to people, so I will keep
///          this function out of the main library for now.  -andrew 2019-3-22)

template<typename Scalar>
void
BlobIntensityProfiles(int const image_size[3], //!< image size
                      Scalar const *const *const *aaafSource,   //!< ignore voxels where mask==0
                      Scalar const *const *const *aaafMask,   //!< ignore voxels where mask==0
                      const vector<array<Scalar,3> > &sphere_centers, //!< coordinates for the center of each sphere (blob)
                      const vector<Scalar> &diameters,         //!< diameter of each sphere (in voxels)
                      BlobCenterCriteria center_criteria,
                      vector<vector<Scalar> > &intensity_profiles //!< store the intensity profiles here
                      )
{
  assert(image_size);
  assert(aaafSource);

  size_t N = sphere_centers.size();
  intensity_profiles.resize(N);
  for (size_t i = 0; i < N; i++) {
    int Rsphere = ceil(diameters[i]/2); //radius of the sphere (surrounding the current blob)
    int ixs = floor(sphere_centers[i][0] + 0.5); //voxel coordinates of the
    int iys = floor(sphere_centers[i][1] + 0.5); //center of the current blob's
    int izs = floor(sphere_centers[i][2] + 0.5); //bounding sphere
    int ix0;  // coordinates of the "center" of each blob.  (Note: Depending
    int iy0;  // on "center_criteria", the "center" of the blob might not be 
    int iz0;  // located at the center of the sphere that encloses it.)
    if (center_criteria == BlobCenterCriteria::CENTER) {
      ix0 = ixs;
      iy0 = iys;
      iz0 = izs;
    }
    else {
      bool first_iter = true;
      Scalar extrema_val = 0.0;
      for (int jz = -Rsphere; jz <= Rsphere; jz++) {
        int izs_jz = izs + jz;
        for (int jy = -Rsphere; jy <= Rsphere; jy++) {
          int iys_jy = iys + jy;
          for (int jx = -Rsphere; jx <= Rsphere; jx++) {
            int ixs_jx = ixs + jx;
            if ((jx*jx + jy*jy + jz*jz) > Rsphere*Rsphere)
              continue;
            if (! ((0 <= ixs_jx) && (ixs_jx <= image_size[0]) &&
                   (0 <= iys_jy) && (iys_jy <= image_size[1]) &&
                   (0 <= izs_jz) && (izs_jz <= image_size[2])))
              continue;
            if (aaafMask && (aaafMask[izs_jz][iys_jy][ixs_jx] == 0.0))
              continue;
            if (first_iter
                ||
                ((center_criteria == BlobCenterCriteria::MAXIMA) && 
                 (aaafSource[izs_jz][iys_jy][ixs_jx] > extrema_val))
                ||
                ((center_criteria == BlobCenterCriteria::MINIMA) && 
                 (aaafSource[izs_jz][iys_jy][ixs_jx] < extrema_val)))
            {
              ix0 = ixs_jx;
              iy0 = iys_jy;
              iz0 = izs_jz;
              extrema_val = aaafSource[izs_jz][iys_jy][ixs_jx];
              first_iter = false;
            }
          }
        }
      }
    } //else clause for "if (center_criteria == BlobCenterCriteria::CENTER)"

    
    int Rsearch = ceil(Rsphere + 
                       sqrt(SQR(ix0-ixs) + SQR(iy0-iys) + SQR(iz0-izs)));
      
    intensity_profiles[i].resize(Rsearch+1);

    vector<Scalar> numerators(Rsearch+1, 0.0);
    vector<Scalar> denominators(Rsearch+1, 0.0);

    for (int jz = -Rsearch; jz <= Rsearch; jz++) {
      int iz0_jz = iz0 + jz;
      for (int jy = -Rsearch; jy <= Rsearch; jy++) {
        int iy0_jy = iy0 + jy;
        for (int jx = -Rsearch; jx <= Rsearch; jx++) {
          int ix0_jx = ix0 + jx;
          if ((jx*jx + jy*jy + jz*jz) > Rsearch*Rsearch)
            continue;
          if (! ((0 <= ix0_jx) && (ix0_jx <= image_size[0]) &&
                 (0 <= iy0_jy) && (iy0_jy <= image_size[1]) &&
                 (0 <= iz0_jz) && (iz0_jz <= image_size[2])))
            continue;
          if (aaafMask && (aaafMask[iz0_jz][iy0_jy][ix0_jx] == 0.0))
            continue;
          // Does this voxel lie within the sphere?
          //    (sphere is centered at ixc, iyc, izc, and has radius Rsphere)
          int Jx = jx + ix0 - ixs;
          int Jy = jy + iy0 - iys;
          int Jz = jz + iz0 - izs;
          int Jr = static_cast<int>(floor(sqrt(Jx*Jx + Jy*Jy + Jz*Jz) + 0.5));
          if (Jr > Rsphere)
            continue; //// If not, ignore it.
          int jr = static_cast<int>(floor(sqrt(jx*jx + jy*jy + jz*jz) + 0.5));
          // Consider all remaining voxels
          numerators[jr] += aaafSource[iz0_jz][iy0_jy][ix0_jx];
          denominators[jr] += 1.0;
        } //for (int jx = -Rsearch; jx <= Rsearch; jx++)
      } //for (int jy = -Rsearch; jy <= Rsearch; jy++)
    } //for (int jz = -Rsearch; jz <= Rsearch; jz++)
    for (int ir = 0; ir <= Rsearch; ir++) {
      if (denominators[ir] == 0.0) {
        intensity_profiles[i].resize(ir);
        break;
      }
      intensity_profiles[i][ir] = numerators[ir] / denominators[ir];
    }
  } //for (size_t i = 0; i < N; i++)
} //BlobIntensityProfiles()

#endif //#ifndef DISABLE_INTENSITY_PROFILES




} //namespace visfd




#endif //#ifndef _FEATURE_UNSUPPORTED_HPP
