///   @file feature.hpp
///   @brief a collection functions for detecting features in images
///   @author Andrew Jewett
///   @date 2019-4-15

#ifndef _FEATURE_HPP
#define _FEATURE_HPP

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
#include <visfd_utils.hpp>    // defines invert_permutation(), SortBlobs(), FindSpheres() ...
#include <alloc2d.hpp>    // defines Alloc2D() and Dealloc2D()
#include <alloc3d.hpp>    // defines Alloc3D() and Dealloc3D()
#include <filter1d.hpp>   // defines "Filter1D" (used in ApplySeparable())
#include <filter3d.hpp>   // defines common 3D image filters
#include <feature_implementation.hpp>



namespace visfd {



/// @brief   Find all of the local minima and local maxima in an image (aaafI).
///          The locations of minima and maxima are stored in the
///          *pv_minima_crds and *pv_maxima_crds arrays, and sorted in
///          increasing and decreasing order respectively.  (IE they are sorted
///          so that the most significant local minima or maxima will appear
///          first in the list.)
///          If either pv_minima_crds or pv_maxima_crds is nullptr, then
///          the minima or maxima will be ignored.
///          The optional pv_minima_scores and pv_maxima_scores store the
///          The caller can automatically discard minima or maxima which
///          are not sufficiently low or high, by supplying thresholds.
///          The optional aaafMask array (if not nullptr) can be used to ignore
///          certain voxels in the image (whose aaafMask entries are zero).
/// @note    Local minima or maxima on the boundaries of the image 
///          (or near the edge of the mask)
///          are not as trustworthy since some of the neighboring voxels
///          will not be available for comparison.  These minima and maxima
///          can be ignored by setting allow_borders=false.  The number of
///          neighbors around every voxel which are considered (eg, 6, 18, 26)
///          can be controlled using the "connectivity" argument.

template<typename Scalar, typename Coordinate, typename IntegerIndex, typename Label>

void
FindExtrema(int const image_size[3],          //!< size of the image in x,y,z directions
            Scalar const *const *const *aaafI,    //!< image array aaafI[iz][iy][ix]
            Scalar const *const *const *aaafMask, //!< optional: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
            vector<array<Coordinate, 3> > *pva_minima_crds, //!< store minima locations (ix,iy,iz) here (if not nullptr)
            vector<array<Coordinate, 3> > *pva_maxima_crds, //!< store maxima locations (ix,iy,iz) here (if not nullptr)
            vector<Scalar> *pv_minima_scores, //!< store corresponding minima aaafI[iz][iy][ix] values here (if not nullptr)
            vector<Scalar> *pv_maxima_scores, //!< store corresponding maxima aaafI[iz][iy][ix] values here (if not nullptr)
            vector<IntegerIndex> *pv_minima_nvoxels, //!< store number of voxels in each minima (usually 1)
            vector<IntegerIndex> *pv_maxima_nvoxels, //!< store number of voxels in each maxima (usually 1)
            Scalar minima_threshold=std::numeric_limits<Scalar>::infinity(), //!< ignore minima with intensities greater than this
            Scalar maxima_threshold=-std::numeric_limits<Scalar>::infinity(), //!< ignore maxima with intensities lessr than this
            int connectivity=3,       //!< square root of search radius around each voxel (1=nearest_neighbors, 2=diagonal2D, 3=diagonal3D)
            bool allow_borders=true,  //!< if true, plateaus that touch the image border (or mask boundary) are valid extrema
            Label ***aaaiDest=nullptr,  //!< optional: create an image showing where the extrema are?
            ostream *pReportProgress=nullptr   //!< optional: print progress to the user?
                   )
{
  bool find_minima = (pva_minima_crds != nullptr);
  bool find_maxima = (pva_maxima_crds != nullptr);

  vector<size_t> minima_indices;
  vector<size_t> maxima_indices;
  vector<Scalar> minima_scores;
  vector<Scalar> maxima_scores;
  vector<size_t> *pv_minima_indices = nullptr;
  vector<size_t> *pv_maxima_indices = nullptr;
  if (find_minima) {
    pv_minima_indices = &minima_indices;
    if (! pv_minima_scores)
      pv_minima_scores = &minima_scores;
  }
  if (find_maxima) {
    pv_maxima_indices = &maxima_indices;
    if (! pv_maxima_scores)
      pv_maxima_scores = &maxima_scores;
  }

  _FindExtrema(image_size,
               aaafI,
               aaafMask,
               pv_minima_indices,
               pv_maxima_indices,
               pv_minima_scores,
               pv_maxima_scores,
               pv_minima_nvoxels,
               pv_maxima_nvoxels,
               minima_threshold,
               maxima_threshold,
               connectivity,
               allow_borders,
               aaaiDest,
               pReportProgress);

  // Now convert the indices back to x,y,z coordinates
  if (pva_minima_crds) {
    size_t N = minima_indices.size();
    assert(pva_minima_crds);
    pva_minima_crds->resize(N);
    for (size_t n = 0; n < N; n++) {
      size_t i = minima_indices[n];
      // convert from a 1D index (i) to 3-D indices (ix, iy, iz)
      size_t ix = i % image_size[0];
      i /= image_size[0];
      size_t iy = i % image_size[1];
      i /= image_size[1];
      size_t iz = i;
      (*pva_minima_crds)[n][0] = ix;
      (*pva_minima_crds)[n][1] = iy;
      (*pva_minima_crds)[n][2] = iz;
    }
  }
  if (pva_maxima_crds) {
    size_t N = maxima_indices.size();
    assert(pva_maxima_crds);
    pva_maxima_crds->resize(N);
    for (size_t n = 0; n < N; n++) {
      size_t i = maxima_indices[n];
      // convert from a 1D index (i) to 3-D indices (ix, iy, iz)
      size_t ix = i % image_size[0];
      i /= image_size[0];
      size_t iy = i % image_size[1];
      i /= image_size[1];
      size_t iz = i;
      (*pva_maxima_crds)[n][0] = ix;
      (*pva_maxima_crds)[n][1] = iy;
      (*pva_maxima_crds)[n][2] = iz;
    }
  }
} //FindExtrema()





/// @brief  The following version of this function seeks either minima
///         or maxima, but not both.  (If you want both, use the other version.
///         That version is faster than invoking this function twice.)
///         See the description of that version for details.

template<typename Scalar, typename Coordinate, typename IntegerIndex, typename Label>

void
FindExtrema(int const image_size[3],            //!< size of input image array
            Scalar const *const *const *aaafI,   //!< input image array
            Scalar const *const *const *aaafMask, //!< if not nullptr then zero entries indicate which voxels to ignore
            vector<array<Coordinate, 3> > &extrema_crds, //!< store the location of each extrema
            vector<Scalar> &extrema_scores, //!< store the brightness of each extrema
            vector<IntegerIndex> &extrema_nvoxels, //!< store number of voxels in each extrema (usually 1)
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
  vector<array<Coordinate, 3> > *null_vai3 = nullptr;  
  vector<Scalar> *null_vf = nullptr;  
  vector<IntegerIndex> *null_vi = nullptr;  

  if (seek_minima) {
    FindExtrema(image_size,
                aaafI,
                aaafMask,
                &extrema_crds,    // store minima locations here
                null_vai3,        // <-- don't search for maxima
                &extrema_scores,  // store minima values here
                null_vf,          // <-- don't search for maxima
                &extrema_nvoxels, // store number of voxels in each minima
                null_vi,          // <-- don't search for maxima
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
    FindExtrema(image_size,
                aaafI,
                aaafMask,
                null_vai3,        // <-- don't search for minima
                &extrema_crds,    // store maxima locations here
                null_vf,          // <-- don't search for minima
                &extrema_scores,  // store maxima values here
                null_vi,          // <-- don't search for minima
                &extrema_nvoxels, // store number of voxels in each maxima
                std::numeric_limits<Scalar>::infinity(),
                threshold,
                connectivity,
                allow_borders,
                aaaiDest,
                pReportProgress);
  }
} // FindExtrema()






/// @brief Find all scale-invariant blobs in the image as a function of Gaussian
///        width (ie. the "sigma" parameter), regardless of overlap.  (A blob's
///        "sigma" parameter represents the optimal width of the Gaussian blur
///        that, when applied to the image, maximizes this blob's visibility.)
///
/// Algorithm described in:
///    Lindeberg,T., Int. J. Comput. Vision., 30(2):77-116, (1998)
///
/// @note  The interpretation of the "sigma" parameter is not always straight-
///        forward.  An alternate version of this function exists ("BlobDogD()")
///        which approximates each blob with a solid sphere, and associates the
///        blob with the diameter of that sphere instead of its "sigma" value.
///
/// @note  This function can find BOTH minima and maxima without additional cost
///        (corresponding to dark blobs on a light background
///         and light blobs on a dark background, respectively).

template<typename Scalar>

void
BlobDog(int const image_size[3], //!< source image size
        Scalar const *const *const *aaafSource,   //!< source image
        Scalar const *const *const *aaafMask,     //!< ignore voxels where mask==0
        const vector<Scalar>& blob_sigma, //!< blob widths to try, ordered
        // optional arguments
        vector<array<Scalar,3> > *pva_minima_crds=nullptr, //!< if not nullptr, stores blob minima x,y,z coords here
        vector<array<Scalar,3> > *pva_maxima_crds=nullptr, //!< if not nullptr, stores blob maxima x,y,z coords here
        vector<Scalar> *pv_minima_sigma=nullptr, //!< if not nullptr, stores the corresponding width for that minima
        vector<Scalar> *pv_maxima_sigma=nullptr, //!< if not nullptr, stores the corresponding width for that maxima
        vector<Scalar> *pv_minima_scores=nullptr, //!< if not nullptr, stores the blob's score?
        vector<Scalar> *pv_maxima_scores=nullptr, //!< (score = intensity after filtering)
        //the following optional parameters are usually left with default values
        Scalar delta_sigma_over_sigma=0.02,//!< δ param for approximating LoG with DoG
        Scalar truncate_ratio=2.8,      //!< how many sigma before truncating?
        Scalar minima_threshold=std::numeric_limits<Scalar>::infinity(), //!< discard blobs with unremarkable scores (disabled by default)
        Scalar maxima_threshold=-std::numeric_limits<Scalar>::infinity(), //!< discard blobs with unremarkable scores (disabled by default)
        bool use_threshold_ratios=true, //!< threshold=ratio*best_score ?
        ostream *pReportProgress = nullptr, //!< optional: report progress to the user?
        Scalar ****aaaafI = nullptr, //!<optional: preallocated memory for filtered images (indexable)
        Scalar **aafI = nullptr      //!<optional: preallocated memory for filtered images (contiguous)
        )

{

  vector<array<Scalar,3> > minima_crds; //store minima blob x,y,z coords here
  vector<array<Scalar,3> > maxima_crds; //store minima blob x,y,z coords here
  vector<Scalar>  minima_sigma;         //corresponding width for that minima
  vector<Scalar>  maxima_sigma;         //corresponding width for that maxima
  vector<Scalar>  minima_scores;        //store the score of each blob minima
  vector<Scalar>  maxima_scores;        //store the score of each blob maxima
  if (pva_minima_crds == nullptr)
    pva_minima_crds = &minima_crds;
  if (pva_maxima_crds == nullptr)
    pva_maxima_crds = &maxima_crds;
  if (pv_minima_sigma == nullptr)
    pv_minima_sigma = &minima_sigma;
  if (pv_maxima_sigma == nullptr)
    pv_maxima_sigma = &maxima_sigma;
  if (pv_minima_scores == nullptr)
    pv_minima_scores = &minima_scores;
  if (pv_maxima_scores == nullptr)
    pv_maxima_scores = &maxima_scores;

  // We need 3 images to store the result of filtering the image
  // using DoG filters with 3 different widths.  Store those images here:

  if (pReportProgress)
    *pReportProgress
      << " -- Attempting to allocate space for 3 more images.        --\n"
      << " -- (If this crashes your computer, find a computer with   --\n"
      << " --  more RAM and use \"ulimit\", OR use a smaller image.)   --\n";

  bool preallocated = ! (aaaafI == nullptr);

  if (! preallocated) {
    assert(aaaafI == nullptr);
    assert(aafI == nullptr);
    aaaafI = new Scalar*** [3];
    aafI = new Scalar* [3];
    Alloc3D(image_size,
            &(aafI[0]),
            &(aaaafI[0]));
    Alloc3D(image_size,
            &(aafI[1]),
            &(aaaafI[1]));
    Alloc3D(image_size,
            &(aafI[2]),
            &(aaaafI[2]));
  }


  Scalar global_min_score = 1.0;  //impossible, minima < 0
  Scalar global_max_score = -1.0; //impossible, maxima > 0

  //REMOVE THIS CRUFT:
  //bool disable_thresholds = ((! use_threshold_ratios) &&
  //                           (minima_threshold >= maxima_threshold));

  if (pReportProgress)
    *pReportProgress << "\n----- Blob detection initiated using "
                     << blob_sigma.size() << " trial Gaussians -----\n\n";

  for (int ir = 0; ir < blob_sigma.size(); ir++) {

    if (pReportProgress)
      *pReportProgress
        << "--- Progress: "<<ir+1<<"/"<<blob_sigma.size() << "\n"
        << "--- Applying DoG filter using sigma[" << ir << "] = "
        << blob_sigma[ir] << " (in voxels) ---\n";

    // Run the image through a DoG filter with the currently selected width.
    // Compare the resulting filtered image with filtered images
    // smaller blob widths.
    // If a given voxel has a value which is larger than it's surrounding 26
    // voxels in this image, -as well as- the same 27 voxels which were filtered
    // using larger and smaller blob widths, respectively...
    //   ... -then- it is a local maximum in 4-dimensional X,Y,Z,R space.
    // Keep track of all of these X,Y,Z,R local maxima voxels and their values
    // (and the local minima as well).

    // We are going to apply a DoG filter to the original image and
    // store it in the aaaaf array.
    // Unfortunately, there probably won't be enough memory 
    // to keep track of all of the filtered versions of the image.
    // Instead we only need to keep track of the last 3 filtered images,
    // So the aaaafI[] array only has 3 entries.  We store the images in:
    //   aaaafI[ j2i[-1]]
    //   aaaafI[ j2i[0] ]
    //   aaaafI[ j2i[1] ]
    // ... and the values of j2i[-1], j2i[0], and j2i[1] cycle between 0,1,2

    int j2i_[3];       // an array of indices.  index into it using 0,1,2
    int *j2i = j2i_+1; // convenient to offset by 1 so we can index using -1,0,1
    j2i[-1] = (ir-2) % 3;
    j2i[0]  = (ir-1) % 3;
    j2i[1]  = ir % 3;

    //Apply the LoG filter (approximated with a DoG filter)
    //      ...using the most recent blob width:
    ApplyLog(image_size,
             aaafSource,
             aaaafI[j2i[1]], //<-store the most recent filtered image
             aaafMask,
             blob_sigma[ir],
             delta_sigma_over_sigma,
             truncate_ratio);

    // We must have filtered at least 3 images using different blob widths
    if (ir < 2)
      continue;

    // The blob widths should be ordered consistently (increasing or decreasing)
    assert(((blob_sigma[ir-2] < blob_sigma[ir-1]) &&   //increasing order
            (blob_sigma[ir-1] < blob_sigma[ir])) ||    
           ((blob_sigma[ir-2] > blob_sigma[ir-1]) &&   //decreasing order
            (blob_sigma[ir-1] > blob_sigma[ir])));


    if (pReportProgress)
      *pReportProgress
        << "--- Searching for local minima & maxima with sigma["<<ir-1<<"] = "
        << blob_sigma[ir-1] << " ---\n";

       //<< "--- Searching for local minima & maxima with width["<<ir-1<<"] = "
       //<< blob_sigma[ir-1]*2.0*sqrt(3) << " ---\n";



    // As we are looking for local minima and maxima we need to
    // keep track of the best scores so far
    // and discard blobs with weak scores as early as possible.
    // If use_threshold_ratios=true (default), then any
    // local maxima which are not > maxima_threshold * highest_score_so_far
    // will be discarded as are minima which are not
    // < maxima_threshold * lowest_score_so_far, which should be < 0.
    // In practice this keep the list of minima and maxima of manageable size.
    // (Otherwise a significant fraction of the voxels in the image could end
    // up being a minima or maxima at some scale, and this could easily use up
    // an enormous amout of memory and be difficult to deal with later.
    // It's better to prevent this list from growing too large at any point.

    #pragma omp parallel
    {

      // (The following variables are private for each thread/processor.
      // Later, the global list of minima and maxima will be updated with the
      // information collected from each processor at the end of this iteration)

      vector<array<Scalar,3> > min_crds_proc; //store minima x,y,z coords here
      vector<array<Scalar,3> > max_crds_proc; //store maxima x,y,z coords here
      vector<Scalar> min_sigma_proc; //corresponding width for that minima
      vector<Scalar> max_sigma_proc; //corresponding width for that maxima
      vector<Scalar> min_scores_proc; //what was the blob's score?
      vector<Scalar> max_scores_proc; //(score = intensity after filtering)
      Scalar global_min_score_proc = global_min_score;
      Scalar global_max_score_proc = global_max_score;

      #pragma omp for collapse(2)
      for (int iz = 0; iz < image_size[2]; iz++) {
        for (int iy = 0; iy < image_size[1]; iy++) {
          for (int ix = 0; ix < image_size[0]; ix++) {
            // Search the 81-1 = 80 surrounding voxels to see if this voxel is
            // either a minima or a maxima in 4-dimensional x,y,z,r space
            bool is_minima = true;
            bool is_maxima = true;
            for (int jr = -1; jr <= 1; jr++) {
              for (int jz = -1; jz <= 1; jz++) {
                for (int jy = -1; jy <= 1; jy++) {
                  for (int jx = -1; jx <= 1; jx++) {
                    if (jx==0 && jy==0 && jz==0 && jr==0)
                      continue; // Skip the central voxel. Check neighbors only
                    int Ix = ix + jx;
                    int Iy = iy + jy;
                    int Iz = iz + jz;
                    // All neighbors must be available for comparison
                    if ((Ix < 0) || (Ix >= image_size[0]) ||
                        (Iy < 0) || (Iy >= image_size[1]) ||
                        (Iz < 0) || (Iz >= image_size[2]) ||
                        (aaafMask && (aaafMask[Iz][Iy][Ix] == 0))) {
                      is_minima = false;
                      is_maxima = false;
                      continue;
                    }
                    Scalar entry    = aaaafI[ j2i[0]  ][ iz ][ iy ][ ix ];
                    Scalar neighbor = aaaafI[ j2i[jr] ][ Iz ][ Iy ][ Ix ];
                    if (neighbor <= entry)
                      is_minima = false;
                    if (neighbor >= entry)
                      is_maxima = false;
                  }
                }
              }
            }

            Scalar score = aaaafI[ j2i[0] ][ iz ][ iy ][ ix ];

            if ((! aaafMask) || (aaafMask[iz][iy][ix] != 0))
            {
              Scalar minima_threshold_so_far = minima_threshold;
              if (use_threshold_ratios)
                minima_threshold_so_far=minima_threshold*global_min_score_proc;
              if (is_minima &&
                  (score < 0.0) &&
                  ((score < minima_threshold_so_far) // || disable_thresholds
                   ))   
              {
                array<Scalar, 3> ixiyiz;
                ixiyiz[0] = ix;
                ixiyiz[1] = iy;
                ixiyiz[2] = iz;
                min_crds_proc.push_back(ixiyiz);
                min_sigma_proc.push_back(blob_sigma[ir-1]);
                min_scores_proc.push_back(score);
                if (score < global_min_score_proc)
                  global_min_score_proc = score;
              }

              Scalar maxima_threshold_so_far = maxima_threshold;
              if (use_threshold_ratios)
                maxima_threshold_so_far=maxima_threshold*global_max_score_proc;
              if (is_maxima &&
                  (score > 0.0) &&
                  ((score > maxima_threshold_so_far) // || disable_thresholds
                   ))   
              {
                array<Scalar, 3> ixiyiz;
                ixiyiz[0] = ix;
                ixiyiz[1] = iy;
                ixiyiz[2] = iz;
                max_crds_proc.push_back(ixiyiz);
                max_sigma_proc.push_back(blob_sigma[ir-1]);
                max_scores_proc.push_back(score);
                if (score > global_max_score_proc)
                  global_max_score_proc = score;
              }
            }
            assert(! (is_minima && is_maxima));
          } //for (int ix=0; ix<image_size[0]; ix++) {
        } //for (int iy=0; iy<image_size[1]; iy++) {
      } //for (int iz=0; iz<image_size[2]; iz++) {

      #pragma omp critical
      {

        // Append the newly minima and maxima discovered by this processor
        // to the global list of minima and maxima:

        pva_minima_crds->insert(pva_minima_crds->end(),
                                min_crds_proc.begin(),
                                min_crds_proc.end());
        min_crds_proc.clear();
        pv_minima_sigma->insert(pv_minima_sigma->end(),
                                min_sigma_proc.begin(),
                                min_sigma_proc.end());
        min_sigma_proc.clear();
        pv_minima_scores->insert(pv_minima_scores->end(),
                                 min_scores_proc.begin(),
                                 min_scores_proc.end());
        min_scores_proc.clear();
        pva_maxima_crds->insert(pva_maxima_crds->end(),
                                max_crds_proc.begin(),
                                max_crds_proc.end());
        max_crds_proc.clear();
        pv_maxima_sigma->insert(pv_maxima_sigma->end(),
                                max_sigma_proc.begin(),
                                max_sigma_proc.end());
        max_sigma_proc.clear();
        pv_maxima_scores->insert(pv_maxima_scores->end(),
                                 max_scores_proc.begin(),
                                 max_scores_proc.end());
        max_scores_proc.clear();

        // Update the global minima and maxima as well:
        if (global_min_score > global_min_score_proc)
          global_min_score = global_min_score_proc;
        if (global_max_score < global_max_score_proc)
          global_max_score = global_max_score_proc;
      } //#pragma omp critical

    } //#pragma omp parallel

    if (pReportProgress)
      *pReportProgress
        << "--- (Found " << pva_minima_crds->size ()
        << " and " << pva_maxima_crds->size()
        << " local minima and maxima, respectively so far) ---\n" << endl;

    assert((pva_minima_crds->size() == pv_minima_sigma->size()) &&
           (pva_minima_crds->size() == pv_minima_scores->size()));
    assert((pva_maxima_crds->size() == pv_maxima_sigma->size()) &&
           (pva_maxima_crds->size() == pv_maxima_scores->size()));

  } //for (ir = 0; ir < blob_sigma.size(); ir++)



  if ((minima_threshold != std::numeric_limits<Scalar>::infinity()) ||
      (maxima_threshold != -std::numeric_limits<Scalar>::infinity()))
  {
    if (pReportProgress)
      *pReportProgress
        << " Discarding poor scoring blobs...\n"
        << endl;

    if (use_threshold_ratios) {
      minima_threshold *= global_min_score;
      maxima_threshold *= global_max_score;
    }

    // Now that we know what the true global minima and maxima are,
    // go back and discard maxima whose scores are not higher than
    // maxima_threshold * global_max_score.
    // (Do the same for local minima as well.)

    { // throw away the poor-scoring minima
      assert(pva_minima_crds && pv_minima_sigma && pv_minima_scores);
      vector<array<Scalar,3> > minima_crds_cpy;
      vector<Scalar> minima_sigma_cpy;
      vector<Scalar> minima_scores_cpy;
      for (int i = 0; i < pv_minima_scores->size(); i++) {
        assert((*pv_minima_scores)[i] < 0.0);
        if ((*pv_minima_scores)[i] <= minima_threshold) {
          minima_crds_cpy.push_back((*pva_minima_crds)[i]);
          minima_sigma_cpy.push_back((*pv_minima_sigma)[i]);
          minima_scores_cpy.push_back((*pv_minima_scores)[i]);
        }
      }
      *pva_minima_crds = minima_crds_cpy;
      *pv_minima_sigma = minima_sigma_cpy;
      *pv_minima_scores = minima_scores_cpy;
    } // throw away the poor-scoring minima

    { // throw away the poor-scoring maxima
      assert(pva_maxima_crds && pv_maxima_sigma && pv_maxima_scores);
      vector<array<Scalar,3> > maxima_crds_cpy;
      vector<Scalar> maxima_sigma_cpy;
      vector<Scalar> maxima_scores_cpy;
      for (int i = 0; i < pv_maxima_scores->size(); i++) {
        assert((*pv_maxima_scores)[i] > 0.0);
        if ((*pv_maxima_scores)[i] >= maxima_threshold) {
          maxima_crds_cpy.push_back((*pva_maxima_crds)[i]);
          maxima_sigma_cpy.push_back((*pv_maxima_sigma)[i]);
          maxima_scores_cpy.push_back((*pv_maxima_scores)[i]);
        }
      }
      *pva_maxima_crds = maxima_crds_cpy;
      *pv_maxima_sigma = maxima_sigma_cpy;
      *pv_maxima_scores = maxima_scores_cpy;
    } // throw away the poor-scoring maxima

    if (pReportProgress)
      *pReportProgress << " ...done.\n" << endl;
  } //if ((minima_threshold != -std::numeric_limits<Scalar>::infinity()) || ...


  if (! preallocated) {
    // Deallocate the temporary arrays we created earlier
    Dealloc3D(image_size,
              &aafI[0],
              &aaaafI[0]);
    Dealloc3D(image_size,
              &aafI[1],
              &aaaafI[1]);
    Dealloc3D(image_size,
              &aafI[2],
              &aaaafI[2]);
    delete [] aaaafI;
    delete [] aafI;
  }

} //BlobDog()




/// @brief Find all scale-invariant blobs in the image as a function of their
///        (effective) diameters.  For conenience, each blob's diameter is 
///        reported instead of its "sigma" parameter.  Diameter is easy
///        to interpret.  (For comparison a blob's "sigma" parameter
///        represents the optimal width of the Gaussian blur that, when
//         applied to the image maximizes this blob's visibility.)
///
/// Algorithm described in:
///    Lindeberg,T., Int. J. Comput. Vision., 30(2):77-116, (1998)
///
/// @note  This function will find BOTH minima and maxima.
///        (corresponding to dark blobs on a light background
///         and light blobs on a dark background, respectively).

template<typename Scalar>

void
BlobDogD(int const image_size[3], //!<source image size
         Scalar const *const *const *aaafSource,   //!< source image
         Scalar const *const *const *aaafMask,     //!< ignore voxels where mask==0
         const vector<Scalar>& blob_diameters, //!< blob widths to try, ordered
         //optional arguments:
         vector<array<Scalar,3> > *pva_minima_crds=nullptr, //!< if not nullptr, stores blob minima x,y,z coords here
         vector<array<Scalar,3> > *pva_maxima_crds=nullptr, //!< if not nullptr, stores blob maxima x,y,z coords here
         vector<Scalar> *pv_minima_diameters=nullptr, //!< if not nullptr, stores the corresponding width for that minima
         vector<Scalar> *pv_maxima_diameters=nullptr, //!< if not nullptr, stores the corresponding width for that maxima
         vector<Scalar> *pv_minima_scores=nullptr, //!< if not nullptr, stores the blob's score?
         vector<Scalar> *pv_maxima_scores=nullptr, //!< (score = intensity after filtering)
         //the following optional parameters are usually left with default values
         Scalar delta_sigma_over_sigma=0.02,//!<param for approximating LoG with DoG
         Scalar truncate_ratio=2.5,    //!<how many sigma before truncating?
         Scalar minima_threshold=std::numeric_limits<Scalar>::infinity(), //!< ignore minima with intensities greater than this
         Scalar maxima_threshold=-std::numeric_limits<Scalar>::infinity(), //!< ignore maxima with intensities lessr than this
         bool    use_threshold_ratios=false, //!<threshold=ratio*best_score?
         ostream *pReportProgress = nullptr, //!<report progress to the user?
         Scalar ****aaaafI = nullptr, //!<preallocated memory for filtered images
         Scalar **aafI = nullptr     //!<preallocated memory for filtered images (conserve memory)
         )
{

  vector<Scalar> minima_sigma;
  vector<Scalar> maxima_sigma;
  vector<Scalar> blob_sigma(blob_diameters.size());
  for (int i=0; i < blob_diameters.size(); i++)
    blob_sigma[i] = blob_diameters[i]/(2.0*sqrt(3));

  BlobDog(image_size,
          aaafSource,
          aaafMask,
          blob_sigma,
          pva_minima_crds,
          pva_maxima_crds,
          &minima_sigma,
          &maxima_sigma,
          pv_minima_scores,
          pv_maxima_scores,
          delta_sigma_over_sigma,
          truncate_ratio,
          minima_threshold,
          maxima_threshold,
          use_threshold_ratios,
          pReportProgress,
          aaaafI,
          aafI);

  if (pv_minima_diameters) {
    pv_minima_diameters->resize(minima_sigma.size());
    for (int i=0; i < minima_sigma.size(); i++)
      (*pv_minima_diameters)[i] = minima_sigma[i] * 2.0 * sqrt(3);
  }

  if (pv_maxima_diameters) {
    pv_maxima_diameters->resize(maxima_sigma.size());
    for (int i=0; i < maxima_sigma.size(); i++)
      (*pv_maxima_diameters)[i] = maxima_sigma[i] * 2.0 * sqrt(3);
  }

} // BlobDogD()







/// @brief  Calculate the volume of overlap between two spheres of radius 
///         Ri and Rj separated by a distance of rij.

template<typename Scalar>

Scalar CalcSphereOverlap(Scalar rij,//!<the distance between the spheres' centers
                         Scalar Ri, //!< the radius of sphere i
                         Scalar Rj  //!< the radius of sphere j
                         )
{
  // WLOG, assume Ri <= Rj.  Insure that below
  if (Ri > Rj) {
    Scalar tmp = Ri;
    Ri = Rj;
    Rj = tmp;
  }
  if (rij <= Ri) {
    return (4*M_PI/3) * Ri*Ri*Ri;
  }

  // "xi" and "xj" are the distances from the sphere centers
  // to the plane where the two spheres intersect.
  Scalar xi = 0.5 * (1.0/rij) * (rij*rij + Ri*Ri - Rj*Rj);
  Scalar xj = 0.5 * (1.0/rij) * (rij*rij + Rj*Rj - Ri*Ri);
  Scalar volume_overlap =
    (M_PI/3)*( Ri*Ri*Ri * (2 - (xi/Ri) * (3 - SQR(xi/Ri))) +
               Rj*Rj*Rj * (2 - (xj/Rj) * (3 - SQR(xj/Rj))) );
  return volume_overlap;
}





/// @brief nonmax suppression for blobs
///
///   Blobs can either be discarded because the distance between their centers
///   exceeds the sum of their (*min_radial_separation_ratio)
///   OR because the overlapping volume exceeds the volume of the larger sphere
///   (after multiplication by max_volume_overlap_large)
///   OR because the overlapping volume exceeds the volume of the smaller sphere
///   (after multiplication by max_volume_overlap_small)
///
/// @note
///   The "sort_blob_criteria" argument can be set to one of these choices:
///     PRIORITIZE_HIGH_SCORES
///     PRIORITIZE_LOW_SCORES
///     PRIORITIZE_HIGH_MAGNITUDE_SCORES
///     PRIORITIZE_LOW_MAGNITUDE_SCORES
///   (As of 2019-6-15, "SortBlobCriteria" is defined in "visfd_utils.hpp")

#include <feature_implementation.hpp> // defines _FindExtrema()

template<typename Scalar>

void
DiscardOverlappingBlobs(vector<array<Scalar,3> >& blob_crds, //!< location of each blob
                        vector<Scalar>& blob_diameters,  //!< diameger of each blob
                        vector<Scalar>& blob_scores, //!< priority of each blob
                        Scalar min_radial_separation_ratio, //!< discard blobs if closer than this (ratio of sum of radii)
                        Scalar max_volume_overlap_large=std::numeric_limits<Scalar>::infinity(), //!< discard blobs which overlap too much with the large blob (disabled by default; 1.0 would also do this)
                        Scalar max_volume_overlap_small=std::numeric_limits<Scalar>::infinity(), //!< discard blobs which overlap too much with the small blob (disabled by default; 1.0 would also do this)
                        SortBlobCriteria sort_blob_criteria=PRIORITIZE_HIGH_MAGNITUDE_SCORES, //!< give priority to high or low scoring blobs? (See explanation above)
                        ostream *pReportProgress = nullptr, //!< report progress back to the user?
                        int scale=6 //!<occupancy_table_size shrunk by this much
                                    //!<relative to source (necessary to reduce memory usage)
                        )
{
  assert(blob_crds.size() == blob_diameters.size());
  assert(blob_crds.size() == blob_scores.size());

  // Strategy:
  // 1) Sort the blobs in order of their scores: from good scores, to bad scores
  // 2) Beginning with the best scoring blobs and proceeding downwards,
  //    check to see if the current blob overlaps with any of the blobs which
  //    came before it (which score better).
  // 3) If there is an overlap, delete it from the list.

  SortBlobs(blob_crds,
            blob_diameters, 
            blob_scores,
            sort_blob_criteria,
            false,
            nullptr,
            pReportProgress);

  if (pReportProgress)
    *pReportProgress
      << "  -- Attempting to allocate space for one more image.       --\n"
      << "  -- (If this crashes your computer, find a computer with   --\n"
      << "  --  more RAM and use \"ulimit\", OR use a smaller image.)   --\n";

  // Occupancy table
  //     (originally named "bool ***aaabOcc")
  // This table stores a list of blobs which occupy a given voxel.

  int occupancy_table_size[3];
  int bounds_min[3] = {0, 0, 0};
  int bounds_max[3] = {-1, -1, -1};
  
  for (int i=0; i < blob_crds.size(); i++) {
    for (int d=0; d < 3; d++) {
      Scalar reff = ceil(blob_diameters[i]/2); //blob radii in units of voxels
      if ((blob_crds[i][d] - reff < bounds_min[d]) ||
          (bounds_min[d] > bounds_max[d]))
        bounds_min[d] = blob_crds[i][d] - reff;
      if ((blob_crds[i][d] + reff > bounds_max[d]) ||
          (bounds_min[d] > bounds_max[d]))
        bounds_max[d] = blob_crds[i][d] + reff;
    }
  }

  // Because it is such an enormous table, the size of the table does not equal
  // the size of the original image. Instead we reduce it by a factor of "scale"
  for (int d=0; d < 3; d++)
    occupancy_table_size[d] = (1 + bounds_max[d] - bounds_min[d]) / scale;


  //vector<vector<vector<vector size_t> > > >  
  //  vvvvOcc(occupancy_table_size[3],
  //          vector<
  //          vector<
  //          <vector size_t>(0)
  //          >(occupancy_table_size[1])
  //          >(occupancy_table_size[2])
  //          );
  //
  // The nested for-loop below is less confusing, so I use this instead:

  vector<vector<vector<vector<size_t> > > > vvvvOcc(occupancy_table_size[2]);
  for(int iz = 0; iz < occupancy_table_size[2]; iz++) {
    vvvvOcc[iz] = vector<vector<vector<size_t> > >(occupancy_table_size[1]);
    for(int iy = 0; iy < occupancy_table_size[1]; iy++) {
      vvvvOcc[iz][iy] = vector<vector<size_t> >(occupancy_table_size[0]);
      //for(int ix = 0; ix < occupancy_table_size[0]; ix++) {
      //  vvvvOcc[iz][iy][ix] = vector<size_t>(0);
      //}
    }
  }

  if (pReportProgress)
    *pReportProgress << "  allocating another blob list copy." << endl;

  vector<array<Scalar,3> > blob_crds_cpy;
  vector<Scalar> blob_diameters_cpy;
  vector<Scalar> blob_scores_cpy;


  if (pReportProgress)
    *pReportProgress
      << "  detecting collisions between "<<blob_crds.size()<<" blobs... ";


  // Loop through all of the blobs and fill the occupancy table.  If a given
  // blob overlaps with an existing blob (presumabely with a better score),
  // then we check to see how much they overlap before deciding whether to
  // discard the new blob.

  for (size_t i = 0; i < blob_crds.size(); i++)
  {
    size_t n_blobs = blob_crds.size();
    assert(n_blobs == blob_diameters.size());
    assert(n_blobs == blob_scores.size());
    bool discard = false;
    Scalar reff_ = blob_diameters[i]/2; //blob radii in units of voxels
    Scalar Reff_ = reff_ / scale; //blob radii expressed in "low rez" units
    int Reff = ceil(Reff_) + 1; //round up (and add 1 for uncertainty in center)
    int Reffsq = Reff*Reff;
    Scalar ix = blob_crds[i][0];      //coordinates of the center of the blob
    Scalar iy = blob_crds[i][1];
    Scalar iz = blob_crds[i][2];
    int Ix = floor((ix - bounds_min[0]) / scale); //blob center coords (lowrez version)
    int Iy = floor((iy - bounds_min[1]) / scale);
    int Iz = floor((iz - bounds_min[2]) / scale);
    for(int Jz = -Reff; Jz <= Reff && (! discard); Jz++) {
      for(int Jy = -Reff; Jy <= Reff && (! discard); Jy++) {
        for(int Jx = -Reff; Jx <= Reff && (! discard); Jx++) {
          int rsq = Jx*Jx + Jy*Jy + Jz*Jz;
          if (! ((0 <= Ix+Jx) && (Ix+Jx < occupancy_table_size[0]) &&
                 (0 <= Iy+Jy) && (Iy+Jy < occupancy_table_size[1]) &&
                 (0 <= Iz+Jz) && (Iz+Jz < occupancy_table_size[2])))
            continue;
          if (rsq > Reffsq)
            continue;

          // loop over the other blobs which occupy this location (if any)
          vector<size_t>& overlapping_blobs = vvvvOcc[Iz+Jz][Iy+Jy][Ix+Jx];
          for (size_t _k=0; _k < overlapping_blobs.size(); _k++) {
            size_t k = overlapping_blobs[_k];
            Scalar kx = blob_crds[k][0];
            Scalar ky = blob_crds[k][1];
            Scalar kz = blob_crds[k][2];
            Scalar rik = sqrt((ix-kx)*(ix-kx)+(iy-ky)*(iy-ky)+(iz-kz)*(iz-kz));
            Scalar ri = blob_diameters[i]/2;
            Scalar rk = blob_diameters[k]/2;
            Scalar vol_overlap = CalcSphereOverlap(rik, ri, rk);
            if (rik < (ri + rk) * min_radial_separation_ratio) {
              discard = true;
            }
            Scalar vi = (4*M_PI/3)*(ri*ri*ri);
            Scalar vk = (4*M_PI/3)*(rk*rk*rk);
            Scalar v_large = vi;
            Scalar v_small = vk;
            if (vk > vi) {
              v_large = vk;
              v_small = vi;
            }
            if ((vol_overlap / v_small > max_volume_overlap_small) ||
                (vol_overlap / v_large > max_volume_overlap_large)) {
              discard = true;
              // REMOVE THIS CRUFT:
              //if (pReportProgress)
              //  *pReportProgress << "discarding blob ("
              //                   << ix << "," << iy << "," << iz
              //                   << "): vover=" << vol_overlap
              //                   << ", vsmall=" << v_small
              //                   << ", vlarge=" << v_large << "\n";
            }
          }
        }
      }
    }
    if (discard) {
      // then don't add an entry to the list of blobs
      continue;
    }
    else {
      // then add an entry to the list of blobs (stored in the following arrays)
      blob_crds_cpy.push_back(blob_crds[i]);
      blob_diameters_cpy.push_back(blob_diameters[i]);
      blob_scores_cpy.push_back(blob_scores[i]);
      // mark the pixels within blob i as occupied by blob i
      for(int Jz = -Reff; Jz <= Reff && (! discard); Jz++) {
        for(int Jy = -Reff; Jy <= Reff && (! discard); Jy++) {
          for(int Jx = -Reff; Jx <= Reff && (! discard); Jx++) {
            int rsq = Jx*Jx + Jy*Jy + Jz*Jz;
            if (! ((0 <= Ix+Jx) && (Ix+Jx < occupancy_table_size[0]) &&
                   (0 <= Iy+Jy) && (Iy+Jy < occupancy_table_size[1]) &&
                   (0 <= Iz+Jz) && (Iz+Jz < occupancy_table_size[2])))
              continue;
            if (rsq > Reffsq)
              continue;
            vvvvOcc[Iz+Jz][Iy+Jy][Ix+Jx].push_back(i);
          }
        }
      }
    } //else clause to "if (discard)"

  } //for (size_t i = 0; i < blob_crds.size(); i++)

  blob_crds = blob_crds_cpy;
  blob_diameters = blob_diameters_cpy;
  blob_scores = blob_scores_cpy;

  if (pReportProgress)
    *pReportProgress << "done.\n";
} //DiscardOverlappingBlobs()







/// @brief  Discard blobs whose centers lie outside the 
///         "masked" region defined by aaafMask.

template<typename Scalar>
void
DiscardMaskedBlobs(vector<array<Scalar,3> >& blob_crds, //!< location of each blob
                   // optional arguments:
                   vector<Scalar> &blob_diameters=nullptr,  //!< diameger of each blob
                   vector<Scalar> &blob_scores=nullptr, //!< priority of each blob
                   Scalar const *const *const *aaafMask = nullptr, //!< if not nullptr then discard blobs whose centers at (ix,iy,iz) satisfy aaafMask[iz][iy][ix] == 0.0
                   ostream *pReportProgress=nullptr)  //!< print progress to the user?

{
  if (pReportProgress)
    *pReportProgress << "  allocating another blob list copy." << endl;

  vector<array<Scalar,3> > blob_crds_cpy;
  vector<Scalar> blob_diameters_cpy;
  vector<Scalar> blob_scores_cpy;

  if (pReportProgress && (aaafMask != nullptr))
    *pReportProgress
      << "  checking which blobs lie within the mask (out of "
      << blob_crds.size() << " blobs)." << endl;

  size_t n_discarded = 0;
  for (size_t i = 0; i < blob_crds.size(); i++)
  {
    int ix = floor(blob_crds[i][0] + 0.5);
    int iy = floor(blob_crds[i][1] + 0.5);
    int iz = floor(blob_crds[i][2] + 0.5);
    if ((aaafMask) && (aaafMask[iz][iy][ix] == 0.0)) {
      n_discarded++;
      continue;
    }
    else {
      blob_crds_cpy.push_back(blob_crds[i]);
      blob_diameters_cpy.push_back(blob_diameters[i]);
      blob_scores_cpy.push_back(blob_scores[i]);
    } 
  } //for (size_t i = 0; i < blob_crds.size(); i++)

  if (pReportProgress && (aaafMask != nullptr))
    *pReportProgress
      << "  discarded " << n_discarded
      << "  blobs lying outside the mask." << endl;

  blob_crds = blob_crds_cpy;
  blob_diameters = blob_diameters_cpy;
  blob_scores = blob_scores_cpy;
} // DiscardMaskedBlobs()




/// @brief find either an upper or a lower bound for a list
///        of 1-D features (scores).
///        (Although both upper and lower bounds are calculated, 
///         only one of them should have a non-infinite value.
///         Otherwise your data is not one-sided.)
///
/// @note  In this variant of the function, the positive and negative training
///        data are saved in separate arrays.
///
/// @return This function returns void.
///         Assuming they are not nullptr, the results are stored in:
///           *pthreshold_lower_bound
///           *pthreshold_uppwer_bound

template<typename Scalar>
void

ChooseBlobScoreThresholds(const vector<array<Scalar,3> >& blob_crds, //!< location of each blob (in voxels, sorted by score in increasing priority)
                          const vector<Scalar>& blob_diameters,  //!< diameger of each blob (sorted by score in increasing priority)
                          const vector<Scalar>& blob_scores, //!< priority of each blob (sorted by score in increasing priority)
                          const vector<array<Scalar,3> >& training_set_pos, //!< locations of blob-like things we are looking for
                          const vector<array<Scalar,3> >& training_set_neg, //!< locations of blob-like things we want to ignore
                          Scalar *pthreshold_lower_bound = nullptr, //!< return threshold to the caller
                          Scalar *pthreshold_upper_bound = nullptr, //!< return threshold to the caller
                          SortBlobCriteria sort_blob_criteria = PRIORITIZE_HIGH_MAGNITUDE_SCORES, //!< give priority to high or low scoring blobs?
                          ostream *pReportProgress = nullptr //!< report progress back to the user?
                          )
{
  assert(blob_crds.size() == blob_diameters.size());
  assert(blob_crds.size() == blob_scores.size());

  size_t Nn = training_set_neg.size();
  size_t Np = training_set_pos.size();

  // Figure out which training_data coordinates lie sufficiently close
  // to one of the blobs to be counted.  Ignore the others.
  // While doing this, also figure out the score of each training datum.
  // (This is the score of the blob that is nearby, if applicable.)

  // Concatinate all of the training data together.
  vector<array<Scalar,3> > training_set_crds = training_set_pos;
  training_set_crds.insert(training_set_crds.end(),
                           training_set_neg.begin(),
                           training_set_neg.end());

  // training_set_accepted = true or false depending on whether it is 
  //                      part of the positive (accepted) training set,
  //                      or the negative (rejected) training set
  vector<bool> training_set_accepted =  vector<bool>(Np + Nn, true);
  for (size_t i = Np; i < Np + Nn; i++)
    training_set_accepted[i] = false;                              

  _ChooseBlobScoreThresholds(blob_crds,
                             blob_diameters,
                             blob_scores,
                             training_set_crds,
                             training_set_accepted,
                             pthreshold_lower_bound,
                             pthreshold_upper_bound,
                             sort_blob_criteria,
                             pReportProgress);

} //ChooseBlobScoreThresholds()





/// @brief find either an upper or a lower bound for a list
///        of 1-D features (scores).
///        (Although both upper and lower bounds are calculated, 
///         only one of them should have a non-infinite value.
///         Otherwise your data is not one-sided.)
///
/// @note  This version of function works with multiple
///        independent training sets.  (That is, training sets corresponding
///        to independent sets of overlapping blob_crds, blob_diameters, ...
///        In practice, these different training sets and blob lists are taken
///        from different images.  The size of the following arguments should
///        equal this number of images:   blob_crds, blob_diameters,
///        blob_scores, training_set_crds, training_set_accepted
///
/// @note  In this variant of the function, the positive and negative training
///        data are saved in separate arrays.
///
/// @return This function returns void.
///         Assuming they are not nullptr, the results are stored in:
///           *pthreshold_lower_bound
///           *pthreshold_uppwer_bound

template<typename Scalar>
void
ChooseBlobScoreThresholdsMulti(const vector<vector<array<Scalar,3> > >& blob_crds, //!< location of each blob (in voxels, sorted by score in increasing priority)
                               const vector<vector<Scalar> >& blob_diameters,  //!< diameger of each blob (sorted by score in increasing priority)
                               const vector<vector<Scalar> >& blob_scores, //!< priority of each blob (sorted by score in increasing priority)
                               const vector<vector<array<Scalar,3> > >& training_set_pos, //!< locations of blob-like things 
                               const vector<vector<array<Scalar,3> > >& training_set_neg, //!< locations of blob-like things
                               Scalar *pthreshold_lower_bound = nullptr, //!< return threshold to the caller
                               Scalar *pthreshold_upper_bound = nullptr, //!< return threshold to the caller
                               SortBlobCriteria sort_blob_criteria = PRIORITIZE_HIGH_MAGNITUDE_SCORES, //!< give priority to high or low scoring blobs?
                               ostream *pReportProgress = nullptr //!< report progress back to the user?
                               )
{
  int Nsets = training_set_pos.size();
  assert(Nsets == training_set_neg.size());
  assert(Nsets == blob_crds.size());
  assert(Nsets == blob_diameters.size());
  assert(Nsets == blob_scores.size());

  vector<vector<array<Scalar,3> > > training_set_crds(Nsets);
  vector<vector<bool> > training_set_accepted(Nsets);

  // Loop over all of the different training sets:
  for (int I = 0; I < Nsets; I++) {

    // Concatinate the positive and negative training data together.
    training_set_crds[I] = training_set_pos[I];
    training_set_crds[I].insert(training_set_crds[I].end(),
                                training_set_neg[I].begin(),
                                training_set_neg[I].end());
    size_t Np = training_set_pos[I].size();
    size_t Nn = training_set_neg[I].size();
    // training_set_accepted = true or false depending on whether it is 
    //                      part of the positive (accepted) training set,
    //                      or the negative (rejected) training set
    training_set_accepted[I] =  vector<bool>(Np + Nn, true);
    for (size_t j = Np; j < Np + Nn; j++)
      training_set_accepted[I][j] = false;
  }

  _ChooseBlobScoreThresholdsMulti(blob_crds,
                                  blob_diameters,
                                  blob_scores,
                                  training_set_crds,
                                  training_set_accepted,
                                  pthreshold_lower_bound,
                                  pthreshold_upper_bound,
                                  sort_blob_criteria,
                                  pReportProgress);

} //ChooseBlobScoreThresholdsMulti()




/// @brief  This function discards blobs according to their scores.
///         It does this by comparing the scores of these blobs with 
///         the scores of blobs provided by the caller which were either
///         accepted or rejected (training_set_pos or training_set_neg).
/// @note:  This function was intended to be used when it is possible to use
///         a single score threshold to distinguish good blobs from bad ones.
///         (It was not not intended to be used when the blobs you are looking
///          for have scores which lie within one or more narrow intervals.)

template<typename Scalar>
void
DiscardBlobsByScoreSupervised(vector<array<Scalar,3> >& blob_crds, //!< location of each blob
                              vector<Scalar>& blob_diameters,  //!< diameger of each blob
                              vector<Scalar>& blob_scores, //!< priority of each blob
                              const vector<array<Scalar,3> >& training_set_pos, //!< locations of blob-like things we are looking for
                              const vector<array<Scalar,3> >& training_set_neg, //!< locations of blob-like things we want to ignore
                              SortBlobCriteria sort_blob_criteria = PRIORITIZE_HIGH_MAGNITUDE_SCORES, //!< give priority to high or low scoring blobs?
                              Scalar *pthreshold_lower_bound = nullptr, //!< optional: return threshold to the caller
                              Scalar *pthreshold_upper_bound = nullptr, //!< optional: return threshold to the caller
                              ostream *pReportProgress = nullptr //!< report progress back to the user?
                              )
{
  assert(blob_crds.size() == blob_diameters.size());
  assert(blob_crds.size() == blob_scores.size());

  Scalar threshold_lower_bound;
  Scalar threshold_upper_bound;

  // Find the range of scores [threshold_lower_bound, threshold_upper_bound]
  // which minimizes the number of incorrectly classified blobs.

  ChooseBlobScoreThresholds(blob_crds,
                            blob_diameters,
                            blob_scores,
                            training_set_pos,
                            training_set_neg,
                            &threshold_lower_bound, //<-store threshold here
                            &threshold_upper_bound, //<-store threshold here
                            sort_blob_criteria,
                            pReportProgress);

  if (pthreshold_lower_bound)
    *pthreshold_lower_bound = threshold_lower_bound;
  if (pthreshold_upper_bound)
    *pthreshold_upper_bound = threshold_upper_bound;
  
  if (pReportProgress)
    *pReportProgress << "  allocating another blob list copy." << endl;

  vector<array<Scalar,3> > blob_crds_cpy;
  vector<Scalar> blob_diameters_cpy;
  vector<Scalar> blob_scores_cpy;

  for (int i = 0; i < blob_crds.size(); i++) {
    if ((blob_scores[i] >= threshold_lower_bound) &&
        (blob_scores[i] <= threshold_upper_bound))
    {
      blob_crds_cpy.push_back(blob_crds[i]);
      blob_diameters_cpy.push_back(blob_diameters[i]);
      blob_scores_cpy.push_back(blob_scores[i]);
    }  
  }

  blob_crds = blob_crds_cpy;
  blob_diameters = blob_diameters_cpy;
  blob_scores = blob_scores_cpy;

} //DiscardBlobsByScoreSupervised()






/// @brief  Calculate matrix of 2nd derivatives (the hessian)
///         as well as the the vector of 1st derivatives (the gradient)
///         of the source image (aaafSource), at every location where aaafMask
///         is non-zero (or everywhere if aaafMask is nullptr)
///         Apply a Gaussian blur to the image (of width sigma) beforehand,
///         (truncating the blur filter at a distance of truncate_ratio*sigma
///          voxels from the center of the Gaussian).
///         Hessians and Gradients are saved in 3-dimensional arrays of
///         either "TensorContainers", or "VectorContainers".
///         Both of these objects must support array subscripting.
///         Each "TensorContainer" object is expected to behave like
///         a one-dimensional array of 6 scalars (hence the "aaaafHessian"
///         argument behaves like a 4-dimensional array)
///         Each "VectorContainer" object is expected to behave like
///         a one-dimensional array of 3 scalars.

template<typename Scalar, typename VectorContainer=Scalar*, typename TensorContainer=Scalar*>

void
CalcHessian(int const image_size[3], //!< source image size
            Scalar const *const *const *aaafSource, //!< source image
            VectorContainer ***aaaafGradient,  //!< save results here (if not nullptr)
            TensorContainer ***aaaafHessian, //!< save results here (if not nullptr)
            Scalar const *const *const *aaafMask,  //!< ignore voxels where mask==0
            Scalar sigma,  //!< Gaussian width in x,y,z drections
            Scalar truncate_ratio=2.5,  //!< how many sigma before truncating?
            ostream *pReportProgress = nullptr  //!< print progress to the user?
            )
{
  assert(aaafSource);
  assert(aaaafHessian);
  
  int truncate_halfwidth = floor(sigma * truncate_ratio);

  // Here we use the fast, sloppy way to compute gradients and Hessians:
  // First smooth the image,
  // Then infer derivatives from finite differences.
  //
  // This only works for large sigma (sigma at least 1.0).
  //
  // A slower, more rigorous approach would be to convolve the image
  // with the derivatives of a Gaussian.  (If we need to analyze images
  // with very thin, closely spaced objects, we can try that way later.)

  // First, apply the Gaussian filter to the image
  if (pReportProgress)
    *pReportProgress
      << " -- Attempting to allocate space for one more image.       --\n"
      << " -- (If this crashes your computer, find a computer with   --\n"
      << " --  more RAM and use \"ulimit\", OR use a smaller image.)   --\n";

  Scalar ***aaafSmoothed;
  Scalar *afSmoothed;
  Alloc3D(image_size,
          &afSmoothed,
          &aaafSmoothed);

  if (pReportProgress)
    *pReportProgress << "\n";

  ApplyGauss(image_size,
             aaafSource,
             aaafSmoothed,
             aaafMask,
             sigma,
             truncate_halfwidth,
             false,
             pReportProgress);

  if (pReportProgress)
    *pReportProgress << "Calculating the Hessian associated with each voxel\n";

  if ((image_size[0] < 3) ||
      (image_size[1] < 3) ||
      (image_size[2] < 3))
    throw VisfdErr("Error: CalcHessian() requires an image that is at least 3 voxels\n"
                   "       wide in the x,y,z directions.\n");
  assert(image_size[0] >= 3);
  assert(image_size[1] >= 3);
  assert(image_size[2] >= 3);

  // Now compute gradients and hessians

  for (int iz = 0; iz < image_size[2]; iz++) {
    #pragma omp parallel for collapse(2)
    for (int iy = 0; iy < image_size[1]; iy++) {
      for (int ix = 0; ix < image_size[0]; ix++) {
        if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
          continue;

        if (aaaafGradient) {
          Scalar gradient[3];

          CalcGradientFiniteDifferences(aaafSmoothed,
                                        ix, iy, iz,
                                        gradient,
                                        image_size);

          // Optional: Insure that the resulting gradient is dimensionless:
          // (Lindeberg 1993 "On Scale Selection for Differential Operators")

          gradient[0] *= sigma;
          gradient[1] *= sigma;
          gradient[2] *= sigma;
          
          aaaafGradient[iz][iy][ix][0] = gradient[0];
          aaaafGradient[iz][iy][ix][1] = gradient[1];
          aaaafGradient[iz][iy][ix][2] = gradient[2];

          #ifndef NDEBUG
          if ((ix==image_size[0]/2) &&
              (iy==image_size[1]/2) &&
              (iz==image_size[2]/2))
          {
            if (pReportProgress)
              *pReportProgress
                 << "[iz][iy][ix]=["<<iz<<"]["<<iy<<"]["<<ix<<"], "
                 << "gradient[0] = " << aaaafGradient[iz][iy][ix][0] << endl;
          }
          #endif  //#ifndef NDEBUG

        }

        if (aaaafHessian) {

          Scalar hessian[3][3];
          CalcHessianFiniteDifferences(aaafSmoothed,
                                       ix, iy, iz,
                                       hessian,
                                       image_size);

          #ifndef NDEBUG
          // DEBUG: REMOVE THE NEXT IF STATMENT AFTER DEBUGGING IS FINISHED
          if ((ix==image_size[0]/2) && //if ((ix==78) && 
              (iy==image_size[1]/2) &&
              (iz==image_size[2]/2))
          {
            if (pReportProgress)
              *pReportProgress
                << "[iz][iy][ix]=["<<iz<<"]["<<iy<<"]["<<ix<<"], "
                << "hessian[0][0] = " << hessian[0][0] << endl;
          }
          #endif  //#ifndef NDEBUG


          // Optional: Insure that the result is dimensionless:
          // (Lindeberg 1993 "On Scale Selection for Differential Operators")
          for (int di=0; di < 3; di++)
            for (int dj=0; dj < 3; dj++)
              hessian[di][dj] *= sigma*sigma;

          // To reduce memory consumption,
          // save the resulting 3x3 matrix in a smaller 1-D array whose index
          // is given by MapIndices_3x3_to_linear[][]
          for (int di = 0; di < 3; di++)
            for (int dj = di; dj < 3; dj++)
              aaaafHessian[iz][iy][ix][ MapIndices_3x3_to_linear[di][dj] ]
                = hessian[di][dj];

        } //if (aaaafHessian)
      } //for (int ix = 1; ix < image_size[0]-1; ix++)
    } //for (int iy = 1; iy < image_size[1]-1; iy++)
  } //for (int iz = 1; iz < image_size[2]-1; iz++)

  Dealloc3D(image_size,
            &afSmoothed,
            &aaafSmoothed);

} //CalcHessian()






/// CalcMomentTensor()
/// This may be algebraically equivalent to CalcHessian()
/// However this version of the function might be more robust for small ridges.
/// (This is because I apply the derivative to the Gaussian filter
///  before applying the filter, ...instead of applying the Gaussian filter
///  first and then taking finite differences afterwards.  If the width of the 
///  object being detected is not much more than 3-voxel wide, the 3-voxel wide
///  differences used in the other implementation are a large source of error.)
/// Unfortunately this version is slower and needs much more memory.
/// Eventually, I might elliminate one of these implementations.

template<typename Scalar, typename FirstMomentContainer, typename SecondMomentContainer>

void
CalcMomentTensor(int const image_size[3], //!< source image size
                 Scalar const *const *const *aaafSource, //!< source image
                 FirstMomentContainer ***aaaaf1stMoment,  //!< save results here (if not nullptr)
                 SecondMomentContainer ***aaaaf2ndMoment, //!< save results here (if not nullptr)
                 Scalar const *const *const *aaafMask,  //!< ignore voxels where mask==0
                 Scalar sigma,  //!< Gaussian width in x,y,z drections
                 Scalar truncate_ratio=2.5,  //!< how many sigma before truncating?
                 ostream *pReportProgress = nullptr  //!< print progress to the user?
                 )
{
  assert(image_size);
  assert(aaafSource);
  
  int truncate_halfwidth = floor(sigma * truncate_ratio);

  //calculate the filter used for the 0'th moment (ordinary Gaussian filter)
  Filter1D<Scalar, int> filter0 = GenFilterGauss1D(sigma,
                                                   truncate_halfwidth);

  //calculate the filter used for the 1st moment (Guassian(x) * x)
  Filter1D<Scalar, int> filter1 = GenFilterGauss1D(sigma,
                                                   truncate_halfwidth);
  for (int i = -truncate_halfwidth; i <= truncate_halfwidth; i++)
    filter1.afW[i] *= i;


  //calculate the filter used for the 2nd moment (Guassian(x) * x^2)
  Filter1D<Scalar, int> filter2 = GenFilterGauss1D(sigma,
                                                   truncate_halfwidth);
  for (int i = -truncate_halfwidth; i <= truncate_halfwidth; i++)
    filter2.afW[i] *= i*i;


  Scalar ***aaafNorm;
  Scalar *afNorm;
  Alloc3D(image_size,
          &afNorm,
          &aaafNorm);

  if (aaafMask) {
    // Calculate the normalization we need by blurring the mask by the Gaussian
    ApplyGauss(image_size,
               aaafMask,   // <-- use the mask as the source image
               aaafNorm,   // <-- save result here
               aaafMask,
               sigma, // width of Gaussian
               truncate_halfwidth,
               true,
               pReportProgress);
  }


  Filter1D<Scalar, int> aFilter[3];


  if (aaaaf1stMoment) {

    if (pReportProgress)
      *pReportProgress
        << " -- Attempting to allocate space for 3 more images.        --\n"
        << " -- (If this crashes your computer, find a computer with   --\n"
        << " --  more RAM and use \"ulimit\", OR use a smaller image.) --\n";

    Scalar ***aaafIx;
    Scalar *afIx;
    Alloc3D(image_size,
            &afIx,
            &aaafIx);

    Scalar ***aaafIy;
    Scalar *afIy;
    Alloc3D(image_size,
            &afIy,
            &aaafIy);

    Scalar ***aaafIz;
    Scalar *afIz;
    Alloc3D(image_size,
            &afIz,
            &aaafIz);

    // calculate the x moment (=x^1 * y^0 * z^0)
    aFilter[0] = filter1;  // x^1
    aFilter[1] = filter0;  // y^0
    aFilter[2] = filter0;  // z^0
    return ApplySeparable(image_size, 
                          aaafSource,
                          aaafIx,
                          aaafMask,
                          aFilter,
                          (aaafMask == nullptr), //don't normalize if theres a mask
                          pReportProgress);

    // calculate the y moment (=x^0 * y^1 * z^0)
    aFilter[0] = filter0;  // x^0
    aFilter[1] = filter1;  // y^1
    aFilter[2] = filter0;  // z^0
    return ApplySeparable(image_size, 
                          aaafSource,
                          aaafIy,
                          aaafMask,
                          aFilter,
                          (aaafMask == nullptr), //don't normalize if theres a mask
                          pReportProgress);

    // calculate the x moment (=x^0 * y^0 * z^1)
    aFilter[0] = filter0;  // x^0
    aFilter[1] = filter0;  // y^0
    aFilter[2] = filter1;  // z^1
    return ApplySeparable(image_size, 
                          aaafSource,
                          aaafIz,
                          aaafMask,
                          aFilter,
                          (aaafMask == nullptr), //don't normalize if theres a mask
                          pReportProgress);

    if (aaafMask) {
      // If a mask was specified, then divide by the contribution from the mask
      // (An infinite sized mask should result in a contribution/weight of 1)
      for(int iz=0; iz<image_size[2]; iz++) {
        for(int iy=0; iy<image_size[1]; iy++) {
          for(int ix=0; ix<image_size[0]; ix++) {
            if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
              continue;
            assert(aaafNorm[iz][iy][ix] > 0.0);
            aaafIx[iz][iy][ix] /= aaafNorm[iz][iy][ix];
            aaafIy[iz][iy][ix] /= aaafNorm[iz][iy][ix];
            aaafIz[iz][iy][ix] /= aaafNorm[iz][iy][ix];
          }
        }
      }
    }

    for (int iz = 1; iz < image_size[2]-1; iz++) {
      for (int iy = 1; iy < image_size[1]-1; iy++) {
        for (int ix = 1; ix < image_size[0]-1; ix++) {
          if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
            continue;
          Scalar first_deriv[3];
          first_deriv[0] = aaafIx[iz][iy][ix];
          first_deriv[1] = aaafIy[iz][iy][ix];
          first_deriv[2] = aaafIz[iz][iy][ix];

          // Optional: Insure that the resulting first_deriv is dimensionless:
          // (Lindeberg 1993 "On Scale Selection for Differential Operators")
          first_deriv[0] *= sigma;
          first_deriv[1] *= sigma;
          first_deriv[2] *= sigma;
          
          // Store the result in aaaaf1stMoment[]
          aaaaf1stMoment[iz][iy][ix][0] = first_deriv[0];
          aaaaf1stMoment[iz][iy][ix][1] = first_deriv[1];
          aaaaf1stMoment[iz][iy][ix][2] = first_deriv[2];
        }
      }
    }
    Dealloc3D(image_size,
              &afIx,
              &aaafIx);
    Dealloc3D(image_size,
              &afIy,
              &aaafIy);
    Dealloc3D(image_size,
              &afIz,
              &aaafIz);
  } //if (pFirstMoment)



  if (aaaaf2ndMoment) {
    if (pReportProgress)
      *pReportProgress << "\n"
        " ------ Calculating the average of nearby voxels: ------\n";
    // P = original image (after subtracting average nearby intensities):

    Scalar ***aaafP;
    Scalar *afP;
    Alloc3D(image_size,
            &afP,
            &aaafP);

    ApplyGauss(image_size,
               aaafSource,
               aaafP,   // <-- save result here
               aaafMask,
               sigma, // width of Gaussian
               truncate_halfwidth,
               true,
               pReportProgress);

    // Subtract the average value from the image intensity, and store in P:
    for(int iz=0; iz<image_size[2]; iz++)
      for(int iy=0; iy<image_size[1]; iy++)
        for(int ix=0; ix<image_size[0]; ix++)
          aaafP[iz][iy][ix] = aaafSource[iz][iy][ix] - aaafP[iz][iy][ix];

    if (pReportProgress)
      *pReportProgress
        << " -- Attempting to allocate space for 6 more images.        --\n"
        << " -- (If this crashes your computer, find a computer with   --\n"
        << " --  more RAM and use \"ulimit\", OR use a smaller image.) --\n";

    Scalar ***aaafIxx;
    Scalar *afIxx;
    Alloc3D(image_size,
            &afIxx,
            &aaafIxx);

    Scalar ***aaafIyy;
    Scalar *afIyy;
    Alloc3D(image_size,
            &afIyy,
            &aaafIyy);

    Scalar ***aaafIzz;
    Scalar *afIzz;
    Alloc3D(image_size,
            &afIzz,
            &aaafIzz);

    Scalar ***aaafIxy;
    Scalar *afIxy;
    Alloc3D(image_size,
            &afIxy,
            &aaafIxy);

    Scalar ***aaafIyz;
    Scalar *afIyz;
    Alloc3D(image_size,
            &afIyz,
            &aaafIyz);

    Scalar ***aaafIxz;
    Scalar *afIxz;
    Alloc3D(image_size,
            &afIxz,
            &aaafIxz);

    // calculate the x*x moment of the innertia (=x^2 * y^0 * z^0)
    aFilter[0] = filter2;  // x^2
    aFilter[1] = filter0;  // y^0
    aFilter[2] = filter0;  // z^0
    return ApplySeparable(image_size, 
                          aaafP,
                          aaafIxx,
                          aaafMask,
                          aFilter,
                          (aaafMask == nullptr), //don't normalize if theres a mask
                          pReportProgress);

    // calculate the y*y moment of the innertia (=x^0 * y^2 * z^0)
    aFilter[0] = filter0;  // x^0
    aFilter[1] = filter2;  // y^2
    aFilter[2] = filter0;  // z^0
    return ApplySeparable(image_size, 
                          aaafP,
                          aaafIyy,
                          aaafMask,
                          aFilter,
                          (aaafMask == nullptr), //don't normalize if theres a mask
                          pReportProgress);

    // calculate the z*z moment of the innertia (=x^0 * y^0 * z^2)
    aFilter[0] = filter0;  // x^0
    aFilter[1] = filter0;  // y^0
    aFilter[2] = filter2;  // z^2
    return ApplySeparable(image_size, 
                          aaafP,
                          aaafIzz,
                          aaafMask,
                          aFilter,
                          (aaafMask==nullptr), //don't normalize if theres a mask
                          pReportProgress);

    // calculate the x*y moment of the innertia (=x^1 * y^1 * z^0)
    aFilter[0] = filter1;  // x^1
    aFilter[1] = filter1;  // y^1
    aFilter[2] = filter0;  // z^0
    return ApplySeparable(image_size, 
                          aaafP,
                          aaafIxy,
                          aaafMask,
                          aFilter,
                          (aaafMask==nullptr), //don't normalize if theres a mask
                          pReportProgress);

    // calculate the y*z moment of the innertia (=x^0 * y^1 * z^1)
    aFilter[0] = filter0;  // x^0
    aFilter[1] = filter1;  // y^1
    aFilter[2] = filter1;  // z^1
    return ApplySeparable(image_size, 
                          aaafP,
                          aaafIyz,
                          aaafMask,
                          aFilter,
                          (aaafMask==nullptr), //don't normalize if theres a mask
                          pReportProgress);

    // calculate the x*z moment of the innertia (=x^1 * y^0 * z^1)
    aFilter[0] = filter1;  // x^1
    aFilter[1] = filter0;  // y^0
    aFilter[2] = filter1;  // z^1
    return ApplySeparable(image_size, 
                          aaafP,
                          aaafIxz,
                          aaafMask,
                          aFilter,
                          (aaafMask==nullptr), //don't normalize if theres a mask
                          pReportProgress);

    if (aaafMask) {
      // If a mask was specified, then divide by the contribution from the mask
      // (An infinite sized mask should result in a contribution/weight of 1)
      for(int iz=0; iz<image_size[2]; iz++) {
        for(int iy=0; iy<image_size[1]; iy++) {
          for(int ix=0; ix<image_size[0]; ix++) {
            if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
              continue;
            assert(aaafNorm[iz][iy][ix] > 0.0);
            aaafIxx[iz][iy][ix] /= aaafNorm[iz][iy][ix];
            aaafIyy[iz][iy][ix] /= aaafNorm[iz][iy][ix];
            aaafIzz[iz][iy][ix] /= aaafNorm[iz][iy][ix];
            aaafIxy[iz][iy][ix] /= aaafNorm[iz][iy][ix];
            aaafIyz[iz][iy][ix] /= aaafNorm[iz][iy][ix];
            aaafIxz[iz][iy][ix] /= aaafNorm[iz][iy][ix];
          }
        }
      }
    }

    for(int iz=0; iz<image_size[2]; iz++) {
      for(int iy=0; iy<image_size[1]; iy++) {
        for(int ix=0; ix<image_size[0]; ix++) {
          if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
            continue;
          // Store the result in aaaaf2ndMoment[].  As usual, to save memory
          // we use the "MapIndices_3x3_to_linear[][]" function to store the
          // entries of the symmetric 3x3 matrix in a 1D array with only 6 entries.
          // (Symmetric 3x3 matrices can have at most 6 unique entries.)
          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[0][0] ] =
            aaafIxx[ix][iy][iz];

          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[0][1] ] =
            aaafIxy[ix][iy][iz];
    
          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[0][2] ] =
            aaafIxz[ix][iy][iz];

          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[1][0] ] =
            aaafIxy[ix][iy][iz];  // (redundant but harmless)

          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[1][1] ] =
            aaafIyy[ix][iy][iz];

          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[1][2] ] =
            aaafIyz[ix][iy][iz];

          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[2][0] ] =
            aaafIxz[ix][iy][iz];  // (redundant but harmless)

          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[2][1] ] =
            aaafIyz[ix][iy][iz];  // (redundant but harmless)

          aaaaf2ndMoment[iz][iy][ix][ MapIndices_3x3_to_linear[2][2] ] =
            aaafIzz[ix][iy][iz];
        }
      }
    }

    if (pReportProgress)
      *pReportProgress << "done ----" << endl;

    Dealloc3D(image_size,
              &afP,
              &aaafP);

    Dealloc3D(image_size,
              &afIxx,
              &aaafIxx);
    Dealloc3D(image_size,
              &afIyy,
              &aaafIyy);
    Dealloc3D(image_size,
              &afIzz,
              &aaafIzz);
    Dealloc3D(image_size,
              &afIxy,
              &aaafIxy);
    Dealloc3D(image_size,
              &afIxz,
              &aaafIxz);
    Dealloc3D(image_size,
              &afIyz,
              &aaafIyz);
  } //if (pSecondMoment)

  Dealloc3D(image_size,
            &afNorm,
            &aaafNorm);

} // CalcMomentTensor()







/// @brief  Convert a volumetric 3D 6-channel image, where each voxel in the
///         image has the (non-redundant) components of a symmetrix 3x3 matrix.
///         The output of this function is another 3D 6-channel image, however
///         each voxel in this image contains the 3-eigenvalues as well as the
///         eigevectors (stored as 3 Shoemake coordinates).
///         If a non-null "aaafMask" argument was specified, voxels in the
///         image are ignored when aaafMask[iz][iy][ix] == 0.
/// @note   The "TensorContainer" object type is expected to behave like
///         a one-dimensional array of 6 scalars.

template<typename Scalar, typename TensorContainer>

void
DiagonalizeHessianImage(int const image_size[3], //!< source image size
                        TensorContainer const *const *const *aaaafSource, //!< input tensor
                        TensorContainer ***aaaafDest, //!< output tensors stored here (can be the same as aaaafSource)
                        Scalar const *const *const *aaafMask,  //!< ignore voxels where mask==0
                        EigenOrderType eival_order = selfadjoint_eigen3::INCREASING_EIVALS, //!< Order of the eigenvalues/eivenvectors.  The default value is typically useful if you are seeking bright objects on a dark background.  Use DECREASING_EIVALS when seeking dark objects on a bright bacground.
                        ostream *pReportProgress = nullptr  //!< print progress to the user?
                        )
{
  assert(aaaafSource);
  if (pReportProgress && aaaafSource)
    *pReportProgress << "\n"
      "---- Diagonalizing the Hessians everywhere (within the mask) ----\n";

  for (int iz = 1; iz < image_size[2]-1; iz++) {
    if (pReportProgress)
      *pReportProgress << "  z="<<iz+1<<" (of "<<image_size[2]<<")" << endl;

    #pragma omp parallel for collapse(2)
    for (int iy = 1; iy < image_size[1]-1; iy++) {
      for (int ix = 1; ix < image_size[0]-1; ix++) {
        if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
          continue;


        #ifndef NDEBUG
        // REMOVE THE NEXT IF STATEMENT AFTER YOU ARE THROUGH DEBUGGING:
        if ((ix==image_size[0]/2) &&
            (iy==image_size[1]/2) &&
            (iz==image_size[2]/2))
        {
          Scalar hessian[3][3];
          for (int di=0; di<3; di++)
            for (int dj=0; dj<3; dj++)
              hessian[di][dj] = aaaafSource[iz][iy][ix][ MapIndices_3x3_to_linear[di][dj] ];
          Scalar quat[4];

          if (pReportProgress)
            *pReportProgress
              << "[iz][iy][ix]=["<<iz<<"]["<<iy<<"]["<<ix<<"]\n"
              << "hessian = \n"
              << "    "<<hessian[0][0]<<","<<hessian[0][1]<<","<<hessian[0][2]<<"\n"
              << "    "<<hessian[1][0]<<","<<hessian[1][1]<<","<<hessian[1][2]<<"\n"
              << "    "<<hessian[2][0]<<","<<hessian[2][1]<<","<<hessian[2][2]<<"\n";

          Scalar eivals[3];
          Scalar eivects[3][3];
          DiagonalizeSym3(hessian,
                          eivals,
                          eivects,
                          eival_order);
          if (Determinant3(eivects) < 0.0) {
            for (int d=0; d<3; d++)
              eivects[0][d] *= -1.0;
          }
          if (pReportProgress)
            *pReportProgress
              << "eivects = \n"
              << "    "<<eivects[0][0]<<","<<eivects[0][1]<<","<<eivects[0][2]<<"\n"
              << "    "<<eivects[1][0]<<","<<eivects[1][1]<<","<<eivects[1][2]<<"\n"
              << "    "<<eivects[2][0]<<","<<eivects[2][1]<<","<<eivects[2][2]<<"\n";

          // Note: Each eigenvector is a currently row-vector in eivects[3][3];
          // It's optional, but I prefer to transpose this, because I think of
          // each eigenvector as a column vector.  Either way should work.
          Transpose3(eivects);
          Matrix2Quaternion(eivects, quat); //convert to 3x3 matrix
          if (pReportProgress)
            *pReportProgress
              << "quat = " << quat[0]<<","<<quat[1]<<","<<quat[2]<<","<<quat[2]<<"\n";
        }
        #endif  //#ifndef NDEBUG




        DiagonalizeFlatSym3(aaaafSource[iz][iy][ix],
                            aaaafDest[iz][iy][ix],
                            eival_order);


        #ifndef NDEBUG
        // REMOVE THE NEXT IF STATEMENT AFTER YOU ARE THROUGH DEBUGGING:
        if ((ix==image_size[0]/2) &&
            (iy==image_size[1]/2) &&
            (iz==image_size[2]/2))
        {
          Scalar shoemake[3]; // <- the eigenvectors stored in "Shoemake" format
          shoemake[0]       = aaaafDest[iz][iy][ix][3];
          shoemake[1]       = aaaafDest[iz][iy][ix][4];
          shoemake[2]       = aaaafDest[iz][iy][ix][5];
          Scalar quat[4];
          Shoemake2Quaternion(shoemake, quat); //convert to quaternion
          Scalar eivects[3][3];
          Quaternion2Matrix(quat, eivects); //convert to 3x3 matrix
          if (pReportProgress)
            *pReportProgress
              << "[iz][iy][ix]=["<<iz<<"]["<<iy<<"]["<<ix<<"]\n"
              << "quat2 = "
              << quat[0]<<","<<quat[1]<<","<<quat[2]<<","<<quat[2]<<"\n"
              << "eivects = \n"
              <<" "<<eivects[0][0]<<","<<eivects[0][1]<<","<<eivects[0][2]<<"\n"
              <<" "<<eivects[1][0]<<","<<eivects[1][1]<<","<<eivects[1][2]<<"\n"
              <<" "<<eivects[2][0]<<","<<eivects[2][1]<<","<<eivects[2][2]<<"\n"
              << endl;
        }
        #endif  //#ifndef NDEBUG

      } //for (int ix = 1; ix < image_size[0]-1; ix++) {
    } //for (int iy = 1; iy < image_size[1]-1; iy++) {
  } //for (int iz = 1; iz < image_size[2]-1; iz++) {
} //DiagonalizeHessianImage()






/// @brief  Read a volumetric 3D 6-channel image, where each voxel in the
///         image has the (non-redundant) components of a symmetrix 3x3 matrix
///         which has been diagonalized using DiagonalizeHessianImage().
///         This function loops over all the voxels and performs the
///         inverse operation.
/// @note   The "TensorContainer" object type is expected to behave like
///         a one-dimensional array of 6 scalars.

template<typename Scalar, typename TensorContainer>

void
UndiagonalizeHessianImage(int const image_size[3],  //!< source image size
                          TensorContainer const *const *const *aaaafSource, //!< input tensor
                          TensorContainer ***aaaafDest, //!< output tensors stored here (can be the same as aaaafSource)
                          Scalar const *const *const *aaafMask,  //!< ignore voxels where mask==0
                          ostream *pReportProgress = nullptr  //!< print progress to the user?
                          )
{
  assert(aaaafSource);

  if (pReportProgress)
    *pReportProgress << "\n"
      "---- Undiagonalizing the Hessians everywhere (within the mask)... "
                     << flush;

  #pragma omp parallel for collapse(2)
  for (int iz = 1; iz < image_size[2]-1; iz++) {
    for (int iy = 1; iy < image_size[1]-1; iy++) {
      for (int ix = 1; ix < image_size[0]-1; ix++) {
        if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
          continue;

        UndiagonalizeFlatSym3(aaaafSource[iz][iy][ix],
                              aaaafDest[iz][iy][ix]);

      } //for (int ix = 1; ix < image_size[0]-1; ix++) {
    } //for (int iy = 1; iy < image_size[1]-1; iy++) {
  } //for (int iz = 1; iz < image_size[2]-1; iz++) {
} //UndiagonalizeHessianImage()





/// @brief: Calculate how "plane"-like a feature along a ridge is
///         from the hessian (the matrix of 2nd derivatives).
///         This function assumes the hessian matrix has been diagonalized
///         and that the first 3 entries of the "diagonalizedHessian"
///         argument are the eigenvalues of the original hessian matrix.

template<typename TensorContainer, typename VectorContainer>

double
ScoreHessianPlanar(TensorContainer diagonalizedHessian,
                   VectorContainer gradient=nullptr)
{
  //typedef decltype(diagonalizedHessian[0]) Scalar;
  double lambda1 = diagonalizedHessian[0];
  double lambda2 = diagonalizedHessian[1];
  double lambda3 = diagonalizedHessian[2];

  // REMOVE THIS CRUFT
  // The "score_ratio" variable is the score function used in Eq(5) of
  // Martinez-Sanchez++Fernandez_JStructBiol2011.
  //double score_ratio;
  //score_ratio = ((abs(lambda1) - sqrt(abs(lambda2*lambda3)))
  //               / SQR(gradient);
  //score_ratio *= score_ratio;

  // REMOVE THIS CRUFT:
  // The following "linear" metric produces interesting results, but
  // the resulting membrane structures that are detected are not well
  // separated from the huge amount of background noise.
  //Scalar Linear_norm = lambda1 - lambda2;
  //score = Linear_norm / SQR(gradient);


  // I decided to try the "Ngamma_norm" metric proposed on p.26 of
  // Lindeberg Int.J.ComputVis.1998,
  // "Edge and ridge detection with automatic scale selection"
  double Nnorm = lambda1*lambda1 - lambda2*lambda2;
  Nnorm *= Nnorm;

  return Nnorm;
} // ScoreHessianPlanar()





/// @brief: Calculate how "plane"-like a feature along a ridge is
///         from the tensor created by the process of tensor voting.
///         This function assumes that "diagonalizedMatrix3x3" has been 
///         diagonalized, and that its first 3 entries 
///         are the eigenvalues of the original hessian matrix.

template<typename TensorContainer>

double
ScoreTensorPlanar(const TensorContainer diagonalizedMatrix3)
{
  //return ScoreHessianPlanar(diagonalizedMatrix3x3,
  //                          nullptr,
  //                          multiplier);
  double lambda1 = diagonalizedMatrix3[0];
  double lambda2 = diagonalizedMatrix3[1];
  return lambda1 - lambda2;  // the "stickness" (See TensorVoting paper)
}



/// @class  TV3D
/// @brief  A class for performing simple tensor-voting image processing 
///         operations in 3D.  Currently only "stick" voting is supported.
///         (Other kinds of tensor voting, such as "plate" and "ball" are not.)
///         This can perform tensor voting for both types of
///         "stick" fields in 3D:
///         (1) Stick-fields corresponding to 2D surface-like features,
///         (2) Stick-fields corresponding to 1D curve-like features.

template<typename Scalar, typename Integer, typename VectorContainer, typename TensorContainer>

class TV3D {

private:

  Scalar sigma;
  Integer exponent;
  Integer halfwidth[3];
  Integer array_size[3];
  Filter3D<Scalar, Integer> radial_decay_lookup;
  array<Scalar, 3> ***aaaafDisplacement;
  array<Scalar, 3> *aafDisplacement;

public:

  TV3D():radial_decay_lookup() {
    Init();
  }

  TV3D(Scalar set_sigma,
       Integer set_exponent,
       Scalar filter_cutoff_ratio=2.5):radial_decay_lookup() {
    Init();
    SetExponent(set_exponent);
    SetSigma(set_sigma, filter_cutoff_ratio);
  }

  TV3D(const Filter3D<Scalar, Integer>& source) {
    Resize(source.halfwidth); // allocates and initializes afH and aaafH
    std::copy(source.aafDisplacement,
              source.aafDisplacement + (array_size[0] * array_size[1] * array_size[2]),
              aafDisplacement);
  }

  ~TV3D() {
    DeallocDisplacement();
  }

  void SetExponent(Scalar set_exponent) {
    exponent = set_exponent;
  }

  void SetSigma(Scalar set_sigma, Scalar filter_cutoff_ratio=2.5)
  {
    sigma = set_sigma;
    Integer halfwidth_single = floor(sigma * filter_cutoff_ratio);
    for (int d=0; d<3; d++)
      halfwidth[d] = halfwidth_single;
    Resize(halfwidth);
  }


  /// @brief  Perform dense stick-voting, using every voxel in the image
  ///         (with non-zero corresponding entries in the aaafMaskSource array)
  ///         as a source, and collecting votes at every voxel
  ///         (with non-zero correspondin entries in the aaafMaskDest array)
  ///         in the aaaafDest array.
  ///         Other kinds of tensor voting ("plate" and "ball")
  ///         are not supported by this function.
  ///         This function can perform tensor voting for both types of
  ///         "stick" fields in 3D:
  ///         (1) Stick-fields corresponding to planar-surface-like features,
  ///         (2) Stick-fields corresponding to curve-like features.
  ///         This function expects a 3D array of "vectors" (aaaafV argument)
  ///         (one vector for each voxel in the original image, 3 numbers each).
  ///         These "vectors" can have user-defined type (implementation),
  ///         however they must support 1-dimensional subscripting (i=0,1,2).
  ///         Optionally, the caller can supply an array of numbers
  ///         ("saliencies") which store the "strength" of each vector.
  ///         If the aaafSaliency[][][] array argument == nullptr, then these
  ///         saliencies will be inferred from the magnitude of the vectors.
  ///         which are stored in the aaaafV array.
  ///         (Otherwise, the vectors are assumed to have been normalized.)
  ///         
  /// After this function is invoked, aaaafDest will store an array of
  /// tensors, one tensor for each voxel in the original image
  /// (unless aaafMaskDest!=nullptr and the corresponding entry there is 0).
  ///
  /// @note:  The computation time for this algorithm is proportional to the 
  ///         number of voxels with non-zero aaafSaliency[][][] values.
  ///         Hence, the speed can be dramatically increased by zeroing
  ///         voxels with low saliency.  For typical cryo-EM images of
  ///         cells, 95% of the voxels can usually be discarded (by setting
  ///         their saliencies to zero), with no effect on the output.

  void
  TVDenseStick(Integer const image_size[3],  //!< source image size
               Scalar const *const *const *aaafSaliency,  //!< optional saliency (score) of each voxel (usually based on Hessian eigenvalues)
               VectorContainer const *const *const *aaaafV,  //!< vector associated with each voxel
               TensorContainer ***aaaafDest,  //!< votes will be collected here
               Scalar const *const *const *aaafMaskSource=nullptr,  //!< ignore voxels in source where mask==0
               Scalar const *const *const *aaafMaskDest=nullptr,  //!< don't cast votes wherever mask==0
               bool detect_curves_not_surfaces=false, //!< do "sticks" represent curve tangents (instead of surface normals)?
               //Scalar saliency_threshold = 0.0,
               bool normalize=true, //!< normalize aaaafDest due to incomplete sums near boundaries?
               bool diagonalize_dest=false, //!< diagonalize each tensor in aaaafDest?
               ostream *pReportProgress=nullptr  //!< print progress to the user?
               )
  {
    assert(aaaafV);

    if (pReportProgress)
      *pReportProgress
        << " -- Attempting to allocate space for one more image.\n"
        << " -- (If this crashes your computer, find a computer with\n"
        << " --  more RAM and use \"ulimit\", OR use a smaller image.)\n";

    Scalar *afDenominator = nullptr;
    Scalar ***aaafDenominator = nullptr;

    if (normalize && aaafMaskSource) {
      Alloc3D(image_size,
              &afDenominator,
              &aaafDenominator);
    }

    // If the user did not specify an aaafSaliency[] array (if nullptr),
    // then we must create our own temporary saliency array.
    Scalar const *const *const *saliency_array = aaafSaliency;
    Scalar *_afSaliency = nullptr;
    Scalar ***_aaafSaliency = nullptr;

    if (! aaafSaliency) {
      // If the caller did not specify an aaafSaliency array, then
      // infer the saliency from the magnitude of aaaafV
      Alloc3D(image_size,
              &(_afSaliency),
              &(_aaafSaliency));
      for (Integer iz=0; iz<image_size[2]; iz++) {
        #pragma omp parallel for collapse(2)
        for (Integer iy=0; iy<image_size[1]; iy++) {
          for (Integer ix=0; ix<image_size[0]; ix++) {
            if ((! aaafMaskDest) || (aaafMaskDest[iz][iy][ix] == 0))
              continue;
            _aaafSaliency[iz][iy][ix] = Length3(aaaafV[iz][iy][ix]);
          }
        }
      }
      saliency_array = _aaafSaliency;
    } //if (! aaafSaliency)


    TVDenseStick(image_size,
                 saliency_array,
                 aaaafV,
                 aaaafDest,
                 aaafMaskSource,
                 aaafMaskDest,
                 detect_curves_not_surfaces,
                 //saliency_threshold,
                 aaafDenominator,
                 pReportProgress);


    // If any of the tensor voting sums were incomplete
    // (due to image boundaries, or mask boundaries).
    // then normalize the resulting magnitudes of the filter
    // at pixels located near the image boundaries (or mask boundaries).
    // OPTIONAL (probably not useful to most users)
    if (normalize) {
      if (pReportProgress)
        *pReportProgress << "  Normalizing the result of tensor voting...";

      if (aaafMaskSource) {

        assert(aaafDenominator);
        for (Integer iz=0; iz<image_size[2]; iz++) {
          #pragma omp parallel for collapse(2)
          for (Integer iy=0; iy<image_size[1]; iy++) {
            for (Integer ix=0; ix<image_size[0]; ix++) {
              if ((! aaafMaskDest) || (aaafMaskDest[iz][iy][ix] == 0))
                continue;
              if (aaafDenominator[iz][iy][ix] > 0.0) {
                assert(aaaafDest[iz][iy][ix]);
                for (int di = 0; di < 3; di++) {
                  for (int dj = di; dj < 3; dj++) {
                    aaaafDest[iz][iy][ix][ MapIndices_3x3_to_linear[di][dj] ]
                      /= aaafDenominator[iz][iy][ix];
                  }
                }
              }
            }
          }
        } //for (Integer iz=0; iz<image_size[2]; iz++)
        Dealloc3D(image_size,
                  &afDenominator,
                  &aaafDenominator);

        if (pReportProgress)
          *pReportProgress << "done." << endl;

      } // if (aaafMask)
      else {
        // THE UGLY CODE BELOW (THE NEXT ELSE-CLAUSE) IS UNNECESSARY,
        // BUT IT MAKES TENSOR-VOTING CODE ABOUT 10% FASTER.
        // IF YOU WANT TO MAKE THE CODE PRETTIER, DELETE THIS ELSE-CLAUSE
        // AND MAKE SURE "aaafDenominator" IS ALLWAYS ALLOCATED.
        assert(aaafDenominator == nullptr);
      
        // When no mask is supplied, 
        // If there is no mask, but the user wants the result to be normalized,
        // then we convolve the filter with the rectangular box. This is cheaper
        // because the convolution of a separable filter with a rectangular box 
        // shaped function is the product of the convolution with three 1-D 
        // functions which are 1 from 0..image_size[d], and 0 everywhere else.
        Integer halfwidth_single = halfwidth[0];
        assert(halfwidth_single == halfwidth[1]);
        assert(halfwidth_single == halfwidth[2]);
        Filter1D<Scalar, Integer> filter1d = 
          GenFilterGauss1D(sigma, halfwidth_single);
        Scalar *aafDenom_precomputed[3];
        for (int d=0; d<3; d++) {
          Scalar *afAllOnes = new Scalar [image_size[d]];
          aafDenom_precomputed[d] = new Scalar [image_size[d]];
          for (Integer i=0; i < image_size[d]; i++)
            afAllOnes[i] = 1.0;
          filter1d.Apply(image_size[d], afAllOnes, aafDenom_precomputed[d]);
          delete [] afAllOnes;
        }
        for (Integer iz = 0; iz < image_size[2]; iz++) {
          #pragma omp parallel for collapse(2)
          for (Integer iy = 0; iy < image_size[1]; iy++) {
            for (Integer ix = 0; ix < image_size[0]; ix++) {
              if ((! aaafMaskDest) || (aaafMaskDest[iz][iy][ix] == 0))
                continue;
              Scalar denominator = (aafDenom_precomputed[0][ix] *
                                    aafDenom_precomputed[1][iy] *
                                    aafDenom_precomputed[2][iz]);
              for (int di=0; di<3; di++) {
                for (int dj=0; dj<3; dj++) {
                  aaaafDest[iz][iy][ix][ MapIndices_3x3_to_linear[di][dj] ]
                    /= denominator;
                }
              }
            }
          }
        }
        // delete the array we created for storing the precomputed denominator:
        for (int d=0; d<3; d++)
          delete [] aafDenom_precomputed[d];
      } // else clause for "if (aaafMaskSource)"

    } //if (normalize)


    if (diagonalize_dest) {
      if (pReportProgress)
        *pReportProgress << "---- Diagonalizing Tensor Voting results ----" << endl;

      // Diagonalize the resulting tensor.
      // The resulting eigenvalues and eigenvectors can be analyzed by the caller.
      DiagonalizeHessianImage(image_size,
                              aaaafDest, //<--undiagonalized tensors (input)
                              aaaafDest, //<--diagonalized tensors (output)
                              aaafMaskDest,
                              selfadjoint_eigen3::DECREASING_EIVALS,
                              pReportProgress);

      // DEBUG: assert that all eigenvalues are nonnegative
      for (Integer iz=0; iz<image_size[2]; iz++) {
        for (Integer iy=0; iy<image_size[1]; iy++) {
          for (Integer ix=0; ix<image_size[0]; ix++) {
            if ((! aaafMaskDest) || (aaafMaskDest[iz][iy][ix] == 0.0))
              continue;
            // The eigenvalues are in the first 3 entries of aaaafDest[iz][iy][ix]
            Scalar *eivals = aaaafDest[iz][iy][ix];
            assert(eivals);
            for (int d=0; d<3; d++)
              assert(eivals[d] >= 0.0);
          }
        }
      }

    } // if (diagonalize_dest)

  } // TVDenseStick()







private:
  /// @brief  Perform dense stick-voting, using every voxel in the image
  ///         as a source, and collecting votes at every voxel
  ///         in the aaaafDest array.
  ///         This version of this function offers the ability to manually.
  ///         manage the aaafDenominator array (used for normalization).
  ///         Most users should use the other version of this function.
  void
  TVDenseStick(Integer const image_size[3],  //!< source image size
               Scalar const *const *const *aaafSaliency,  //!< saliency (score) of each voxel (usually based on Hessian eigenvalues)
               VectorContainer const *const *const *aaaafV,  //!< vector associated with each voxel
               TensorContainer ***aaaafDest,  //!< votes will be collected here
               Scalar const *const *const *aaafMaskSource=nullptr,  //!< ignore voxels in source where mask==0
               Scalar const *const *const *aaafMaskDest=nullptr,  //!< don't cast votes wherever mask==0
               bool detect_curves_not_surfaces=false,
               //Scalar saliency_threshold = 0.0,
               Scalar ***aaafDenominator=nullptr,
               ostream *pReportProgress=nullptr  //!< print progress to the user?
               )
  {
    assert(aaafSaliency);
    assert(aaaafV);
    assert(aaaafDest);

    //optional: count the number of voxels which can
    //          cast votes (useful for benchmarking)
    if (pReportProgress) {
      size_t n_salient = 0;
      size_t n_all = 0;
      for (Integer iz=0; iz<image_size[2]; iz++) {
        for (Integer iy=0; iy<image_size[1]; iy++) {
          for (Integer ix=0; ix<image_size[0]; ix++) {
            if (aaafMaskSource && (aaafMaskSource[iz][iy][ix] == 0))
              continue;
            n_all++;
            //if (aaafSaliency[iz][iy][ix] > saliency_threshold)
            if (aaafSaliency[iz][iy][ix] != 0.0)
              n_salient++;
          }
        }
      }
      *pReportProgress << "  (fraction of salient voxels = "
                       << (static_cast<double>(n_salient) / n_all)
                       << ")\n"
                       << "  (use aggressive thresholding to improve speed)"
                       << endl;
    }


    // First, initialize the arrays which will store the results with zeros.
    for (Integer iz=0; iz<image_size[2]; iz++)
      for (Integer iy=0; iy<image_size[1]; iy++)
        for (Integer ix=0; ix<image_size[0]; ix++)
          if (aaaafDest[iz][iy][ix])
            for (int di=0; di<3; di++)
              for (int dj=0; dj<3; dj++)
                aaaafDest[iz][iy][ix][ MapIndices_3x3_to_linear[di][dj] ] = 0.0;

    if (aaafDenominator) {
      // aaafDenominator[][][] keeps track of how much of the sum of
      // tensor-voting contributions was available at each voxel location.
      // (If some voxels were unavailable or outside the boundaries of 
      //  the image, they cannot cast votes at this voxel location.)
      //  This array keeps track of that.)  We should initialize this array too.
      for (Integer iz=0; iz<image_size[2]; iz++)
        for (Integer iy=0; iy<image_size[1]; iy++)
          for (Integer ix=0; ix<image_size[0]; ix++)
            aaafDenominator[iz][iy][ix] = 0.0;
    }

    // REMOVE THIS CRUFT
    //assert(pv);
    //assert(pV->nchannels() == 3);
    //pV->Resize(image_size, aaafMaskSource, pReportProgress);

    if (pReportProgress)
      *pReportProgress << "---- Begin Tensor Voting (dense, stick) ----\n"
                       << "  progress: processing plane#" << endl;

    // The mask should be 1 everywhere we want to consider, and 0 elsewhere.
    // Multiplying the density in the tomogram by the mask removes some of 
    // the voxels from consideration later on when we do the filtering.
    // (Later, we will adjust the weight of the average we compute when we
    //  apply the filter in order to account for the voxels we deleted now.)

    for (Integer iz=0; iz<image_size[2]; iz++) {

      if (pReportProgress)
        *pReportProgress << "  " << iz+1 << " / " << image_size[2] << "\n";

      #pragma omp parallel for collapse(2)
      for (Integer iy=0; iy<image_size[1]; iy++) {

        for (Integer ix=0; ix<image_size[0]; ix++) {

          // ------------ Uncomment whichever is faster: -----------
          // EITHER uncomment TVCastStickVotes() or TVReceiveStickVotes()
          // (VERSION1 OR VERSION2), BUT NOT BOTH

          // VERSION1:
          //           Have the voxel at ix,iy,iz cast votes at nearby voxels:

          TVCastStickVotes(ix, iy, iz,
                           image_size,
                           aaafSaliency,
                           aaaafV,
                           aaaafDest,
                           aaafMaskSource,
                           aaafMaskDest,
                           detect_curves_not_surfaces,
                           //saliency_threshold,
                           aaafDenominator);

          // VERSION 2: Have the voxel at ix,iy,iz receive votes
          //            from nearby voxels
          //
          //TVReceiveStickVotes(ix, iy, iz,
          //                    image_size,
          //                    aaafSaliency,
          //                    aaaafV,
          //                    aaaafDest,
          //                    aaafMaskSource,
          //                    aaafMaskDest,
          //                    detect_curves_not_surfaces,
          //                    //saliency_threshold,
          //                    (aaafDenominator
          //                     ? &(aaafDenominator[iz][iy][ix])
          //                     : nullptr));

        }
      }
    }

  } //TVDenseStick()


  /// @brief  Cast stick votes from one voxel to all nearby voxels
  void
  TVCastStickVotes(Integer ix,  //!< coordinates of the voter
                   Integer iy,  //!< coordinates of the voter
                   Integer iz,  //!< coordinates of the voter
                   Integer const image_size[3],
                   Scalar const *const *const *aaafSaliency, //!< saliency (score) of each voxel (usually calculated from Hessian eigenvalues)
                   VectorContainer const *const *const *aaaafV,  //!< vector associated with each voxel
                   TensorContainer ***aaaafDest,  //!< votes will be collected here
                   Scalar const *const *const *aaafMaskSource,  //!< ignore voxels in source where mask==0
                   Scalar const *const *const *aaafMaskDest,  //!< ignore voxels in dest where mask==0
                   bool detect_curves_not_surfaces = false,
                   //Scalar saliency_threshold = 0.0,
                   Scalar ***aaafDenominator = nullptr) const
  {
    assert(aaafSaliency);
    assert(aaaafV);
    assert(aaaafDest);

    Scalar saliency = aaafSaliency[iz][iy][ix];
    if (saliency == 0.0)
    //if (saliency <= saliency_threshold)
      return;

    Scalar mask_val = 1.0;
    if (aaafMaskSource) {
      mask_val = aaafMaskSource[iz][iy][ix];
      if (mask_val == 0.0)
        return;
    }

    Scalar n[3]; // the direction of the stick tensor (see below)
    for (int d=0; d<3; d++)
      n[d] = aaaafV[iz][iy][ix][d];

    //Note: The "filter_val" also is needed to calculate
    //      the denominator used in normalization.
    //      It is unusual to use a mask unless you intend
    //      to normalize the result later, but I don't enforce this

    for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++) {
      Integer iz_jz = iz+jz;
      if ((iz_jz < 0) || (image_size[2] <= iz_jz))
        continue;

      for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {
        Integer iy_jy = iy+jy;
        if ((iy_jy < 0) || (image_size[1] <= iy_jy))
          continue;

        for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {
          Integer ix_jx = ix+jx;
          if ((ix_jx < 0) || (image_size[0] <= ix_jx))
            continue;

          if ((aaafMaskDest) && (aaafMaskDest[iz_jz][iy_jy][ix_jx] == 0.0))
            continue;
          assert(aaaafDest[iz_jz][iy_jy][ix_jx]);

          // The function describing how the vote-strength falls off with
          // distance has been precomputed and is stored in
          // radial_decay_lookup.aaafH[jz][jy][jx];
          // In most tensor-voting implementations, this is a Gaussian.
          Scalar filter_val = radial_decay_lookup.aaafH[jz][jy][jx];
          if (aaafMaskSource)
            filter_val *= mask_val;

          Scalar decay_radial = filter_val;
          if (decay_radial == 0.0)
            continue;

          Scalar r[3];
          for (int d=0; d<3; d++)
            r[d] = aaaafDisplacement[jz][jy][jx][d];

          //
          //   .
          //   :\
          //   : \
          // ->:  \<- theta
          //   :   \
          //   :    \
          //   :     \
          //   :      \
          //   :       \
          //   ^        \
          //   |         \
          // n |          \
          //   |        .-' r
          //   |,,..--''
          //  0
          //
          // "n" = the direction of the stick tensor.  If the object being 
          //       detected is a surface, n is perpendicular to that surface.
          //       If the object is a curve, then n points tangent to the curve.
          // "r" = the position of the vote receiver relative to the voter
          //
          // "theta" = the angle of r relative to the plane perpendicular to n.
          //           (NOT the angle of r relative to n.  This is a confusing
          //            convention, but this is how it is normally defined.)

          Scalar sintheta = DotProduct3(r, n); //(sin() not cos(), see diagram)
          Scalar sinx2 = sintheta * 2.0;
          Scalar sin2 = sintheta * sintheta;
          Scalar cos2 = 1.0 - sin2;
          Scalar angle_dependence2 = cos2;
          if (detect_curves_not_surfaces) {
            // If we are detecting 1-dimensional curves (polymers, etc...)
            // instead of 2-dimensional surfaces (membranes, ...)
            // then the stick direction is assumed to be along the 
            // of the curve, (as opposed to perpendicular to the 2D surface).
            // As such
            angle_dependence2 = sin2;
          }

          Scalar decay_angular;

          switch(exponent) {
          case 2:
            decay_angular = angle_dependence2;
            break;
          case 4:
            decay_angular = angle_dependence2 * angle_dependence2;
            break;
          default:
            //Scalar angle_dependence = sqrt(angle_dependence2);
            decay_angular = pow(angle_dependence2, 0.5*exponent);
            break;
          }

          Scalar n_rotated[3];
          for (int d=0; d<3; d++) {
            if (detect_curves_not_surfaces)
              n_rotated[d] = n[d] - sinx2*r[d];
            else
              n_rotated[d] = sinx2*r[d] - n[d];
          }

          Scalar tensor_vote[3][3];
          for (int di = 0; di < 3; di++) {
            for (int dj = di; dj < 3; dj++) {
              tensor_vote[di][dj] = (saliency *
                                     decay_radial *
                                     decay_angular *
                                     n_rotated[di] * n_rotated[dj]);

              // OLD CODE: I used to implement aaaafDest as a 5-D array.
              //
              //aaaafDest[iz_jz][iy_jy][ix_jx][di][dj] += tensor_vote[di][dj];
              //
              // NEW CODE:
              // The ix_jx,iy_jy,iz_jz'th entry in aaaafDest should be
              // a 3x3 matrix. Since this matrix is symmentric, it contains
              // only 6 non-redundant entries.  Consequently there is no
              // need to store 9 numbers if only 6 are needed.
              // So I implemented a version of the 3x3 matrix which has
              // only 6 entries, arranged in a 1-D array of size 6.
              // To access these entries, use "MapIndices_3x3_to_linear[][]".

              aaaafDest[iz_jz][iy_jy][ix_jx][MapIndices_3x3_to_linear[di][dj]]
                += tensor_vote[di][dj];
            }
          }

          if (aaafDenominator)
            aaafDenominator[iz_jz][iy_jy][ix_jx] += filter_val;

        } // for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++)
      } // for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++)
    } // for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++)


  } // TVCastStickVotes()




  /// @brief  Receive stick votes from all voxels near a chosen "receiver" voxel
  void
  TVReceiveStickVotes(Integer ix,  //!< coordinates of the receiver voxel
                      Integer iy,  //!< coordinates of the receiver voxel
                      Integer iz,  //!< coordinates of the receiver voxel
                      Integer const image_size[3],
                      Scalar const *const *const *aaafSaliency, //!< saliency (score) of each voxel (usually calculated from Hessian eigenvalues)
                      VectorContainer const *const *const *aaaafV,  //!< vector associated with each voxel
                      TensorContainer ***aaaafDest,  //!< votes will be collected here
                      Scalar const *const *const *aaafMaskSource,  //!< ignore voxels in source where mask==0
                      Scalar const *const *const *aaafMaskDest,  //!< ignore voxels in dest where mask==0
                      bool detect_curves_not_surfaces = false,
                      //Scalar saliency_threshold = 0.0,
                      Scalar *pDenominator = nullptr) const
  {
    assert(aaafSaliency);
    assert(aaaafV);
    assert(aaaafDest);

    if (aaafMaskDest && (aaafMaskDest[iz][iy][ix] == 0.0))
      return;

    Scalar denominator = 0.0;

    for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++) {
      Integer iz_jz = iz-jz;
      if ((iz_jz < 0) || (image_size[2] <= iz_jz))
        continue;

      for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {
        Integer iy_jy = iy-jy;
        if ((iy_jy < 0) || (image_size[1] <= iy_jy))
          continue;

        for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {
          Integer ix_jx = ix-jx;
          if ((ix_jx < 0) || (image_size[0] <= ix_jx))
            continue;

          // The function describing how the vote-strength falls off with
          // distance has been precomputed and is stored in
          // radial_decay_lookup.aaafH[jz][jy][jx];
          // In most tensor-voting implementations, this is a Gaussian.
          Scalar filter_val = radial_decay_lookup.aaafH[jz][jy][jx];

          if (aaafMaskSource) {
            Scalar mask_val = aaafMaskSource[iz_jz][iy_jy][ix_jx];
            if (mask_val == 0.0)
              continue;
            filter_val *= mask_val;
          }
          //Note: The "filter_val" also is needed to calculate
          //      the denominator used in normalization.
          //      It is unusual to use a mask unless you intend
          //      to normalize the result later, but I don't enforce this

          Scalar saliency = aaafSaliency[iz_jz][iy_jy][ix_jx];
          //if (saliency <= saliency_threshold)
          if (saliency == 0.0)
            continue;

          Scalar decay_radial = filter_val;
          if (decay_radial == 0.0)
            continue;

          Scalar r[3];
          Scalar n[3];
          for (int d=0; d<3; d++) {
            r[d] = aaaafDisplacement[jz][jy][jx][d];
            n[d] = aaaafV[iz_jz][iy_jy][ix_jx][d];
          }

          //
          //   .
          //   :\
          //   : \
          // ->:  \<- theta
          //   :   \
          //   :    \
          //   :     \
          //   :      \
          //   :       \
          //   ^        \
          //   |         \
          // n |          \
          //   |        .-' r
          //   |,,..--''
          //  0
          //
          // "n" = the direction of the stick tensor.  If the object being 
          //       detected is a surface, n is perpendicular to that surface.
          //       If the object is a curve, then n points tangent to the curve.
          // "r" = the position of the vote receiver relative to the voter
          //
          // "theta" = the angle of r relative to the plane perpendicular to n.
          //           (NOT the angle of r relative to n.  This is a confusing
          //            convention, but this is how it is normally defined.)

          Scalar sintheta = DotProduct3(r, n); //(sin() not cos(), see diagram)
          Scalar sinx2 = sintheta * 2.0;
          Scalar sin2 = sintheta * sintheta;
          Scalar cos2 = 1.0 - sin2;
          Scalar angle_dependence2 = cos2;
          if (detect_curves_not_surfaces) {
            // If we are detecting 1-dimensional curves (polymers, etc...)
            // instead of 2-dimensional surfaces (membranes, ...)
            // then the stick direction is assumed to be along the 
            // of the curve, (as opposed to perpendicular to the 2D surface).
            // As such
            angle_dependence2 = sin2;
          }

          Scalar decay_angular;

          switch(exponent) {
          case 2:
            decay_angular = angle_dependence2;
            break;
          case 4:
            decay_angular = angle_dependence2 * angle_dependence2;
            break;
          default:
            //Scalar angle_dependence = sqrt(angle_dependence2);
            decay_angular = pow(angle_dependence2, 0.5*exponent);
            break;
          }

          Scalar n_rotated[3];
          for (int d=0; d<3; d++) {
            if (detect_curves_not_surfaces)
              n_rotated[d] = n[d] - sinx2*r[d];
            else
              n_rotated[d] = sinx2*r[d] - n[d];
          }

          Scalar tensor_vote[3][3];
          for (int di = 0; di < 3; di++) {
            for (int dj = di; dj < 3; dj++) {
              tensor_vote[di][dj] = (saliency *
                                     decay_radial *
                                     decay_angular *
                                     n_rotated[di] * n_rotated[dj]);

              // OLD CODE: I used to implement aaaafDest as a 5-D array.
              //
              //aaaafDest[iz][iy][ix][di][dj] += tensor_vote[di][dj];
              //
              // NEW CODE:
              // The ix,iy,iz'th entry in aaaafDest should be a 3x3 matrix.
              // Since this matrix is symmentric, it contains only 6
              // non-redundant entries.  Consequently there is no need to
              // store 9 numbers if only 6 are needed.
              // So I implemented a version of the 3x3 matrix which has
              // only 6 entries, arranged in a 1-D array of size 6.
              // To access these entries, use "MapIndices_3x3_to_linear[][]".

              aaaafDest[iz][iy][ix][ MapIndices_3x3_to_linear[di][dj] ]
                += tensor_vote[di][dj];
            }
          }

          if (pDenominator)
            denominator += filter_val;
        } // for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++)
      } // for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++)
    } // for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++)

    if (pDenominator)
      *pDenominator = denominator;

  } // TVReceiveStickVotes()



  void swap(TV3D<Scalar,Integer,VectorContainer,TensorContainer> &other)
  {
    std::swap(sigma, other.sigma);
    std::swap(exponent, other.exponent);
    std::swap(radial_decay_lookup, other.radial_decay_lookup);
    std::swap(aafDisplacement, other.aafDisplacement);
    std::swap(aaaafDisplacement, other.aaaafDisplacement);
    std::swap(halfwidth, other.halfwidth);
    std::swap(array_size, other.array_size);
  }


  TV3D<Scalar,Integer,VectorContainer,TensorContainer>&
    operator = (TV3D<Scalar,Integer,VectorContainer,TensorContainer> source)
  {
    this->swap(source);
    return *this;
  }

private:

  void Init() {
    sigma = 0.0;
    halfwidth[0] = -1;
    halfwidth[1] = -1;
    halfwidth[2] = -1;
    array_size[0] = -1;
    array_size[1] = -1;
    array_size[2] = -1;
    aafDisplacement = nullptr;
    aaaafDisplacement = nullptr;
  }


  void Resize(Integer set_halfwidth[3]) {
    for (int d=0; d<3; d++) {
      array_size[d] = 2*halfwidth[d] + 1;
    }

    Scalar sigmas[3] = {sigma, sigma, sigma};
    radial_decay_lookup =
      GenFilterGenGauss3D(sigmas,
                          static_cast<Scalar>(2.0),
                          halfwidth);

    AllocDisplacement();
    PrecalcDisplacement(aaaafDisplacement);
  }


  void DeallocDisplacement() {
    if (aaaafDisplacement) {
      for (int d=0; d<3; d++) {
        array_size[d] = 2*halfwidth[d] + 1;
      }
      //shift pointers back to normal
      for (Integer iz = -halfwidth[2]; iz <= halfwidth[2]; iz++) {
        for (Integer iy = -halfwidth[1]; iy <= halfwidth[1]; iy++) {
          aaaafDisplacement[iz][iy] -= halfwidth[0];
        }
        aaaafDisplacement[iz] -= halfwidth[1];
      }
      aaaafDisplacement -= halfwidth[2];
      Dealloc3D(array_size,
                &aafDisplacement,
                &aaaafDisplacement);
    }
  }

  void AllocDisplacement() {
    if (aaaafDisplacement)
      Dealloc3D(array_size,
                &aafDisplacement,
                &aaaafDisplacement);
    Alloc3D(array_size,
            &aafDisplacement,
            &aaaafDisplacement);
    //shift pointers to enable indexing from i = -halfwidth .. +halfwidth
    for (Integer iz = 0; iz < array_size[2]; iz++) {
      for (Integer iy = 0; iy < array_size[1]; iy++) {
        aaaafDisplacement[iz][iy] += halfwidth[0];
      }
      aaaafDisplacement[iz] += halfwidth[1];
    }
    aaaafDisplacement += halfwidth[2];
  }

  void PrecalcDisplacement(array<Scalar, 3> ***aaaafDisplacement) {
    // pre-compute the normalized radius unit vector
    for (Integer iz = -halfwidth[2]; iz <= halfwidth[2]; iz++) {
      for (Integer iy = -halfwidth[1]; iy <= halfwidth[1]; iy++) {
        for (Integer ix = -halfwidth[0]; ix <= halfwidth[0]; ix++) {
          Scalar length = sqrt(ix*ix + iy*iy + iz*iz);
          if (length == 0)
            length = 1.0;
          aaaafDisplacement[iz][iy][ix][0] = ix / length;
          aaaafDisplacement[iz][iy][ix][1] = iy / length;
          aaaafDisplacement[iz][iy][ix][2] = iz / length;
        }
      }
    }
  }
}; // class TV3D



} //namespace visfd



#endif //#ifndef _FEATURE_HPP
