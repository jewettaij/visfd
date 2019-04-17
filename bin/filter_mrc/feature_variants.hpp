#ifndef _FEATURE_VARIANTS_HPP
#define _FEATURE_VARIANTS_HPP

/// THE CODE IN THIS FILE WAS NOT INTENDED FOR PUBLIC USE
/// This file contains several variants of functions found in <feature.hpp>.
/// These versions differ only in that they accept slightly different format
/// of arguments, usually providing multiple (often redundant) ways to
/// specify parameters (such as the filter-window-width).
/// The goal was to give the end user multiple different convenient options for
/// specifying these parameters, along with some reasonable default values.
/// (I don't expect end users to ever want to mess with these defaults.)
/// Unfortunately, these new functions have uglier, clumsier argument lists
/// (which is why I don't include them in <feature.hpp>).


#include <cstring>
#include <ostream>
#include <vector>
#include <tuple>
#include <set>
#include <cassert>
#include <queue>
#include <cmath>
#include <limits>
using namespace std;
#include <visfd.hpp>
#include "filter3d_variants.hpp"



namespace visfd {



/// @brief Find scale-invariant blobs in the image as a function of diameter.
///        This variant detects blobs and performs non-max suppression.
///
/// Algorithm described in:
///    Lindeberg,T., Int. J. Comput. Vision., 30(2):77-116, (1998)
/// This version refers to blobs by their diameter (instead of "sigma").
/// This version can discard blobs which overlap with existing blobs.
/// (This is sometimes called "non-max suppression".)
/// This function is not intended for public use.

template<class Scalar>
static void
BlobDogNM(int const image_size[3], //!<source image size
          Scalar const *const *const *aaafSource,   //!<source image
          Scalar const *const *const *aaafMask,     //!<ignore voxels where mask==0
          const vector<Scalar>& blob_diameters, //!< blob widths to try, ordered
          vector<array<Scalar,3> >& minima_crds, //!<store minima x,y,z coords here
          vector<array<Scalar,3> >& maxima_crds, //!<store maxima x,y,z coords here
          vector<Scalar>& minima_diameters, //!< corresponding width for that minima
          vector<Scalar>& maxima_diameters, //!< corresponding width for that maxima
          vector<Scalar>& minima_scores, //!< what was the blob's score?
          vector<Scalar>& maxima_scores, //!< (score = intensity after filtering)
          Scalar delta_sigma_over_sigma,//!< param for approximating LOG with DOG
          Scalar truncate_ratio,      //!< how many sigma before truncating?
          Scalar minima_threshold=0.5,  //!< discard blobs with unremarkable scores
          Scalar maxima_threshold=0.5,  //!< discard blobs with unremarkable scores
          bool   use_threshold_ratios=true, //!< threshold=ratio*best_score?
          Scalar sep_ratio_thresh=1.0,          //!< minimum radial separation between blobs
          Scalar nonmax_max_overlap_large=1.0,  //!< maximum volume overlap with larger blob
          Scalar nonmax_max_overlap_small=1.0,  //!< maximum volume overlap with smaller blob
          // optional arguments
          ostream *pReportProgress = NULL, //!< report progress to the user?
          Scalar ****aaaafI = NULL, //!<preallocated memory for filtered images
          Scalar **aafI = NULL     //!<preallocated memory for filtered images
          )
{

  BlobDogD(image_size,
           aaafSource,
           aaafMask,
           blob_diameters,
           minima_crds,
           maxima_crds,
           minima_diameters,
           maxima_diameters,
           minima_scores,
           maxima_scores,
           delta_sigma_over_sigma,
           truncate_ratio,
           minima_threshold,
           maxima_threshold,
           use_threshold_ratios,
           pReportProgress,
           aaaafI,
           aafI);

  bool discard_overlapping_blobs =
    ((sep_ratio_thresh > 0.0) ||
     (nonmax_max_overlap_small < 1.0) ||
     (nonmax_max_overlap_large < 1.0));   

  if (! discard_overlapping_blobs)
    return;

  if (pReportProgress)
    *pReportProgress << "----------- Removing overlapping blobs -----------\n" << endl;


  if (pReportProgress)
    *pReportProgress << "--- Discarding overlapping minima blobs ---\n";

  DiscardOverlappingBlobs(minima_crds,
                          minima_diameters,
                          minima_scores,
                          PRIORITIZE_LOW_SCORES,
                          sep_ratio_thresh,
                          nonmax_max_overlap_large,
                          nonmax_max_overlap_small,
                          pReportProgress);

  if (pReportProgress)
    *pReportProgress << "done --\n"
                     << "--- Discarding overlapping maxima blobs ---\n";

  DiscardOverlappingBlobs(maxima_crds,
                          maxima_diameters,
                          maxima_scores,
                          PRIORITIZE_HIGH_SCORES,
                          sep_ratio_thresh,
                          nonmax_max_overlap_large,
                          nonmax_max_overlap_small,
                          pReportProgress);

} //BlobDogNM()


/// @brief Find scale-invariant blobs in the image as a function of diameter.
///        In this minor variant, the user can specify the filter window width
///        either in units of sigma, or by specifying the decay threshold.
///        This function was not intended for public use.
///
/// Algorithm described in:
///    Lindeberg,T., Int. J. Comput. Vision., 30(2):77-116, (1998)
/// This version refers to blobs by their diameter (instead of "sigma").
/// This version can discard blobs which overlap with existing blobs.
/// (This is sometimes called "non-max suppression".)

template<class Scalar>
static
void
_BlobDogNM(int const image_size[3], //!<source image size
           Scalar const *const *const *aaafSource,   //!<source image
           Scalar const *const *const *aaafMask,     //!<ignore voxels where mask==0
           const vector<Scalar>& blob_diameters, //!<list of diameters to try (ordered)
           vector<array<Scalar,3> >& minima_crds, //!<store minima x,y,z coords here
           vector<array<Scalar,3> >& maxima_crds, //!<store maxima x,y,z coords here
           vector<Scalar>& minima_diameters, //!<corresponding radius for that minima
           vector<Scalar>& maxima_diameters, //!<corresponding radius for that maxima
           vector<Scalar>& minima_scores, //!<what was the blobs score?
           vector<Scalar>& maxima_scores, //!<(score = intensity after filtering)
           Scalar delta_sigma_over_sigma, //!<difference in Gauss widths parameter
           Scalar filter_truncate_ratio,     //!<how many sigma before truncating?
           Scalar filter_truncate_threshold, //!<decay in filter before truncating
           Scalar minima_threshold,       //!<discard unremarkable minima
           Scalar maxima_threshold,       //!<discard unremarkable maxima
           bool use_threshold_ratios=true, //!<threshold=ratio*best_score ?
           Scalar sep_ratio_thresh=1.0,          //!<minimum radial separation between blobs
           Scalar nonmax_max_overlap_large=1.0,  //!<maximum volume overlap with larger blob
           Scalar nonmax_max_overlap_small=1.0,  //!<maximum volume overlap with smaller blob
           ostream *pReportProgress = NULL,
           Scalar ****aaaafI = NULL, //!<preallocated memory for filtered images
           Scalar **aafI = NULL     //!<preallocated memory for filtered images
           )
{
  
  if (filter_truncate_ratio <= 0) {
    assert(filter_truncate_threshold > 0.0);
    //    filter_truncate_threshold = exp(-(1/2)*filter_truncate_ratio^2);
    //    -> filter_truncate_ratio^2 = -2*log(filter_truncate_threshold)
    filter_truncate_ratio = sqrt(-2*log(filter_truncate_threshold));
  }

  BlobDogNM(image_size,
            aaafSource,
            aaafMask,
            blob_diameters,
            minima_crds,
            maxima_crds,
            minima_diameters,
            maxima_diameters,
            minima_scores,
            maxima_scores,
            delta_sigma_over_sigma,
            filter_truncate_ratio,
            minima_threshold,
            maxima_threshold,
            use_threshold_ratios,
            sep_ratio_thresh,
            nonmax_max_overlap_large,
            nonmax_max_overlap_small,
            pReportProgress,
            aaaafI,
            aafI);

} //_BlobDogNM(...,filter_truncate_ratio,filter_truncate_threshold,...)




} //namespace visfd



#endif //#ifndef _FEATURE_VARIANTS_HPP
