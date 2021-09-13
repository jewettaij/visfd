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
#include <feature.hpp>
#include "filter3d_variants.hpp"



namespace visfd {


/// @brief   The "CalcMomentTensor()" function was intended to be a variant
/// of the "CalcHessian()" function.  I HAVE NOT TESTED THIS VERSION.
/// (In this version, I apply the derivative to the Gaussian filter
///  before applying the filter, ...instead of applying the Gaussian filter
///  first and then taking finite differences afterwards.  Speculation:
///  If the width of the object being detected is not much more than
///  3-voxel wide, the 3-voxel wide differences used in the other
///  implementation are a large source of error.  I'm not sure this approach
///  is better though.)
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


  Scalar ***aaafNorm = Alloc3D<Scalar>(image_size);

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

    Scalar ***aaafIx = Alloc3D<Scalar>(image_size);
    Scalar ***aaafIy = Alloc3D<Scalar>(image_size);
    Scalar ***aaafIz = Alloc3D<Scalar>(image_size);

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
    Dealloc3D(aaafIx);
    Dealloc3D(aaafIy);
    Dealloc3D(aaafIz);
  } //if (pFirstMoment)



  if (aaaaf2ndMoment) {
    if (pReportProgress)
      *pReportProgress << "\n"
        " ------ Calculating the average of nearby voxels: ------\n";
    // P = original image (after subtracting average nearby intensities):

    Scalar ***aaafP = Alloc3D<Scalar>(image_size);

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

    Scalar ***aaafIxx = Alloc3D<Scalar>(image_size);
    Scalar ***aaafIyy = Alloc3D<Scalar>(image_size);
    Scalar ***aaafIzz = Alloc3D<Scalar>(image_size);
    Scalar ***aaafIxy = Alloc3D<Scalar>(image_size);
    Scalar ***aaafIyz = Alloc3D<Scalar>(image_size);
    Scalar ***aaafIxz = Alloc3D<Scalar>(image_size);

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

    Dealloc3D(aaafP);

    Dealloc3D(aaafIxx);
    Dealloc3D(aaafIyy);
    Dealloc3D(aaafIzz);
    Dealloc3D(aaafIxy);
    Dealloc3D(aaafIxz);
    Dealloc3D(aaafIyz);
  } //if (pSecondMoment)

  Dealloc3D(aaafNorm);

} // CalcMomentTensor()




/// @brief Find scale-invariant blobs in the image as a function of diameter.
///        This variant detects blobs and performs non-max suppression.
///
/// Algorithm described in:
///    Lindeberg,T., Int. J. Comput. Vision., 30(2):77-116, (1998)
/// This version refers to blobs by their diameter (instead of "sigma").
/// This version can discard blobs which overlap with existing blobs.
/// (This is sometimes called "non-max suppression".)
/// This function is not intended for public use.

template<typename Scalar>
static void
BlobDogNM(int const image_size[3], //!<source image size
          Scalar const *const *const *aaafSource,   //!<source image
          Scalar const *const *const *aaafMask,     //!<ignore voxels where mask==0
          const vector<Scalar>& blob_diameters, //!< blob widths to try, ordered
          //optional arguments:
          vector<array<Scalar,3> > *pva_minima_crds=nullptr, //!< if not nullptr, stores blob minima x,y,z coords here
          vector<array<Scalar,3> > *pva_maxima_crds=nullptr, //!< if not nullptr, stores blob maxima x,y,z coords here
          vector<Scalar> *pv_minima_diameters=nullptr, //!< if not nullptr, stores the corresponding width for that minima
          vector<Scalar> *pv_maxima_diameters=nullptr, //!< if not nullptr, stores the corresponding width for that maxima
          vector<Scalar> *pv_minima_scores=nullptr, //!< if not nullptr, stores the blob's score?
          vector<Scalar> *pv_maxima_scores=nullptr, //!< (score = intensity after filtering)
          const Scalar aspect_ratio[3]=nullptr, //!<multiply blob_sigma by different numbers in the X,Y,Z directions (default:1,1,1)
          //the following optional parameters are usually left with default values
          Scalar delta_sigma_over_sigma=0.02,//!< param for approximating LOG with DOG
          Scalar truncate_ratio=2.5,      //!< how many sigma before truncating?
          Scalar minima_threshold=0.5,  //!< discard blobs with unremarkable scores
          Scalar maxima_threshold=0.5,  //!< discard blobs with unremarkable scores
          bool   use_threshold_ratios=true, //!< threshold=ratio*best_score?
          Scalar sep_ratio_thresh=1.0,          //!< minimum radial separation between blobs
          Scalar nonmax_max_overlap_large=1.0,  //!< maximum volume overlap with larger blob
          Scalar nonmax_max_overlap_small=1.0,  //!< maximum volume overlap with smaller blob
          // optional arguments
          ostream *pReportProgress = nullptr, //!< report progress to the user?
          Scalar ****aaaafI = nullptr //!<optional: preallocated memory for filtered images
          )
{

  vector<array<Scalar,3> > minima_crds; //store minima blob x,y,z coords here
  vector<array<Scalar,3> > maxima_crds; //store minima blob x,y,z coords here
  vector<Scalar>  minima_diameters;     //corresponding width for that minima
  vector<Scalar>  maxima_diameters;     //corresponding width for that maxima
  vector<Scalar>  minima_scores;        //store the score of each blob minima
  vector<Scalar>  maxima_scores;        //store the score of each blob maxima
  if (pva_minima_crds == nullptr)
    pva_minima_crds = &minima_crds;
  if (pva_maxima_crds == nullptr)
    pva_maxima_crds = &maxima_crds;
  if (pv_minima_diameters == nullptr)
    pv_minima_diameters = &minima_diameters;
  if (pv_maxima_diameters == nullptr)
    pv_maxima_diameters = &maxima_diameters;
  if (pv_minima_scores == nullptr)
    pv_minima_scores = &minima_scores;
  if (pv_maxima_scores == nullptr)
    pv_maxima_scores = &maxima_scores;

  Scalar default_aspect_ratio[3] = {1.0, 1.0, 1.0};
  const Scalar *_aspect_ratio = aspect_ratio;
  if (_aspect_ratio == nullptr)
    _aspect_ratio = default_aspect_ratio;

  BlobDogD(image_size,
           aaafSource,
           aaafMask,
           blob_diameters,
           pva_minima_crds,
           pva_maxima_crds,
           pv_minima_diameters,
           pv_maxima_diameters,
           pv_minima_scores,
           pv_maxima_scores,
           _aspect_ratio,
           delta_sigma_over_sigma,
           truncate_ratio,
           minima_threshold,
           maxima_threshold,
           use_threshold_ratios,
           pReportProgress,
           aaaafI);

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

  DiscardOverlappingBlobs(*pva_minima_crds,
                          *pv_minima_diameters,
                          *pv_minima_scores,
                          sep_ratio_thresh,
                          nonmax_max_overlap_large,
                          nonmax_max_overlap_small,
                          SORT_INCREASING,
                          pReportProgress);

  if (pReportProgress)
    *pReportProgress << "done --\n"
                     << "--- Discarding overlapping maxima blobs ---\n";

  DiscardOverlappingBlobs(*pva_maxima_crds,
                          *pv_maxima_diameters,
                          *pv_maxima_scores,
                          sep_ratio_thresh,
                          nonmax_max_overlap_large,
                          nonmax_max_overlap_small,
                          SORT_DECREASING,
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

template<typename Scalar>
static
void
_BlobDogNM(int const image_size[3], //!<source image size
           Scalar const *const *const *aaafSource,   //!<source image
           Scalar const *const *const *aaafMask,     //!<ignore voxels where mask==0
           const vector<Scalar>& blob_diameters, //!<list of diameters to try (ordered)
           //optional arguments:
           vector<array<Scalar,3> > *pva_minima_crds=nullptr, //!< if not nullptr, stores blob minima x,y,z coords here
           vector<array<Scalar,3> > *pva_maxima_crds=nullptr, //!< if not nullptr, stores blob maxima x,y,z coords here
           vector<Scalar> *pv_minima_diameters=nullptr, //!< if not nullptr, stores the corresponding width for that minima
           vector<Scalar> *pv_maxima_diameters=nullptr, //!< if not nullptr, stores the corresponding width for that maxima
           vector<Scalar> *pv_minima_scores=nullptr, //!< if not nullptr, stores the blob's score?
           vector<Scalar> *pv_maxima_scores=nullptr, //!< (score = intensity after filtering)
           const Scalar aspect_ratio[3]=nullptr, //!<multiply blob_sigma by different numbers in the X,Y,Z directions (default:1,1,1)
           //the following optional parameters are usually left with default values
           Scalar delta_sigma_over_sigma=0.02, //!<difference in Gauss widths parameter
           Scalar filter_truncate_ratio=2.5,   //!<how many sigma before truncating?
           Scalar filter_truncate_threshold=0.02, //!<decay in filter before truncating
           Scalar minima_threshold=0.0,    //!<discard unremarkable minima
           Scalar maxima_threshold=0.0,    //!<discard unremarkable maxima
           bool use_threshold_ratios=true, //!<threshold=ratio*best_score ?
           Scalar sep_ratio_thresh=1.0,          //!<minimum radial separation between blobs
           Scalar nonmax_max_overlap_large=1.0,  //!<maximum volume overlap with larger blob
           Scalar nonmax_max_overlap_small=1.0,  //!<maximum volume overlap with smaller blob
           ostream *pReportProgress = nullptr,
           Scalar ****aaaafI = nullptr //!<optional: preallocated memory for filtered images
           )
{
  
  if (filter_truncate_ratio <= 0) {
    assert(filter_truncate_threshold > 0.0);
    //    filter_truncate_threshold = exp(-(1/2)*filter_truncate_ratio^2);
    //    -> filter_truncate_ratio^2 = -2*log(filter_truncate_threshold)
    filter_truncate_ratio = sqrt(-2*log(filter_truncate_threshold));
  }

  Scalar default_aspect_ratio[3] = {1.0, 1.0, 1.0};
  const Scalar *_aspect_ratio = aspect_ratio;
  if (_aspect_ratio == nullptr)
    _aspect_ratio = default_aspect_ratio;

  BlobDogNM(image_size,
            aaafSource,
            aaafMask,
            blob_diameters,
            pva_minima_crds,
            pva_maxima_crds,
            pv_minima_diameters,
            pv_maxima_diameters,
            pv_minima_scores,
            pv_maxima_scores,
            _aspect_ratio,
            delta_sigma_over_sigma,
            filter_truncate_ratio,
            minima_threshold,
            maxima_threshold,
            use_threshold_ratios,
            sep_ratio_thresh,
            nonmax_max_overlap_large,
            nonmax_max_overlap_small,
            pReportProgress,
            aaaafI);

} //_BlobDogNM(...,filter_truncate_ratio,filter_truncate_threshold,...)




} //namespace visfd



#endif //#ifndef _FEATURE_VARIANTS_HPP
