#ifndef _FILTER1D_H
#define _FILTER1D_H

#include <cstring>
using namespace std;


template<class RealNum >
inline RealNum ABS(RealNum x) { return ((x<0.0) ? -x: x); }

template<class RealNum >
inline RealNum SQR(RealNum x) { return x*x; }

template<class RealNum >
inline RealNum MAX(RealNum x, RealNum y) { return ((x<y) ? y : x); }



template<class RealNum, class Integer>

class Filter1D {
public:
  RealNum *afH;
  Integer halfwidth; //distance from the filter center to its edge, in pixels
  Integer array_size; //size of the array in x,y directions (in voxels)


  // Apply the filter to tomographic data in the "afSource" array
  // Save the results in the "afDest" array.  (A "mask" is optional.)
  // All arrays are 1D and assumed to be the same size
  void Apply(Integer const size_source,
             RealNum *afSource,
             RealNum *afDest,
             RealNum const *afMask = NULL,
             bool normalize = false) const
             //bool precompute_mask_times_source = true) const
  {

    // Apply the filter to the original tomogram data. (Store in afSource)
    //        ___
    //        \
    // g(i) = /__  h(j) * f(i-j)
    //         j
    //
    // where: f(i) is the original density of the tomogram at position ix,iy
    //        h(j) is the filter (smoothing function)
    //       
    //  Note on mask functions:
    //     When summing over j, we ignore contributions from voxels where
    //     mask(i-j) is zero.  We don't count them in the average.
    //     Because h(j) is not necessarily normalized, g(i) is later divided by
    //     the area under the curve h(j)*mask(i-j)    (as a function of j)
    //     

    // The mask should be 1 everywhere we want to consider, and 0 elsewhere.
    // Multiplying the density in the tomogram by the mask removes some of 
    // the voxels from consideration later on when we do the filtering.
    // (Later, we will adjust the weight of the average we compute when we
    //  apply the filter in order to account for the voxels we deleted now.)
    // Precomupting the mask product is faster but modifies the source image.

    //if (afMask && precompute_mask_times_source)
    if (afMask) 
      for (int ix=0; ix<size_source; ix++)
        afSource[ix] *= afMask[ix];

    for (Integer ix=0; ix<size_source; ix++) {

      if ((afMask) && (afMask[ix] == 0.0)) {
        afDest[ix] = 0.0;
        continue;
      }
          
      RealNum g = 0.0;
      RealNum denominator = 0.0;

      for (Integer jx=-halfwidth; jx<=halfwidth; jx++) {

        if ((ix-jx < 0) || (size_source <= ix-jx))
          continue;

        RealNum delta_g = afH[jx+halfwidth] * afSource[ix-jx];

        //if (! precompute_mask_times_source)
        //  delta_g *= afMask[ix-jx];

        g += delta_g;
          // Note: We previously applied the mask by multiplying
          //          aDensity[ix]
          //         by afMask[ix]   (if present)
        if (normalize) {
          if (afMask)
            denominator += afMask[ix-jx] * afH[jx+halfwidth];
          else
            denominator += afH[jx+halfwidth];
        }
                                          
        // Note: If there were no mask, and if the filter is normalized
        // then denominator=1 always, and we could skip the line above.
      }

      if (normalize) {
        if (denominator > 0.0)
          g /= denominator;
        else
          //Otherwise, this position lies outside the mask region.
          g = 0.0;
      }

      afDest[ix] = g;
    }
  } //Apply()

  void Alloc(Integer set_halfwidth) {
    halfwidth = set_halfwidth;
    array_size = 1 + 2*halfwidth;
    afH = new RealNum [array_size];
    for (int i = 0; i < array_size; i++)
      afH[i] = -1.0e38; //(if uninitiliazed memory read, we will know)
  }

  void Dealloc() {
    if (afH)
      delete [] afH;
    halfwidth = -1;
    array_size = -1;
  }


  void Resize(Integer set_halfwidth) {
    Dealloc();
    Alloc(set_halfwidth);
  }


  inline Filter1D<RealNum, Integer>&
    operator = (const Filter1D<RealNum, Integer>& source) {
    Resize(source.halfwidth); // allocates and initializes af and afH
    //for(Int ix=-halfwidth; ix<=halfwidth; ix++)
    //  afH[ix] = source.afH[ix];
    // Use memcpy() instead:
    memcpy(afH,
           source.afH,
           (1+2*halfwidth) * sizeof(RealNum));
  } // operator = ()


  void Init() {
    halfwidth = -1;
    afH = NULL;
  }

  Filter1D() {
    Init();
  }

  Filter1D(Integer halfwidth) {
    Init();
    Resize(halfwidth);
  }

  Filter1D(const Filter1D<RealNum, Integer>& f) {
    Init();
    *this = f;
  }

  ~Filter1D() {
    Dealloc();
  }


  void Normalize() {
    // Make sure the sum of the filter weights is 1
    RealNum total = 0.0;
    for (Integer ix=-halfwidth; ix<=halfwidth; ix++)
      total += afH[ix+halfwidth];
    for (Integer ix=-halfwidth; ix<=halfwidth; ix++)
      afH[ix+halfwidth] /= total;
  }


}; // class Filter1D




template<class RealNum>
// GenFilterGauss1D generates a 1-D filter and fills its array with values
// corresponding to a normalized Gaussian evaluated at evenly spaced intervals.
// The caller must specify the "sigma" parameter (width of the Gaussian,
// in units of pixels/voxels), in addition to the "halfwidth" parameter, which
// indicates the number of entries in the array (in units of pixels/voxels).
Filter1D<RealNum, int>
GenFilterGauss1D(RealNum sigma,  // The "sigma" paramgeter in the Gaussian
                 int halfwidth,  // number of entries in the filter array / 2
                 ostream *pReportProgress = NULL)
{
  Filter1D<RealNum, int> filter(halfwidth);

  RealNum sum = 0.0;
  for (int i=-halfwidth; i<=halfwidth; i++) {
    if (sigma == 0.0) //(don't crash when sigma=0)
      filter.afH[i+halfwidth] = ((i == 0) ? 1.0 : 0.0);
    else
      //filter.afH[i+halfwidth] = exp(-0.5*(i*i)/(sigma*sigma));
      filter.afH[i+halfwidth] = exp(-(i*i)/(2.0*sigma*sigma));
    sum += filter.afH[i+halfwidth];
  }

  //Normalize:
  for (int i=-halfwidth; i<=halfwidth; i++) {
    filter.afH[i+halfwidth] /= sum;

    //FOR DEBUGGING REMOVE EVENTUALLY:
    if (pReportProgress)
      *pReportProgress <<"Gauss1D: afH["<<i<<"] = "
                       << filter.afH[i+halfwidth] << endl;
  }
  return filter;
} //GenFilterGauss1D(sigma, halfwidth)




// REMOVE THIS CRUFT EVENTUALLY:
//
//template<class RealNum>
//// This function generates a 1-D filter and fills its array with values
//// corresponding to a normalized Gaussian evaluated at evenly spaced intervals
//// The caller must specify the "sigma" parameter (width of the Gaussian).
//// and a fractional number between 0 and 1 which indicate how far the
//// filter can decay.  Only voxels whose Gaussian intensity decays by less
//// than this threshold (relative to the central peak) will be kept.
//Filter1D<RealNum, int>
//GenFilterGauss1DThresh(RealNum sigma,
//                       RealNum truncate_threshold,
//                       ostream *pReportProgress = NULL)
//{
//  // How wide should the filter be?
//  int halfwidth;
//  assert(truncate_threshold > 0.0);
//
//  // Choose the filter domain window based on the "truncate_threshold"
//
//  //    truncate_threshold = exp(-0.5*(halfwidth/sigma)^2);
//  //    -> (halfwidth/sigma)^2 = -2*log(truncate_threshold)
//  halfwidth = floor(sigma * sqrt(-2*log(truncate_threshold)));
//
//  return GenFilterGauss1D(sigma/sqrt(2.0), halfwidth, pReportProgress);
//} //GenFilterGauss1D(sigma, truncate_threshold)





#endif //#ifndef _FILTER1D_H
