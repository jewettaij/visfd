#ifndef _FILTER3D_H
#define _FILTER3D_H

#include <cstring>
#include <ostream>
#include <vector>
using namespace std;
#include <filter1d.h>  // defines "Filter1D" (used in "ApplyGauss()")



template<class RealNum, class Integer>
RealNum AverageArr(Integer const array_size[3],
                   RealNum ***aaafH,
                   RealNum ***aaafW = NULL) 
{
  double total = 0.0;
  double denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        double h = aaafH[iz][iy][ix];
        if (aaafW) {
          double w = aaafW[iz][iy][ix];
          h *= w;
          denom += w;
        }
        else
          denom += 1.0;
        total += h;
      }
    }
  }
  return static_cast<RealNum>(total / denom);
} //void AverageArr()



template<class RealNum, class Integer>
RealNum AveSqrArr(Integer const array_size[3],
                  RealNum ***aaafH,
                  RealNum ***aaafW = NULL) 
{
  double total = 0.0;
  double denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        double h = aaafH[iz][iy][ix];
        h *= h;

        if (aaafW) {
          double w = aaafW[iz][iy][ix];
          h *= w;
          denom += w;
        }
        else
          denom += 1.0;
        total += h;
      }
    }
  }
  return static_cast<RealNum>(sqrt(total / denom));
} //void AveSqrArr()


template<class RealNum, class Integer>
RealNum StdDevArr(Integer const array_size[3],
                  RealNum ***aaafH,
                  RealNum ***aaafW = NULL) 
{
  double ave = AverageArr(array_size, aaafH, aaafW);
  double total = 0.0;
  double denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        double h = (aaafH[iz][iy][ix] - ave);
        h *= h;

        if (aaafW) {
          double w = aaafW[iz][iy][ix];
          h *= w;
          denom += w;
        }
        else
          denom += 1.0;
        total += h;
      }
    }
  }
  return static_cast<RealNum>(sqrt(total / denom));
} //void StdDev()



template<class RealNum, class Integer>
RealNum SumSqrArr(Integer const array_size[3],
                  RealNum ***aaafH,
                  RealNum ***aaafW = NULL) 
{
  double total = 0.0;
  double denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        double h = aaafH[iz][iy][ix];
        h *= h;
        if (aaafW)
          h *= aaafW[iz][iy][ix];
        total += h;
      }
    }
  }
  return static_cast<RealNum>(total);
} //void SumSqrArr()



template<class RealNum, class Integer>
RealNum SumArr(Integer const array_size[3],
               RealNum ***aaafH,
               RealNum ***aaafW = NULL) 
{
  double total = 0.0;
  double denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        double h = aaafH[iz][iy][ix];
        total += h;
        if (aaafW)
          h *= aaafW[iz][iy][ix];
      }
    }
  }
  return static_cast<RealNum>(total);
} //void SumArr()



template<class RealNum, class Integer>
void AddScalarArr(RealNum offset,
                  Integer const array_size[3],
                  RealNum ***aaafH)
{
  for (Integer iz = 0; iz < array_size[2]; iz++)
    for (Integer iy = 0; iy < array_size[1]; iy++)
      for (Integer ix = 0; ix < array_size[0]; ix++)
        aaafH[iz][iy][ix] += offset;
}



template<class RealNum, class Integer>
void MultiplyScalarArr(RealNum scale,
                       Integer const array_size[3],
                       RealNum ***aaafH)
{
  for (Integer iz = 0; iz < array_size[2]; iz++)
    for (Integer iy = 0; iy < array_size[1]; iy++)
      for (Integer ix = 0; ix < array_size[0]; ix++)
        aaafH[iz][iy][ix] *= scale;
}



template<class RealNum, class Integer>
void MinMaxArr(RealNum& min,
               RealNum& max,
               Integer const array_size[3],
               RealNum ***aaafH,
               RealNum ***aaafMask = NULL) 
{
  bool first = true;
  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        if ((aaafMask) and (aaafMask[iz][iy][ix] == 0.0))
          continue;
        if (first || (aaafH[iz][iy][ix] > max))
          max = aaafH[iz][iy][ix];
        if (first || (aaafH[iz][iy][ix] < min))
          min = aaafH[iz][iy][ix];
        first = false;
      }
    }
  }
} //void MinMaxArr()




template<class RealNum, class Integer>
void
IntensityHistogramArr(RealNum **paHistX, //pointer to an array of intensities
                      long long **paHistY,//number of voxels with that intensity
                      Integer &nbins,     //specify number of bins (if positive)
                      RealNum &bin_width, //alternatively, suggest bin_width
                      Integer const array_size[3],
                      RealNum ***aaafI,
                      RealNum ***aaafW = NULL) 
{
  RealNum hmin = 0.0;
  RealNum hmax = -1.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        if (hmin > hmax) {
          hmin = aaafI[iz][iy][ix];
          hmax = aaafI[iz][iy][ix];
        }
        if (aaafI[iz][iy][ix] < hmin)
          hmin = aaafI[iz][iy][ix];
        if (aaafI[iz][iy][ix] > hmax)
          hmax = aaafI[iz][iy][ix];
      }
    }
  }

  if (nbins > 1) {
    bin_width = (hmax - hmin) / (nbins - 1);
  }
  else if (bin_width > 0) {
    nbins = 1 + ceil((hmax - hmin)/bin_width);
  }

  assert(paHistX);
  *paHistX = new RealNum [nbins];
  for (int i=0; i<nbins; i++)
    (*paHistX)[i] = hmin + i*bin_width;

  assert(paHistY);
  *paHistY = new long long [nbins];
  for (int i=0; i < nbins; i++)
    (*paHistY)[i] = 0;
  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        RealNum h = aaafI[iz][iy][ix];
        if ((! aaafW) || (aaafW[iz][iy][ix] != 0.0)) {
          RealNum h = aaafI[iz][iy][ix];
          Integer i = floor( (h - hmin)/bin_width );
          assert((0 <= i) && (i < nbins));
          (*paHistY)[i]++;
        }
      }
    }
  }
} //void IntensityHistogramArr()





template<class RealNum, class Integer>

class Filter3D {
public:
  RealNum *af;
  RealNum ***aaafH;
  Integer halfwidth[3]; //num pixels from filter center to edge in x,y,z directions
  Integer array_size[3]; //size of the array in x,y,z directions (in voxels)

  RealNum ApplyToVoxel(Integer ix,  // voxel's position in the x,y,z directions
                       Integer iy,
                       Integer iz,
                       Integer const size_source[3], // size of the source image
                       RealNum ***aaafSource,        // source image
                       RealNum ***aaafMask = NULL,   // consider only some voxels?
                       bool normalize = false) const // normalize incomplete sums?

  // Apply a filter to a single voxel from the source image (aaafSource)
  // Arguments:
  // ix,iy,iz the voxel's position in that 3-D image
  // size_source[] contains size of the source image (in the x,y,z directions)
  // aaafSource[][][] is the source image
  // aaafMask[][][] = 0 for voxels we should exclude from consideration
  // normalize = true if we want to normalize the weights considered in the sum
  // (useful if the sum was not complete due to some voxels being masked out).
  {
    RealNum g = 0.0;
    RealNum denominator = 0.0;

    for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++) {

      if ((iz-jz < 0) || (size_source[2] <= iz-jz))
        continue;

      for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {

        if ((iy-jy < 0) || (size_source[1] <= iy-jy))
          continue;

        for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {

          if ((ix-jx < 0) || (size_source[0] <= ix-jx))
            continue;

          RealNum delta_g = 
            aaafH[jz+halfwidth[2]]
                 [jy+halfwidth[1]]
                 [jx+halfwidth[0]]
            *
            aaafSource[iz-jz][iy-jy][ix-jx];

          //if (! precompute_mask_times_source)
          //  delta_g *= aaafMask[iz-jz][iy-jy][ix-jx];

          g += delta_g;
                       // Note: We previously applied the mask by multiplying
                       //          aaafDensity[iz][iy][ix]
                       //          by aaafMask[iz][iy][ix]   (if present)
          if (normalize) {
            if (aaafMask)
              denominator += aaafMask[iz-jz][iy-jy][ix-jx] *
                aaafH[jz+halfwidth[2]]
                     [jy+halfwidth[1]]
                     [jx+halfwidth[0]];
            else
              denominator += aaafH[jz+halfwidth[2]]
                                  [jy+halfwidth[1]]
                                  [jx+halfwidth[0]];
          }
          // Note: If there is mask, and if the filter is normalized
          // then denominator=1 always (except in the edge of the image),
          // and we could probably skip the "if (normalize) {...}" code above.
        }
      }
    }

    if (normalize) {
      if (denominator > 0.0)
        g /= denominator;
      else
        //Otherwise, this position lies outside the mask region.
        g = 0.0;
    }
    return g;
  } // ApplyToVoxel()





  // Apply the filter to tomographic data in the "aaafSource" array
  // Save the results in the "aaafDest" array.  (A "mask" is optional.)
  // All arrays are 3D and assumed to be the same size
  void Apply(Integer const size_source[3],
             RealNum ***aaafSource,
             RealNum ***aaafDest,
             RealNum ***aaafMask = NULL,
             bool normalize = false,
             //bool precompute_mask_times_source = true,
             ostream *pReportProgress = NULL
             ) const
  {

    // Apply the filter to the original tomogram data. (Store in aaafSource)
    //        ___
    //        \
    // g(i) = /__  h(j) * f(i-j)
    //         j
    //
    // where: f(i) is the original density of the tomogram at position ix,iy,iz
    //        h(j) is the filter (smoothing function)
    //       
    //  Note on mask functions:
    //     When summing over j, we ignore contributions from voxels where
    //     mask(i-j) is zero.  We don't count them in the average.
    //     Because h(j) is not necessarily normalized, g(i) is later divided by
    //     the area under the curve h(j)*mask(i-j)    (as a function of j)
    //     

    if (pReportProgress)
      *pReportProgress << "  progress: processing plane#" << endl;

    // The mask should be 1 everywhere we want to consider, and 0 elsewhere.
    // Multiplying the density in the tomogram by the mask removes some of 
    // the voxels from consideration later on when we do the filtering.
    // (Later, we will adjust the weight of the average we compute when we
    //  apply the filter in order to account for the voxels we deleted now.)
    // Precomupting the mask product is faster but modifies the source image.

    //if (aaafMask && precompute_mask_times_source)
    if (aaafMask) 
      for (int iz=0; iz<size_source[2]; iz++)
        for (int iy=0; iy<size_source[1]; iy++)
          for (int ix=0; ix<size_source[0]; ix++)
            aaafSource[iz][iy][ix] *= aaafMask[iz][iy][ix];

    for (Integer iz=0; iz<size_source[2]; iz++) {

      if (pReportProgress)
        *pReportProgress << "  " << iz+1 << " / " << size_source[2] << "\n";

      for (Integer iy=0; iy<size_source[1]; iy++) {

        for (Integer ix=0; ix<size_source[0]; ix++) {

          // Calculate the effect of the filter on
          // the voxel located at position ix,iy,iz

          if ((aaafMask) && (aaafMask[iz][iy][ix] == 0.0)) {
            aaafDest[iz][iy][ix] = 0.0;
            continue;
          }

          aaafDest[iz][iy][ix] =
            ApplyToVoxel(ix, iy, iz,
                         size_source,
                         aaafSource,
                         aaafMask,
                         normalize);
        }
      }
    }
  } // Apply()



  #ifndef DISABLE_TEMPLATE_MATCHING

  RealNum
  TemplateError(Integer ix, //voxel's position in the x,y,z directions
                Integer iy,
                Integer iz,
                Integer const size_source[3], // size of the source image
                RealNum ***aaafSource,
                RealNum intensity_ratio, //how much to scale template intensities
                RealNum ***aaafW, //weights used in all summations
                RealNum template_compare_exponent,
                RealNum ***aaafMask = NULL,
                //bool normalize = false,
                //bool precompute_mask_times_source = true,
                ostream *pReportProgress = NULL
                ) const 
  {
      RealNum g = 0.0;
      RealNum denominator = 0.0;

      for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++) {

      if ((iz-jz < 0) || (size_source[2] <= iz-jz))
        continue;

      for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {

        if ((iy-jy < 0) || (size_source[1] <= iy-jy))
          continue;

        for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {

          if ((ix-jx < 0) || (size_source[0] <= ix-jx))
            continue;

          // DELETE THIS CRUFT:
          //// Kluge:
          //// If the filter h = 0 here, then it probably means we are outside
          //// the radius we were supposed to consider.  In that case,
          //// ignore that voxel. The goal is to avoid consideration of "corner"
          //// voxels which lie outside this radius of consideration.
          //// The voxel brightnesses here are not relevant and would 
          //// screw up (and overwhelm) the sum of differences so far
          ////
          //// (Perhaps I should do something more robust here
          //// but in practice, voxels where h=0 occur inside this radius 
          //// so rarely that I think it's safe to ignore them when they are 0.
          //// I may need to handle this in a better way sometime in the future
          //// especially if we give the user more control over the filter.)
          //if (aaafH[jz+halfwidth[2]]
          //         [jy+halfwidth[1]]
          //         [jx+halfwidth[0]] == 0)
          //  continue;


          RealNum delta_g = 
            (intensity_ratio * aaafH[jz+halfwidth[2]]
                                    [jy+halfwidth[1]]
                                    [jx+halfwidth[0]]
             -
             aaafSource[iz-jz][iy-jy][ix-jx]);

          if (template_compare_exponent == 2.0)
            delta_g *= delta_g;
          else
            delta_g = pow(delta_g, template_compare_exponent);

          //if (! precompute_mask_times_source)
          delta_g *= 
              aaafW[jz+halfwidth[2]]
                   [jy+halfwidth[1]]
                   [jx+halfwidth[0]];


          if (aaafMask)     //<--rarely used. I may delete "aaafMask" later
            delta_g *= aaafMask[iz-jz][iy-jy][ix-jx];

          g += delta_g;

          //if (normalize) {
          if (aaafMask)     //<--rarely used. I may delete "aaafMask" later
            denominator +=
              aaafMask[iz-jz][iy-jy][ix-jx]
              *
              aaafW[jz+halfwidth[2]]
                   [jy+halfwidth[1]]
                   [jx+halfwidth[0]];

            // DELETE THIS CRUFT:
            //*SQR(aaafH[jz+halfwidth[2]]
            //          [jy+halfwidth[1]]
            //          [jx+halfwidth[0]]);

          else
            //denominator += 1.0;
            denominator += aaafW[jz+halfwidth[2]]
                                [jy+halfwidth[1]]
                                [jx+halfwidth[0]];

            // DELETE THIS CRUFT:
            // += SQR(aaafH[jz+halfwidth[2]]
            //             [jy+halfwidth[1]]
            //             [jx+halfwidth[0]]);

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
  } // TemplateError()

  

  void
  ScanTemplateError(Integer const size_source[3], // size of the source image
                    RealNum ***aaafSource,
                    RealNum ***aaafDest,
                    RealNum ***aaafC, //how much to scale template intensities
                    RealNum ***aaafW, //weights used in all summations
                    RealNum template_compare_exponent,
                    RealNum ***aaafMask = NULL,
                    //bool normalize = false,
                    //bool precompute_mask_times_source = true,
                    ostream *pReportProgress = NULL
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
                          template_compare_exponent,
                          NULL); //aaafMask);
                          //normalize);
        }
      }
    }
  } //ScanTemplateError()

  #endif //#ifndef DISABLE_TEMPLATE_MATCHING



  void Alloc(Integer const set_halfwidth[3]) {
    //Integer array_size[3];
    for(Integer d=0; d < 3; d++) {
      halfwidth[d] = set_halfwidth[d];
      array_size[d] = 1 + 2*halfwidth[d];
    }
    Alloc3D(array_size, &af, &aaafH);
    for (int iz = 0; iz < array_size[2]; iz++)
      for (int iy = 0; iy < array_size[1]; iy++)
        for (int ix = 0; ix < array_size[0]; ix++)
          aaafH[iz][iy][ix] = -1.0e38; //(if uninitiliazed memory read, we will know)
  }

  void Dealloc() {
    //Integer array_size[3];
    for(Integer d=0; d < 3; d++) {
      array_size[d] = 1 + 2*halfwidth[d];
      halfwidth[d] = -1;
    }
    Dealloc3D(array_size, &af, &aaafH);
    for(Integer d=0; d < 3; d++)
      array_size[d] = -1;
  }

  void Resize(Integer const set_halfwidth[3]) {
    Dealloc();
    Alloc(set_halfwidth);
  }

  inline Filter3D<RealNum, Integer>&
    operator = (const Filter3D<RealNum, Integer>& source) {
    Resize(source.halfwidth); // allocates and initializes af and aaafH
    //for(Int iz=-halfwidth[2]; iz<=halfwidth[2]; iz++)
    //  for(Int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
    //    for(Int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
    //      aaafH[iz][iy][ix] = source.aaafH[iz][iy][ix];
    // Use memcpy() instead:
    memcpy(af,
           source.af,
           ((1+2*halfwidth[0]) * (1+2*halfwidth[1]) * (1+2*halfwidth[2]))
           *sizeof(RealNum));
  } // operator = ()

  void Init() {
    halfwidth[0] = -1;
    halfwidth[1] = -1;
    halfwidth[2] = -1;
    af = NULL;
    aaafH = NULL;
  }
  

  Filter3D(Integer const set_halfwidth[3]) {
    Init();
    Resize(set_halfwidth);
  }

  Filter3D() {
    Init();
  }

  Filter3D(const Filter3D<RealNum, Integer>& f) {
    Init();
    *this = f;
  }
  

  ~Filter3D() {
    Dealloc();
  }

  RealNum Average(RealNum ***aaafW) {
    return AverageArr(array_size, aaafH, aaafW);
  }

  RealNum AverageSqr(RealNum ***aaafW) {
    return AveSqrArr(array_size, aaafH, aaafW);
  }

  RealNum StdDev(RealNum ***aaafW) {
    return StdDevArr(array_size, aaafH, aaafW);
  }

  RealNum SumSqr(RealNum ***aaafW) {
    return SumSqrArr(array_size, aaafH, aaafW);
  }

  RealNum Sum(RealNum ***aaafW) {
    return SumArr(array_size, aaafH, aaafW);
  }

  void AddScalar(RealNum offset) {
    AddScalarArr(offset, array_size, aaafH);
  }

  void MultiplyScalar(RealNum scale) {
    MultiplyScalarArr(scale, array_size, aaafH);
  }


  void Normalize() {
    // Make sure the sum of the filter weights is 1
    RealNum total = 0.0;
    for (Integer iz=-halfwidth[2]; iz<=halfwidth[2]; iz++)
      for (Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
        for (Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
          total +=
            aaafH[iz+halfwidth[2]][iy+halfwidth[1]][ix+halfwidth[0]];
    assert(total > 0.0);
    for (Integer iz=-halfwidth[2]; iz<=halfwidth[2]; iz++)
      for (Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
        for (Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
          aaafH[iz+halfwidth[2]][iy+halfwidth[1]][ix+halfwidth[0]]
            /= total;
  }

}; // class Filter







template<class RealNum>
// Create a 3-D filter and fill it with a "generalized Gaussian" function:
//   h(x,y,z) = A*exp(-r^m)
// where   r  = sqrt((x/s_x)^2 + (y/s_y)^2 + (z/s_z)^2)
//   and   A  is determined by normalization of the discrete sum
// Note: "A" is equal to the value stored in the middle of the array,
//       The caller can determine what "A" is by looking at this value.
Filter3D<RealNum, int>
GenFilterGenGauss3D(RealNum width[3],    //"s_x", "s_y", "s_z" parameters
                    RealNum m_exp,       //"m" exponent parameter
                    int truncate_halfwidth[3],
                    RealNum *pA=NULL,    //optional:report A coeff to user
                    ostream *pReportEquation = NULL)
{
  RealNum window_threshold = 1.0;
  for (int d=0; d<3; d++) {
    RealNum h = ((width[d]>0) ? exp(-pow(truncate_halfwidth[d]/width[d], m_exp)) : 1.0);
    if (h < window_threshold)
      window_threshold = h;
  }
  Filter3D<RealNum, int> filter(truncate_halfwidth);
  RealNum total = 0;
  for (int iz=-filter.halfwidth[2]; iz<=filter.halfwidth[2]; iz++) {
    for (int iy=-filter.halfwidth[1]; iy<=filter.halfwidth[1]; iy++) {
      for (int ix=-filter.halfwidth[0]; ix<=filter.halfwidth[0]; ix++) {
        RealNum r = sqrt(SQR(ix/width[0])+SQR(iy/width[1])+SQR(iz/width[2]));
        RealNum h = ((r>0) ? exp(-pow(r, m_exp)) : 1.0);
        if (ABS(h) < window_threshold)
          h = 0.0; //This eliminates corner entries which fall below threshold
                   //(and eliminates anisotropic artifacts due to these corners)
                   //There's no reason to keep any entries less than min value.
        filter.aaafH[iz+filter.halfwidth[2]]
                    [iy+filter.halfwidth[1]]
                    [ix+filter.halfwidth[0]] = h;
        total += h;
      }
    }
  }
  // normalize:
  for (int iz=-filter.halfwidth[2]; iz<=filter.halfwidth[2]; iz++) {
    for (int iy=-filter.halfwidth[1]; iy<=filter.halfwidth[1]; iy++) {
      for (int ix=-filter.halfwidth[0]; ix<=filter.halfwidth[0]; ix++) {
        filter.aaafH[iz+filter.halfwidth[2]]
                    [iy+filter.halfwidth[1]]
                    [ix+filter.halfwidth[0]] /= total;

        //FOR DEBUGGING REMOVE EVENTUALLY:
        if (pReportEquation)
          *pReportEquation << "GenGauss3D:" //<< window_threshold
                           << " aaafH["<<iz<<"]["<<iy<<"]["<<ix<<"] = "
                           << filter.aaafH[iz+filter.halfwidth[2]]
                                          [iy+filter.halfwidth[1]]
                                          [ix+filter.halfwidth[0]] << endl;
      }
    }
  }

  // The coefficient in front of the Gaussian ("A")
  // equals the height of its central peak,
  // which is located in the middle of the array
  RealNum A = filter.aaafH[filter.halfwidth[0]]
                          [filter.halfwidth[1]]
                          [filter.halfwidth[2]];

  if (pA) {
    *pA = A;
  }

  if (pReportEquation) {
    *pReportEquation << " Filter Used:\n"
      " h(x,y,z)   = A*exp(-((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2)^(m/2))\n"
      "  ... where      A = " << A << "\n"
      "                 m = " << m_exp << "\n" 
      "   (a_x, a_y, a_z) = "
                     << "(" << width[0]
                     << " " << width[1]
                     << " " << width[2] << ")\n";
    *pReportEquation << " You can plot a slice of this function\n"
                     << "     in the X direction using:\n"
      " draw_filter_1D.py -gauss " << A << " " << width[0] << endl;
    *pReportEquation << " and in the Y direction using:\n"
      " draw_filter_1D.py -gauss " << A << " " << width[1] << endl;
    *pReportEquation << " and in the Z direction using:\n"
      " draw_filter_1D.py -gauss " << A << " " << width[2] << endl;
    *pReportEquation << " You can plot a slice of this function\n"
                     << "     in the X direction using:\n"
      " draw_filter_1D.py -ggauss " << A
                     << " " << width[0]
                     << " " << m_exp << endl;
    *pReportEquation << " and in the Y direction using:\n"
      " draw_filter_1D.py -ggauss " << A
                     << " " << width[1]
                     << " " << m_exp << endl;
    *pReportEquation << " and in the Z direction using:\n"
      " draw_filter_1D.py -ggauss " << A
                     << " " << width[2]
                     << " " << m_exp << endl;
  }
  return filter;
} //GenFilterGenGauss3D(width, m_exp, truncate_halfwidth)




template<class RealNum>
Filter3D<RealNum, int>
GenFilterGenGauss3D(RealNum width[3],          //"s_x", "s_y", "s_z" parameters
                    RealNum m_exp,             //"m" parameter in formula
                    RealNum filter_cutoff_ratio,
                    RealNum *pA=NULL,    //optional:report A coeff to user
                    ostream *pReportProgress = NULL)
{
  // choose the width of the filter window based on the filter_cutoff_ratio
  int truncate_halfwidth[3];
  int ix = 0;

  for (int d=0; d<3; d++) {
    truncate_halfwidth[d] = floor(width[d]*filter_cutoff_ratio);
  }
  return GenFilterGenGauss3D(width,
                             m_exp,
                             truncate_halfwidth,
                             pA,
                             pReportProgress);
} //GenFilterGenGauss3D(width, m_exp, filter_cutoff_ratio)






// (The next function was not intended for public use)

// Create a 3-D filter and fill it with a difference of (generalized) Gaussians:
// This version requires that the caller has already created individual
// filters for the two gaussians.
// All this function does is subtract one filter from the other (and rescale).
template<class RealNum>
Filter3D<RealNum, int> 
_GenFilterDogg3D(RealNum width_a[3],  //"a" parameter in formula
                 RealNum width_b[3],  //"b" parameter in formula
                 RealNum m_exp,  //"m" parameter in formula
                 RealNum n_exp,  //"n" parameter in formula
                 Filter3D<RealNum, int>& filter_A, //filters for the two
                 Filter3D<RealNum, int>& filter_B, //gaussians
                 RealNum *pA=NULL, //optional:report A,B coeffs to user
                 RealNum *pB=NULL,
                 ostream *pReportEquation = NULL)
{
  RealNum A, B;
  //A, B = height of the central peak
  //       The central peak is located in the middle of the filter's array
  //       (at position "halfwidth")
  A = filter_A.aaafH[filter_A.halfwidth[0]]
                    [filter_A.halfwidth[1]]
                    [filter_A.halfwidth[2]];
  B = filter_B.aaafH[filter_B.halfwidth[0]]
                    [filter_B.halfwidth[1]]
                    [filter_B.halfwidth[2]];


  // The "difference of gaussians" filter is the difference between
  // these two (generalized) gaussian filters.
  int halfwidth[3];
  halfwidth[0] = MAX(filter_A.halfwidth[0], filter_B.halfwidth[0]);
  halfwidth[1] = MAX(filter_A.halfwidth[1], filter_B.halfwidth[1]);
  halfwidth[2] = MAX(filter_A.halfwidth[2], filter_B.halfwidth[2]);
  Filter3D<RealNum, int> filter(halfwidth);

  //FOR DEBUGGING REMOVE EVENTUALLY
  if (pReportEquation)
    *pReportEquation << "Array of 3D filter entries:" << endl;


  for (int iz=-halfwidth[2]; iz<=halfwidth[2]; iz++) {
    for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
      for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
        filter.aaafH[iz+halfwidth[2]]
                    [iy+halfwidth[1]]
                    [ix+halfwidth[0]] = 0.0;

        // The two filters may have different widths, so we have to check
        // that ix,iy and iz lie within the domain of these two filters before
        // adding or subtracting their values from the final GDOG filter.
        if (((-filter_A.halfwidth[0]<=ix) && (ix<=filter_A.halfwidth[0])) &&
            ((-filter_A.halfwidth[1]<=iy) && (iy<=filter_A.halfwidth[1])) &&
            ((-filter_A.halfwidth[2]<=iz) && (iz<=filter_A.halfwidth[2])))

          filter.aaafH[iz+halfwidth[2]]
                      [iy+halfwidth[1]]
                      [ix+halfwidth[0]]
            +=
            filter_A.aaafH[iz+filter_A.halfwidth[2]]
                          [iy+filter_A.halfwidth[1]]
                          [ix+filter_A.halfwidth[0]];  //   /  (A-B);  COMMENTING OUT

        // COMMENTING OUT: (The factor of 1/(A-B) insures that the central peak has height 1)


        if (((-filter_B.halfwidth[0]<=ix) && (ix<=filter_B.halfwidth[0])) &&
            ((-filter_B.halfwidth[1]<=iy) && (iy<=filter_B.halfwidth[1])) &&
            ((-filter_B.halfwidth[2]<=iz) && (iz<=filter_B.halfwidth[2])))

          filter.aaafH[iz+halfwidth[2]]
                      [iy+halfwidth[1]]
                      [ix+halfwidth[0]]
            -=
            filter_B.aaafH[iz+filter_B.halfwidth[2]]
                          [iy+filter_B.halfwidth[1]]
                          [ix+filter_B.halfwidth[0]];  //   /  (A-B);  COMMENTING OUT


        //FOR DEBUGGING REMOVE EVENTUALLY
        if (pReportEquation)
          *pReportEquation << "GenDogg: aaafH["<<iz<<"]["<<iy<<"]["<<ix<<"] = "
                           << filter.aaafH[iz+halfwidth[2]]
                                          [iy+halfwidth[1]]
                                          [ix+halfwidth[0]] << endl;
        //*pReportEquation << aaafH[iz+halfwidth[2]]
        //                         [iy+halfwidth[1]]
        //                         [ix+halfwidth[0]];
        //if (ix == halfwidth[0]) pReportEquation << "\n"; else pReportEquation << " ";

      } // for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
    } // for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
  } // for (int iz=-halfwidth[2]; iz<=halfwidth[2]; iz++) {


  //FOR DEBUGGING REMOVE EVENTUALLY
  if (pReportEquation)
    *pReportEquation << "\n";

  // COMMENTING OUT the factor of 1/(A-B):
  //A = A/(A-B);
  //B = B/(A-B);

  if (pA && pB) {
    *pA = A; // Rescale A and B numbers returned to the caller
    *pB = B; // (because we divided the array entries by (A-B) earlier)
  }

  if (pReportEquation) {
    *pReportEquation << " Filter Used:\n"
      " h(x,y,z)   = h_a(x,y,z) - h_b(x,y,z)\n"
      " h_a(x,y,z) = A*exp(-((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2)^(m/2))\n"
      " h_b(x,y,z) = B*exp(-((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2)^(n/2))\n"
      "  ... where      A = " << A << "\n"
      "                 B = " << B << "\n" 
      "                 m = " << m_exp << "\n"
      "                 n = " << n_exp << "\n" 
      "   (a_x, a_y, a_z) = "
                     << "(" << width_a[0]
                     << " " << width_a[1]
                     << " " << width_a[2] << ")\n"
      "   (b_x, b_y, b_z) = "
                     << "(" << width_b[0]
                     << " " << width_b[1]
                     << " " << width_b[2] << ")\n";
    *pReportEquation << " You can plot a slice of this function\n"
                     << "     in the X direction using:\n"
      " draw_filter_1D.py -dogg " << A << " " << B
                     << " " << width_a[0] << " " << width_b[0]
                     << " " << m_exp << " " << n_exp << endl;
    *pReportEquation << " and in the Y direction using:\n"
      " draw_filter_1D.py -dogg " << A << " " << B
                     << " " << width_a[1] << " " << width_b[1]
                     << " " << m_exp << " " << n_exp << endl;
    *pReportEquation << " and in the Z direction using:\n"
      " draw_filter_1D.py -dogg " << A << " " << B
                     << " " << width_a[2] << " " << width_b[2]
                     << " " << m_exp << " " << n_exp << endl;
  } //if (pReportEquation)
  
  return filter;
} //_GenFilterDogg3D()




template<class RealNum>
// Create a 3-D filter and fill it with a difference of (generalized) Gaussians:
Filter3D<RealNum, int> 
GenFilterDogg3D(RealNum width_a[3],   //"a" parameter in formula
                RealNum width_b[3],   //"b" parameter in formula
                RealNum m_exp,        //"m" parameter in formula
                RealNum n_exp,        //"n" parameter in formula
                int halfwidth[3],
                RealNum *pA=NULL,     //optional:report A,B coeffs to user
                RealNum *pB=NULL,     //optional:report A,B coeffs to user
                ostream *pReportEquation = NULL) //optional: print params used?
{
  Filter3D<RealNum, int> filter_A =
    GenFilterGenGauss3D(width_a,      //"a_x", "a_y" gaussian width parameters
                        m_exp,        //"m" exponent parameter
                        halfwidth);

  Filter3D<RealNum, int> filter_B =
    GenFilterGenGauss3D(width_b,      //"b_x", "b_y" gaussian width parameters
                        n_exp,        //"n" exponent parameter
                        halfwidth);

  return _GenFilterDogg3D(width_a,
                          width_b,
                          m_exp,
                          n_exp,
                          filter_A, filter_B,
                          pA,
                          pB,
                          pReportEquation);
} //GenFilterDogg3D(...halfwidth...)





// (The next function was not intended for public use)
template<class RealNum>
RealNum
_ApplyGauss3D(int const image_size[3], 
              RealNum ***aaafSource,
              RealNum ***aaafDest,
              RealNum ***aaafMask,
              Filter1D<float, int> aFilter[3],
              bool normalize = true,
              ostream *pReportProgress = NULL)
{
  assert(aaafSource);
  assert(aaafDest);

  // Optional
  // The "A" 3-D Gaussian coefficient is the product of the
  // 1-D Gaussian coefficients in the X,Y,Z directions.
  // Those coefficients happen to equal the value of the 1-D
  // Gaussian evaluated at its peak, which is stored in the central entry at
  // "halfwidth".  (The 1-D filter arrays have size equal to 2*halfwidth+1)
  RealNum A_coeff = (aFilter[0].afH[aFilter[0].halfwidth] *
                     aFilter[1].afH[aFilter[1].halfwidth] *
                     aFilter[2].afH[aFilter[2].halfwidth]);

  // Create temporary arrays to perform the filter in each direction:
  // (I suppose if I really cared about speed, I could alter Filter1D::Apply() 
  //  function to use pointer arithmatic and allow large strides.
  //  This would also eliminate the need for these temporary arrays.)
  float *aafSource[3];
  float *aafDest[3];
  float *aafMask[3];
  for (int d=0; d < 3; d++) {
    aafDest[d]   = new float [image_size[d]];
    aafSource[d] = new float [image_size[d]];
    aafMask[d]   = NULL;
    if (aaafMask)
      aafMask[d] = new float [image_size[d]];
  }

  // This is a "separable" filter.
  // You can apply the filter sequentially, in the X, Y, Z directions
  // (instead of applying the filter simultaneously in all 3 directions,
  //  which would be much slower)

  // Initially copy aaafSource into aaafDest
  // (We don't want to have to allocate temporary array to 
  //  store the result of each successive filter operation. 
  //  Instead just store the most recent filter operation in aaafDest,
  //  and perform each operation on whatever's currently in aaafDest.)
  for (int iz = 0; iz < image_size[2]; iz++)
    for (int iy = 0; iy < image_size[1]; iy++)
      for (int ix = 0; ix < image_size[0]; ix++)
        aaafDest[iz][iy][ix] = aaafSource[iz][iy][ix];


  int d; //direction where we are applying the filter (x<==>0, y<==>1, z<==>2)


  // Apply the filter in the Z direction (d=2):
  d = 2;
  if (pReportProgress)
    *pReportProgress << "  progress: Applying Z filter. Processing Y plane#"
                     << endl;
  for (int iy = 0; iy < image_size[1]; iy++) {
    if (pReportProgress)
      *pReportProgress << "  " << iy+1 << " / " << image_size[1] << "\n";
    for (int ix = 0; ix < image_size[0]; ix++) {
      // copy the data we need to the temporary array
      for (int iz = 0; iz < image_size[2]; iz++) {
        aafSource[d][iz] = aaafDest[iz][iy][ix];  // copy from aaafDest
        if (aaafMask)
          aafMask[d][iz] = aaafMask[iz][iy][ix];
      }

      // apply the filter to the 1-D temporary array
      aFilter[d].Apply(image_size[d],
                       aafSource[d],
                       aafDest[d],
                       aafMask[d],
                       true);
                       //precompute_mask_times_source);
      for (int iz = 0; iz < image_size[d]; iz++)
        aaafDest[iz][iy][ix] = aafDest[d][iz];  // copy back into aaafDest
    } //for (int ix = 0; ix < image_size[0]; ix++)
  } //for (int iy = 0; iy < image_size[1]; iy++)


  // Apply the filter in the Y direction:
  d=1;
  if (pReportProgress)
    *pReportProgress << "  progress: Applying Y filter. Processing Z plane#"
                     << endl;
  for (int iz = 0; iz < image_size[2]; iz++) {
    if (pReportProgress)
      *pReportProgress << "  " << iz+1 << " / " << image_size[2] << "\n";
    for (int ix = 0; ix < image_size[0]; ix++) {
      // copy the data we need to the temporary array
      for (int iy = 0; iy < image_size[1]; iy++) {
        aafSource[d][iy] = aaafDest[iz][iy][ix]; //data from previous aaafDest
        if (aaafMask)
          aafMask[d][iy] = aaafMask[iz][iy][ix];
      }
      // apply the filter to the 1-D temporary array
      aFilter[d].Apply(image_size[d],
                       aafSource[d],
                       aafDest[d],
                       aafMask[d],
                       true);
                       //precompute_mask_times_source);
      for (int iy = 0; iy < image_size[d]; iy++)
        aaafDest[iz][iy][ix] = aafDest[d][iy];  // copy back into aaafDest
    } //for (int ix = 0; ix < image_size[0]; ix++) {
  } //for (int iz = 0; iz < image_size[2]; iz++)

  // Apply the filter in the X direction:
  d=0;
  if (pReportProgress)
    *pReportProgress << "  progress: Applying X filter. Processing Z plane#"
                     << endl;
  for (int iz = 0; iz < image_size[2]; iz++) {
    if (pReportProgress)
      *pReportProgress << "  " << iz+1 << " / " << image_size[2] << "\n";
    for (int iy = 0; iy < image_size[1]; iy++) {
      // copy the data we need to the temporary array
      for (int ix = 0; ix < image_size[0]; ix++) {
        aafSource[d][ix] = aaafDest[iz][iy][ix]; //data from previous aaafDest
        if (aaafMask)
          aafMask[d][ix] = aaafMask[iz][iy][ix];
      }
      // apply the filter to the 1-D temporary array
      aFilter[d].Apply(image_size[d],
                       aafSource[d],
                       aafDest[d],
                       aafMask[d],
                       true);
                       //precompute_mask_times_source);
      for (int ix = 0; ix < image_size[d]; ix++)
        aaafDest[iz][iy][ix] = aafDest[d][ix];  // copy back into aaafDest
    } //for (int iy = 0; iy < image_size[1]; iy++)
  } //for (int iz = 0; iz < image_size[2]; iz++)

  // delete the temporary arrays
  for (int d=0; d<3; d++) {
    delete [] aafSource[d];
    delete [] aafDest[d];
    if (aafMask[d])
      delete [] aafMask[d];
  }

  if (! normalize) {
    // By default the convolutions used here are all normalized.
    // To counteract that, we can divide by the A_coeff.
  for (int iz = 0; iz < image_size[2]; iz++)
    for (int iy = 0; iy < image_size[1]; iy++)
      for (int ix = 0; ix < image_size[0]; ix++)
        aaafDest[iz][iy][ix] /= A_coeff;
    A_coeff = 1.0;
  }

  return A_coeff;

} //_ApplyGauss3D(aFilter)




template<class RealNum>
RealNum
ApplyGauss3D(int const image_size[3], 
             RealNum ***aaafSource,
             RealNum ***aaafDest,
             RealNum ***aaafMask,
             RealNum const sigma[3],
             int const truncate_halfwidth[3],
             bool normalize = true,
             ostream *pReportProgress = NULL)
{
  assert(aaafSource);
  assert(aaafDest);
  //assert(aaafMask);
  //Allocate filters in all 3 directions.  (Later apply them sequentially.)
  Filter1D<float, int> aFilter[3];
  for (int d=0; d < 3; d++)
    aFilter[d] = GenFilterGauss1D(sigma[d], truncate_halfwidth[d]);

  return _ApplyGauss3D(image_size, 
                       aaafSource,
                       aaafDest,
                       aaafMask,
                       aFilter,
                       normalize,
                       pReportProgress);

} //ApplyGauss3D(sigma, truncate_halfwidth, ...)





template<class RealNum>
RealNum
ApplyGauss3D(int const image_size[3], 
             RealNum ***aaafSource,
             RealNum ***aaafDest,
             RealNum ***aaafMask,
             RealNum const sigma[3],
             RealNum truncate_ratio,
             bool normalize = true,
             ostream *pReportProgress = NULL)
{
  RealNum A;
  int truncate_halfwidth[3];
  if (truncate_ratio > 0)
    for (int d=0; d < 3; d++)
      truncate_halfwidth[d] = floor(sigma[d] * truncate_ratio);

  A = ApplyGauss3D(image_size, 
                   aaafSource,
                   aaafDest,
                   aaafMask,
                   sigma,
                   truncate_halfwidth,
                   normalize,
                   pReportProgress);
  return A;
} // ApplyGauss3D(..., truncate_ratio, ...)






template<class RealNum>
void
ApplyDog3D(int const image_size[3],
           RealNum ***aaafSource,
           RealNum ***aaafDest,
           RealNum ***aaafMask,
           RealNum const sigma_a[3],
           RealNum const sigma_b[3],
           int const truncate_halfwidth[3],//(half the filter width along x,y,z)
           RealNum *pA = NULL,
           RealNum *pB = NULL,
           ostream *pReportProgress = NULL)
{
  RealNum ***aaafTemp; //temporary array to store partially processed tomogram
  RealNum *afTemp;     //temporary array to store partially processed tomogram

  if (pReportProgress)
    *pReportProgress
      << " -- Attempting to allocate space for one more image. --\n"
      << " -- (If this crashes your computer, find a computer  --\n"
      << " --  with more RAM and use \"ulimit\")                 --\n";

  Alloc3D(image_size,
          &afTemp,
          &aaafTemp);

  RealNum A, B;        // let the user know what A B coefficients were used

  // Convolve the original source with the 1st Gaussian
  A = ApplyGauss3D(image_size,
                   aaafSource,
                   aaafDest,   // <-- save result here
                   aaafMask,
                   sigma_a,
                   truncate_halfwidth,
                   true,
                   pReportProgress);

  // Convolve the original source with the 2nd Gaussian
  B = ApplyGauss3D(image_size,
                   aaafSource,
                   aaafTemp,   // <-- save result here
                   aaafMask,
                   sigma_b,
                   truncate_halfwidth,
                   true,
                   pReportProgress);

  // Subtract the second convolved signal from the first
  for (int iz = 0; iz < image_size[2]; iz++)
    for (int iy = 0; iy < image_size[1]; iy++)
      for (int ix = 0; ix < image_size[0]; ix++)
        aaafDest[iz][iy][ix] -= aaafTemp[iz][iy][ix];

  // Deallocate the temporary array
  Dealloc3D(image_size,
            &afTemp,
            &aaafTemp);

} // ApplyDog3D()





template<class RealNum>
void
ApplyDogScaleFree3D(int const image_size[3], //source image size
                    RealNum ***aaafSource,   //source image
                    RealNum ***aaafDest,     //save results here
                    RealNum ***aaafMask,     //ignore voxels where mask==0
                    RealNum const sigma[3],  //Gaussian width in x,y,z drections
                    RealNum delta_sigma_over_sigma, //difference in Gauss widths
                    RealNum truncate_ratio, //how many sigma before truncating?
                    RealNum *pA = NULL,
                    RealNum *pB = NULL,
                    ostream *pReportProgress = NULL)
{
  // "-log" approximates to the "Laplacian of a Gaussian" ("DOG") filter
  // with the difference of two Gaussians ("DOG") filter.
  // (The two Gaussians have widths which are slightly above and
  //  slightly below the width parameters specified by the user.)
  // For background details see:
  // https://en.wikipedia.org/wiki/Blob_detection
  // https://en.wikipedia.org/wiki/Difference_of_Gaussians
  // https://en.wikipedia.org/wiki/Mexican_hat_wavelet
  RealNum sigma_a[3];
  RealNum sigma_b[3];
  sigma_a[0] = sigma[0] * (1.0 - 0.5*delta_sigma_over_sigma);
  sigma_a[1] = sigma[1] * (1.0 - 0.5*delta_sigma_over_sigma);
  sigma_a[2] = sigma[2] * (1.0 - 0.5*delta_sigma_over_sigma);
  sigma_b[0] = sigma[0] * (1.0 + 0.5*delta_sigma_over_sigma);
  sigma_b[1] = sigma[1] * (1.0 + 0.5*delta_sigma_over_sigma);
  sigma_b[2] = sigma[2] * (1.0 + 0.5*delta_sigma_over_sigma);

  int truncate_halfwidth[3];
  for (int d=0; d < 3; d++)
    truncate_halfwidth[d] = floor(truncate_ratio *
                                  MAX(sigma_a[d],
                                      sigma_b[d]));

  ApplyDog3D(image_size,
             aaafSource,
             aaafDest,
             aaafMask,
             sigma_a,
             sigma_b,
             truncate_halfwidth,
             pA,
             pB,
             pReportProgress);

  // Do we care about the overall magnitude of the result?
  //In order to approximate the magnitude of the
  //  "Scale-Normalized Laplacian of a Gaussian"   (a.k.a. "Lnorm")
  //...it is necessary to multiply our result so far by "t / delta_t"
  //
  //https://en.wikipedia.org/wiki/Blob_detection#The_Laplacian_of_Gaussian
  //
  //...where "t" = (1/2) * log_width^2
  //         because filter_mrc includes factor of (1/2) in log_width
  //   and "delta_t" = t * delta_sigma_over_sigma^2
  //   (Note: The "(1/2)" is present in both "t" and "delta_t" so the ratio 
  //          between them is unaffected by its presence)

  RealNum t_over_delta_t = 1.0 / SQR(delta_sigma_over_sigma);

  for (int iz = 0; iz < image_size[2]; iz++)
    for (int iy = 0; iy < image_size[1]; iy++)
      for (int ix = 0; ix < image_size[0]; ix++)
        aaafDest[iz][iy][ix] *= t_over_delta_t;

  if (pA)
    *pA *= t_over_delta_t;
  if (pB)
    *pB *= t_over_delta_t;

  // Optional: The laplacian is the 2nd derivative w respect to voxel location
  //           Should we instead make it the 2nd derivative w respect to
  //           physical distance?  (This makes it easier to compare
  //           images taken at different resolutions.)
  //           If so, multiply A, B, and aaafDest result by 1/voxel_width^2
  //           (Don't worry about this for now.  Let the calller deal with this)
  //            

} // ApplyDogScaleFree3D()





template<class RealNum>
void
ApplyDogScaleFree3D(int const image_size[3], //source image size
                    RealNum ***aaafSource,   //source image
                    RealNum ***aaafDest,     //save results here
                    RealNum ***aaafMask,     //ignore voxels where mask==0
                    RealNum sigma,  //Gaussian width in x,y,z drections
                    RealNum delta_sigma_over_sigma, //difference in Gauss widths
                    RealNum truncate_ratio,  //how many sigma before truncating?
                    RealNum *pA = NULL,
                    RealNum *pB = NULL,
                    ostream *pReportProgress = NULL)
{
  RealNum sigma_xyz[3] = {sigma, sigma, sigma};
  ApplyDogScaleFree3D(image_size,
                      aaafSource,
                      aaafDest,
                      aaafMask,
                      sigma_xyz,
                      delta_sigma_over_sigma,
                      truncate_ratio,
                      pA,
                      pB,
                      pReportProgress);
}



// Find all scale-invariant blobs in the image:
template<class RealNum>
void
BlobDog3D(int const image_size[3], //source image size
          RealNum ***aaafSource,   //source image
          RealNum ***aaafMask,     //ignore voxels where mask==0
          const vector<RealNum>& blob_widths, // blob widths to try, ordered
          vector<array<int,3> >& minima_crds, // store minima x,y,z coords here
          vector<array<int,3> >& maxima_crds, // store maxima x,y,z coords here
          vector<RealNum>& minima_sigma, // corresponding width for that minima
          vector<RealNum>& maxima_sigma, // corresponding width for that maxima
          vector<RealNum>& minima_scores, // what was the blob's score?
          vector<RealNum>& maxima_scores, // (score = intensity after filtering)
          RealNum delta_sigma_over_sigma,// param for approximating LOG with DOG
          RealNum minima_threshold,    // discard blobs with unremarkable scores
          RealNum maxima_threshold,    // discard blobs with unremarkable scores
          RealNum truncate_ratio,      // how many sigma before truncating?
          // optional arguments
          ostream *pReportProgress = NULL, // report progress to the user?
          RealNum ****aaaafI = NULL, //preallocated memory for filtered images
          RealNum **aafI = NULL)     //preallocated memory for filtered images
                                     //(sometimes the caller wants to access)
{

  // We need 3 images to store the result of filtering the image
  // using DOG filters with 3 different widths.  Store those images here:

  cerr << " -- Attempting to allocate space for 3 more images.  --\n"
       << " -- (If this crashes your computer, find a computer  --\n"
       << " --  with more RAM and use \"ulimit\")               --\n";

  bool preallocated = ! (aaaafI == NULL);

  if (! preallocated) {
    assert(aaaafI == NULL);
    assert(aafI == NULL);
    aaaafI = new RealNum*** [3];
    aafI = new RealNum* [3];
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

  for (int ir = 0; ir < blob_widths.size(); ir++) {

    cerr << "---------- sigma[" << ir << "] = "
         << blob_widths[ir]
         << " (in voxels)----------\n";

    // Run the image through a DOG filter with the currently selected width.
    // Compare the resulting filtered image with filtered images
    // smaller blob widths.
    // If a given voxel has a value which is larger than it's surrounding 26
    // voxels in this image, -as well as- the same 27 voxels which were filtered
    // using larger and smaller blob widths, respectively...
    //   ... -then- it is a local maximum in 4-dimensional X,Y,Z,R space.
    // Keep track of all of these X,Y,Z,R local maxima voxels and their values
    // (and the local minima as well).

    cerr << "--- Applying DOG filter using sigma[" << ir << "] = "
         << blob_widths[ir] << " (in voxels) ---\n";

    // We are going to apply a DOG filter to the original image and
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

    //Apply the DOG filter using the most recent blob width:
    ApplyDogScaleFree3D(image_size,
                        aaafSource,
                        aaaafI[j2i[1]], //<-store the most recent filtered image
                        aaafMask,
                        blob_widths[ir],
                        delta_sigma_over_sigma,
                        truncate_ratio);

    // We must have filtered at least 3 images using different blob widths
    if (ir < 2)
      continue;

    // The blob widths should be ordered consistently (increasing or decreasing)
    assert(((blob_widths[ir-2] < blob_widths[ir-1]) &&   //increasing order
            (blob_widths[ir-1] < blob_widths[ir])) ||    
           ((blob_widths[ir-2] > blob_widths[ir-1]) &&   //decreasing order
            (blob_widths[ir-1] > blob_widths[ir])));


    cerr << "--- searching for local minima & maxima at sigma[" << ir-1 << "] = "
         << blob_widths[ir-1] << " ---\n";

    for (int iz = 0; iz < image_size[2]; iz++) {
      cerr << "  " << iz+1 << " / " << image_size[2] << "\n";
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
                  if (aaaafI[ j2i[jr] ][ Iz ][ Iy ][ Ix ] <=
                      aaaafI[ j2i[0]  ][ iz ][ iy ][ ix ])
                    is_minima = false;
                  if (aaaafI[ j2i[jr] ][ Iz ][ Iy ][ Ix ] >=
                      aaaafI[ j2i[0]  ][ iz ][ iy ][ ix ])
                    is_maxima = false;
                }
              }
            }
          }

          RealNum score = aaaafI[ j2i[0] ][ iz ][ iy ][ ix ];

          if ((! aaafMask) || (aaafMask[iz][iy][ix] != 0))
          {
            if (is_minima &&
                ((score < minima_threshold) ||
                 (minima_threshold > maxima_threshold)))
            {
              array<int, 3> ixiyiz;
              ixiyiz[0] = ix;
              ixiyiz[1] = iy;
              ixiyiz[2] = iz;
              minima_crds.push_back(ixiyiz);
              minima_sigma.push_back(blob_widths[ir-1]);
              minima_scores.push_back(score);
            }
            if (is_maxima &&
                ((score > maxima_threshold) ||
                 (minima_threshold > maxima_threshold)))
            {
              array<int, 3> ixiyiz;
              ixiyiz[0] = ix;
              ixiyiz[1] = iy;
              ixiyiz[2] = iz;
              maxima_crds.push_back(ixiyiz);
              maxima_sigma.push_back(blob_widths[ir-1]);
              maxima_scores.push_back(score);
            }
          }
          assert(! (is_minima && is_maxima));
        } //for (int ix=0; ix<image_size[0]; ix++) {
      } //for (int iy=0; iy<image_size[1]; iy++) {
    } //for (int iz=0; iz<image_size[2]; iz++) {

    cerr << "--- found " << minima_crds.size ()
         << " and " << maxima_crds.size()
         << " local minima and maxima, respectively ---\n" << endl;

    assert((minima_crds.size() == minima_sigma.size()) &&
           (minima_crds.size() == minima_scores.size()));
    assert((maxima_crds.size() == maxima_sigma.size()) &&
           (maxima_crds.size() == maxima_scores.size()));

  } //for (ir = 0; ir < blob_widths.size(); ir++)


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

} //BlobDog3D()


//apply permutation in-place, credit: Raymond Chen
template<typename T>
void
apply_permutation(vector<T>& v,
                  vector<int>& indices)
{
  using swap;
  for (size_t i = 0; i < indices.size(); i++) {
    auto current = i;
    while (i != indices[current]) {
      auto next = indices[current];
      swap(v[current], v[next]);
      indices[current] = current;
      current = next;
    }
    indices[current] = current;
  }
} //apply_permutation()


          
          
template<class RealNum>
void
DiscardOverlappingBlobs3D(vector<array<int,3> >& blob_crds,
                          vector<RealNum>& blob_sigma, 
                          vector<RealNum>& blob_scores,
                          occupancy_table_size) {
  cerr << "     -- Attempting to allocate space for one more image map --\n"
       << "     -- (If this crashes your computer, find a computer     --\n"
       << "     --  with more RAM and use \"ulimit\")                    --\n";

  bool *abOccupied;
  bool *aaabOccupied;
  Alloc3D(occupancy_table_size,
          &abOccupied,
          &aaabOccupied);
  for(int iz = 0; iz < occupancy_table_size[2]; iz++)
    for(int iy = 0; iy < occupancy_table_size[1]; iy++)
      for(int ix = 0; ix < occupancy_table_size[0]; ix++)
        aaabOccupied[iz][iy][ix] = false;

  vector<long long> i_keep_blob;
  for (long long i=0; i < nmin; i++) {
    bool discard = false;
    float reff = (1.0-nonmax_suppression_max_overlap)*blob_sigma[i]*sqrt(2.0);
    int Reff = ceil(reff);
    int Reffsq = ceil(reff*reff);
    ix = blob_crds[i][0];
    iy = blob_crds[i][1];
    iz = blob_crds[i][2];
    for(int jz = -Reff; jz <= Reff && (! discard); jz++) {
      for(int jy = -Reff; jy <= Reff && (! discard); jy++) {
        for(int jx = -Reff; jx <= Reff && (! discard); jx++) {
          rsq = jx*jx + jy*jy + jz*jz;
          if (rsq > Reffsq)
            continue;
          if (aaabOccupied[iz+jz][iy+jy][ix+jx])
            discard = true;
          aaabOccupied[iz+jz][iy+jy][ix+jx] = true;
        }
      }
    }
    if (discard) {


    }
  }

  Dealloc3D(occupancy_table_size,
            &abOccupied,
            &aaabOccupied);
} //DiscardOverlappingBlobs3D()




// Find all scale-invariant blobs in the image,
// This version can discard blobs which overlap with existing blobs
// using the "max_overlap" parameter (which varies between 0 and 1).
// (This is sometimes called "non-max suppression")
// If two blobs lie within a distance of
//     max_overlap * (R1 + R2)
// from each other, then the blob with the lower score is discarded.
template<class RealNum>
void
BlobDog3D(int const image_size[3], //source image size
          RealNum ***aaafSource,   //source image
          RealNum ***aaafMask,     //ignore voxels where mask==0
          const vector<RealNum>& blob_widths, // blob widths to try, ordered
          vector<array<int,3> >& minima_crds, // store minima x,y,z coords here
          vector<array<int,3> >& maxima_crds, // store maxima x,y,z coords here
          vector<RealNum>& minima_sigma, // corresponding width for that minima
          vector<RealNum>& maxima_sigma, // corresponding width for that maxima
          vector<RealNum>& minima_scores, // what was the blob's score?
          vector<RealNum>& maxima_scores, // (score = intensity after filtering)
          RealNum delta_sigma_over_sigma,// param for approximating LOG with DOG
          RealNum minima_threshold,    // discard blobs with unremarkable scores
          RealNum maxima_threshold,    // discard blobs with unremarkable scores
          RealNum truncate_ratio,      // how many sigma before truncating?
          RealNum max_overlap,
          // optional arguments
          ostream *pReportProgress = NULL, // report progress to the user?
          RealNum ****aaaafI = NULL, //preallocated memory for filtered images
          RealNum **aafI = NULL)     //preallocated memory for filtered images
                                     //(sometimes the caller wants to access)
{

  BlobDog3D(image_size,
            aaafSource,   //source image
            aaafMask,     //ignore voxels where mask==0
            blob_widths, // blob widths to try, ordered
            minima_crds, // store minima x,y,z coords here
            maxima_crds, // store maxima x,y,z coords here
            minima_sigma, // corresponding width for that minima
            maxima_sigma, // corresponding width for that maxima
            minima_scores, // what was the blob's score?
            maxima_scores, // (score = intensity after filtering)
            delta_sigma_over_sigma,// param for approximating LOG with DOG
            minima_threshold,    // discard blobs with unremarkable scores
            maxima_threshold,    // discard blobs with unremarkable scores
            truncate_ratio,      // how many sigma before truncating?
            // optional arguments
            pReportProgress, // report progress to the user?
            aaaafI, //preallocated memory for filtered images
            aafI);     //preallocated memory for filtered images

  cerr << "--- Removing overlapping blobs ---" << endl;

  long long n_min = minima_crds.size();
  assert(n_min == minima_sigma.size());
  assert(n_min == minima_scores.size());
  long long n_max = maxima_crds.size();
  assert(n_max == maxima_sigma.size());
  assert(n_max == maxima_scores.size());
  vector<tuple<float, long long> > score_index_minima(n_min);
  vector<tuple<float, long long> > score_index_maxima(n_max);
  for (long long i = 0; i < n_min; i++)
    score_index_minima[i] = tuple<float, long>(minima_scores[i], i);
  for (long long i = 0; i < n_max; i++)
    score_index_maxima[i] = tuple<float, long>(maxima_scores[i], i);

  if (n_min > 0) {
    cerr << "-- Sorting minima blobs according to their scores in descending order...";
    sort(score_index_minima.rbegin(),
         score_index_minima.rend(),
         lessthan());
    vector<long long> permute_min;
    for (i = 0; i < score_index_minima.size(); i++)
      permute_min[i] = score_index_minima[i];
    score_index_minima.clear();
    apply_permutation(minima_crds,  permute_min);
    apply_permutation(minima_sigma, permute_min);
    apply_permutation(minima_scores, permute_min);
    cerr << "done --" << endl;
  }
  if (n_max > 0) {
    cerr << "-- Sorting maxima blobs according to their scores in ascending order...";
    sort(score_index_maxima.begin(),
         score_index_maxima.end(),
         lessthan());
    vector<long long> permute_max;
    for (i = 0; i < score_index_maxima.size(); i++)
      permute_max[i] = score_index_maxima[i];
    score_index_maxima.clear();
    apply_permutation(maxima_crds,  permute_max);
    apply_permutation(maxima_sigma, permute_max);
    apply_permutation(maxima_scores, permute_max);
    cerr << "done --" << endl;
  }

  
  if (max_overlap < 1.0) {
    cerr << "-- Discarding overlapping minima blobs (volume overlap > "
         << max_overlap << ")...";

    DiscardOverlappingBlobs3D(minima_crds,
                              minima_sigma,
                              minima_scores,
                              image_size);

    cerr << "done --" << endl;
    cerr << "-- Discarding overlapping maxima blobs (volume overlap > "
         << max_overlap << ")...";
    DiscardOverlappingBlobs3D(maxima_crds,
                              maxima_sigma,
                              maxima_scores,
                              image_size);
    cerr << "done --" << endl;
  }

} //BlobDog3D()




#endif //#ifndef _FILTER3D_H
