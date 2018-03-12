#ifndef _FILTER2D_H
#define _FILTER2D_H

#include <cstring>
#include <ostream>
using namespace std;





template<class RealNum, class Integer>

class Filter2D {
public:
  RealNum *af;
  RealNum **aafH;
  Integer halfwidth[2]; //num pixels from filter center to edge in x,y directions
  Integer array_size[2]; //size of the array in x,y directions (in voxels)



  // Apply the filter to tomographic data in the "aafSource" array
  // Save the results in the "aafDest" array.  (A "mask" is optional.)
  // All arrays are 2D and assumed to be the same size
  void Apply(Integer const size_source[2],
             RealNum **aafSource,
             RealNum **aafDest,
             RealNum **aafMask = NULL,
             bool normalize = false) const
             //bool precompute_mask_times_source = true) const
  {

    // Apply the filter to the original tomogram data. (Store in aafSource)
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

    //if (aafMask && precompute_mask_times_source)
    if (aafMask) 
      for (int iy=0; iy<size_source[1]; iy++)
        for (int ix=0; ix<size_source[0]; ix++)
          aafSource[iy][ix] *= aafMask[iy][ix];

    for (Integer iy=0; iy<size_source[1]; iy++) {

      for (Integer ix=0; ix<size_source[0]; ix++) {

        if ((aafMask) && (aafMask[iy][ix] == 0.0)) {
          aafDest[iy][ix] = 0.0;
          continue;
        }
          
        RealNum g = 0.0;
        RealNum denominator = 0.0;

        for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {

          if ((iy-jy < 0) || (size_source[1] <= iy-jy))
            continue;

          for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {

            if ((ix-jx < 0) || (size_source[0] <= ix-jx))
              continue;

            if ((jx == 0) && (jy == 0)) {
              g += 1.0;
              g -= 1.0;
            }
            RealNum delta_g = 
              aafH[jy+halfwidth[1]][jx+halfwidth[0]] * aafSource[iy-jy][ix-jx];

            //if (! precompute_mask_times_source)
            //  delta_g *= aafMask[iy-jy][ix-jx];

            g += delta_g;
              // Note: We previously applied the mask by multiplying
              //          aaDensity[iy][ix]
              //         by aafMask[iy][ix]   (if present)
            if (normalize) {
              if (aafMask)
                denominator +=
                  aafMask[iy-jy][ix-jx] * aafH[jy+halfwidth[1]][jx+halfwidth[0]];
                                               
              else
                denominator += aafH[jy+halfwidth[1]][jx+halfwidth[0]];
            }
                                          
            // Note: If there were no mask, and if the filter is normalized
            // then denominator=1 always, and we could skip the line above.
          }
        }

        if (normalize) {
          if (denominator != 0.0)
            g /= denominator;
          else
            //Otherwise, this position lies outside the mask region.
            g = 0.0;
        }

        aafDest[iy][ix] = g;
      }
    }
  } //Apply()


  void Alloc(Integer const set_halfwidth[2]) {
    //Integer array_size[2];
    for(Integer d=0; d < 2; d++) {
      halfwidth[d] = set_halfwidth[d];
      array_size[d] = 1 + 2*halfwidth[d];
    }
    Alloc2D(array_size, &af, &aafH);
    for (int iy = 0; iy < array_size[1]; iy++)
      for (int ix = 0; ix < array_size[0]; ix++)
        aafH[iy][ix] = -1.0e38; //(if uninitiliazed memory read, we will know)
  }

  void Dealloc() {
    //Integer array_size[2];
    for(Integer d=0; d < 2; d++) {
      array_size[d] = 1 + 2*halfwidth[d];
      halfwidth[d] = -1;
    }
    Dealloc2D(array_size, &af, &aafH);
    for(Integer d=0; d < 2; d++)
      array_size[d] = -1;
  }


  void Resize(Integer const set_halfwidth[2]) {
    //Integer array_size[2];
    if (af && aafH) {
      for(Integer d=0; d < 2; d++)
        array_size[d] = 1 + 2*halfwidth[d];
      Dealloc2D(array_size, &af, &aafH);
    }
    for(Integer d=0; d < 2; d++) {
      halfwidth[d] = set_halfwidth[d];
      array_size[d] = 1 + 2*halfwidth[d];
    }
    Alloc2D(array_size, &af, &aafH);
  }

  inline Filter2D<RealNum, Integer>&
    operator = (const Filter2D<RealNum, Integer>& source) {
    Resize(source.halfwidth); // allocates and initializes af and aaafH
    //for(Int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
    //  for(Int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
    //    aafH[iy][ix] = source.aafH[iy][ix];
    // Use memcpy() instead:
    memcpy(af,
           source.af,
           ((1+2*halfwidth[0]) * (1+2*halfwidth[1]))
           *sizeof(RealNum));
  } // operator = ()

  void Init() {
    halfwidth[0] = -1;
    halfwidth[1] = -1;
    af = NULL;
    aafH = NULL;
  }

  Filter2D(Integer const set_halfwidth[2]) {
    Init();
    Resize(set_halfwidth);
  }

  Filter2D() {
    Init();
  }

  Filter2D(const Filter2D<RealNum, Integer>& f) {
    Init();
    *this = f;
  }

  ~Filter2D() {
    Dealloc();
  }

  void Normalize() {
    // Make sure the sum of the filter weights is 1
    RealNum total = 0.0;
    for (Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
      for (Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
        total += aafH[iy+halfwidth[1]][ix+halfwidth[0]];
    for (Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
      for (Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
        aafH[iy+halfwidth[1]][ix+halfwidth[0]] /= total;
  }

}; // class Filter2D





template<class RealNum>
// Create a 2-D filter and fill it with a "generalized Gaussian" function:
//    h_xy(r) = A*exp(-r^m)
// where   r  = sqrt((x/s_x)^2 + (y/s_y)^2)
//   and   A  is determined by normalization of the discrete sum
// Note: "A" is equal to the value stored in the middle of the array,
//       The caller can determine what "A" is by looking at this value.
Filter2D<RealNum, int>
GenFilterGenGauss2D(RealNum width[2],    //"s_x", "s_y" parameters
                    RealNum m_exp,       //"m" exponent parameter
                    int halfwidth[2],
                    RealNum *pA=NULL,    //optional:report A coeff to user
                    ostream *pReportProgress = NULL)
{
  RealNum window_threshold = 1.0;
  for (int d=0; d<2; d++) {
    RealNum h = ((width[d]>0) ? exp(-pow(halfwidth[d]/width[d], m_exp)) : 1.0);
    if (h < window_threshold)
      window_threshold = h;
  }
  Filter2D<RealNum, int> filter(halfwidth);
  RealNum total = 0;
  for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
    for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
      RealNum r = sqrt(SQR(ix/width[0]) + SQR(iy/width[1]));
      RealNum h = ((r>0) ? exp(-pow(r, m_exp)) : 1.0);
      if (ABS(h) < window_threshold)
        h = 0.0; // this eliminates corner entries which fall below threshold
                 // (and eliminates anisotropic artifacts due to these corners)
                 // There's no reason to keep any entries less than min value.
      filter.aafH[iy+halfwidth[1]]
                 [ix+halfwidth[0]] = h;
      total += h;
    }
  }
  // normalize:
  for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
    for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
      filter.aafH[iy+halfwidth[1]]
                 [ix+halfwidth[0]] /= total;
      //FOR DEBUGGING REMOVE EVENTUALLY
      if (pReportProgress)
        *pReportProgress << "GenGauss2D:" //<< window_threshold
                         <<" aafH["<<iy<<"]["<<ix<<"] = "
                         << filter.aafH[iy+halfwidth[1]]
                                       [ix+halfwidth[0]] << endl;
    }
  }
  return filter;
} //GenFilterGenGauss2D(width, m_exp, halfwidth)






template<class RealNum>
Filter2D<RealNum, int>
GenFilterGenGauss2D(RealNum width[2],            //"s_x", "s_y" parameters
                    RealNum m_exp,               //"m" parameter in formula
                    RealNum filter_cutoff_ratio,
                    RealNum *pA=NULL,    //optional:report A coeff to user
                    ostream *pReportProgress = NULL)
{
  // choose the width of the filter window based on the filter_cutoff_ratio
  int halfwidth[2];
  int ix = 0;

  for (int d=0; d<2; d++) {
    halfwidth[d] = floor(width[d]*filter_cutoff_ratio);
  }
  return GenFilterGenGauss2D(width,
                             m_exp,
                             halfwidth,
                             pA,
                             pReportProgress);
} //GenFilterGenGauss2D(width, m_exp, filter_cutoff_ratio)





template<class RealNum>
// Create a 2-D filter and fill it with a difference of (generalized) Gaussians:
// This version requires that the caller has already created individual
// filters for the two gaussians.
// All this function does is subtract one filter from the other (and rescale).
Filter2D<RealNum, int> 
_GenFilterDogg2D(RealNum width_a[2],  //"a" parameter in formula
                 RealNum width_b[2],  //"b" parameter in formula
                 RealNum m_exp,  //"m" parameter in formula
                 RealNum n_exp,  //"n" parameter in formula
                 Filter2D<RealNum, int>& filterXY_A, //filters for the two
                 Filter2D<RealNum, int>& filterXY_B, //gaussians
                 RealNum *pA=NULL, //optional:report A,B coeffs to user
                 RealNum *pB=NULL, //optional:report A,B coeffs to user
                 ostream *pReportProgress = NULL)
{
  RealNum A, B;
  //A, B = height of the central peak
  //       The central peak is located in the middle of the filter's array
  //       (at position "halfwidth")
  A = filterXY_A.aafH[filterXY_A.halfwidth[0]]
                     [filterXY_A.halfwidth[1]];
  B = filterXY_B.aafH[filterXY_B.halfwidth[0]]
                     [filterXY_B.halfwidth[1]];


  // The "difference of gaussians" filter is the difference between
  // these two (generalized) gaussian filters.
  int halfwidth[2];
  halfwidth[0] = MAX(filterXY_A.halfwidth[0], filterXY_B.halfwidth[0]);
  halfwidth[1] = MAX(filterXY_A.halfwidth[1], filterXY_B.halfwidth[1]);
  Filter2D<RealNum, int> filter(halfwidth);

  if (pReportProgress)
    *pReportProgress  << "Array of 2D filter entries:" << endl;

  for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
    for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
      filter.aafH[iy+halfwidth[1]]
                 [ix+halfwidth[0]] = 0.0;

      // The two filters may have different widths, so we have to check
      // that ix and iy lie within the domain of these two filters before
      // adding or subtracting their values from the final GDOG filter.
      if (((-filterXY_A.halfwidth[0]<=ix) && (ix<=filterXY_A.halfwidth[0])) &&
          ((-filterXY_A.halfwidth[1]<=iy) && (iy<=filterXY_A.halfwidth[1])))

        filter.aafH[iy+halfwidth[1]]
                   [ix+halfwidth[0]] +=
          filterXY_A.aafH[iy+filterXY_A.halfwidth[1]]
                         [ix+filterXY_A.halfwidth[0]];   //  /  (A-B); COMMENTING OUT
      // COMMENTING OUT: (The factor of 1/(A-B) insures that the central peak has height 1)

      if (((-filterXY_B.halfwidth[0]<=ix) && (ix<=filterXY_B.halfwidth[0])) &&
          ((-filterXY_B.halfwidth[1]<=iy) && (iy<=filterXY_B.halfwidth[1])))

        filter.aafH[iy+halfwidth[1]]
                   [ix+halfwidth[0]] -=
          filterXY_B.aafH[iy+filterXY_B.halfwidth[1]]
                         [ix+filterXY_B.halfwidth[0]];   //  /  (A-B); COMMENTING OUT


      //FOR DEBUGGING REMOVE EVENTUALLY
      if (pReportProgress)
        *pReportProgress << "GenDogg2D: aafH["<<iy<<"]["<<ix<<"] = "
                         << filter.aafH[iy+halfwidth[1]]
                                       [ix+halfwidth[0]] << endl;

      // *pReportProgress  << aafH[iy+halfwidth[1]][ix+halfwidth[0]];
      //if (ix == halfwidth[0])
      //  *pReportProgress << "\n";
      //else
      //  *pReportProgress << " ";

    } // for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
  } // for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {

  if (pReportProgress)
    *pReportProgress << "\n";

  // COMMENTING OUT 1/(A-B)
  //if (pA && pB) {
  //  *pA = A/(A-B); // Rescale A and B numbers returned to the caller
  //  *pB = B/(A-B); // (because we divided the array entries by (A-B) earlier)
  //}

  return filter;

} //_GenFilterDogg2D()



template<class RealNum>
// Create a 2-D filter and fill it with a difference of (generalized) Gaussians:
Filter2D<RealNum, int> 
GenFilterDogg2D(RealNum width_a[2],  //"a" parameter in formula
                RealNum width_b[2],  //"b" parameter in formula
                RealNum m_exp,       //"m" parameter in formula
                RealNum n_exp,       //"n" parameter in formula
                int halfwidth[2],
                RealNum *pA = NULL,  //optional:report A,B coeffs to user
                RealNum *pB = NULL,  //optional:report A,B coeffs to user
                ostream *pReportProgress = NULL)
{
  Filter2D<RealNum, int> filterXY_A =
    GenFilterGenGauss2D(width_a,      //"a_x", "a_y" gaussian width parameters
                        m_exp,        //"n" exponent parameter
                        halfwidth);
                        //pReportProgress);

  Filter2D<RealNum, int> filterXY_B =
    GenFilterGenGauss2D(width_b,      //"b_x", "b_y" gaussian width parameters
                        n_exp,        //"n" exponent parameter
                        halfwidth);
                        //pReportProgress);

  return _GenFilterDogg2D(width_a,
                          width_b,
                          m_exp,
                          n_exp,
                          filterXY_A, filterXY_B,
                          pA,
                          pB,
                          pReportProgress);
} //GenFilterDogg2D(...halfwidth...)



#endif //#ifndef _FILTER2D_H
