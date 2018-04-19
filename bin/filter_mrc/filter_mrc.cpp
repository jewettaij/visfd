// (Note: For gcc version 4.8.3, you must compile using: g++ -std=c++11)


#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <array>
//#include <fftw3.h>  not needed yet
using namespace std;
#include <err_report.h>
#include <alloc2d.h>
#include <alloc3d.h>
#include <filter1d.h>
#include <filter2d.h>
#include <filter3d.h>
#include <threshold.h>
#include <mrc_simple.h>
#include <random_gen.h>
#include "settings.h"





template<class RealNum>
Filter2D<RealNum, int>
GenFilterGenGauss2D(RealNum width[2],    //"s_x", "s_y" parameters
                    RealNum m_exp,       //"m" parameter in formula
                    RealNum filter_truncate_ratio, //controls window width
                    RealNum filter_truncate_threshold, //controls window width
                    RealNum *pA=NULL,    //optional:report A coeff to user
                    ostream *pReportEquation = NULL)
{
  // choose the filter window width based on the filter_truncate_threshold
  int halfwidth[2];
  int ix = 0;

  if (filter_truncate_ratio < 0.0)
    // Choose the filter domain window based on the "window_threshold"
    //    filter_truncate_threshold = exp(-(filter_truncate_ratio)^m_exp);
    //    -> filter__truncate_ratio^m_exp = -log(filter_window_threshold)
    filter_truncate_ratio = pow(-log(filter_truncate_threshold), 1.0/m_exp);

  return GenFilterGenGauss2D(width,
                             m_exp,
                             filter_truncate_ratio,
                             pA,
                             pReportEquation);

} //GenFilterGenGauss2D(...,filter_truncate_ratio, filter_truncate_thresh,...)





template<class RealNum>
Filter3D<RealNum, int>
GenFilterGenGauss3D(RealNum width[3],    //"s_x", "s_y", "s_z", parameters
                    RealNum m_exp,       //"m" parameter in formula
                    RealNum filter_truncate_ratio, //controls window width
                    RealNum filter_truncate_threshold, //controls window width
                    RealNum *pA=NULL,    //optional:report A coeff to user
                    ostream *pReportEquation = NULL)
{
  // choose the filter window width based on the filter_truncate_threshold
  int halfwidth[3];
  int ix = 0;

  if (filter_truncate_ratio < 0.0)
    // Choose the filter domain window based on the "window_threshold"
    //    filter_truncate_threshold = exp(-(filter_truncate_ratio)^m_exp);
    //    -> filter__truncate_ratio^m_exp = -log(filter_window_threshold)
    filter_truncate_ratio = pow(-log(filter_truncate_threshold), 1.0/m_exp);

  return GenFilterGenGauss3D(width,
                             m_exp,
                             filter_truncate_ratio,
                             pA,
                             pReportEquation);
} //GenFilterGenGauss2D(...,filter_truncate_ratio, filter_truncate_thresh,...)





template<class RealNum>
// Create a 2-D filter and fill it with a difference of (generalized) Gaussians:
Filter2D<RealNum, int> 
GenFilterDogg2D(RealNum width_a[2],  //"a" parameter in formula
                RealNum width_b[2],  //"b" parameter in formula
                RealNum m_exp,       //"m" parameter in formula
                RealNum n_exp,       //"n" parameter in formula
                RealNum filter_truncate_ratio,//how many sigma before cutoff?
                RealNum filter_truncate_threshold,//cutoff when decay below threshold
                RealNum *pA=NULL, //optional:report A,B coeffs to user
                RealNum *pB=NULL, //optional:report A,B coeffs to user
                ostream *pReportProgress = NULL)
{
  RealNum filter_truncate_ratio_A = filter_truncate_ratio;
  RealNum filter_truncate_ratio_B = filter_truncate_ratio;
  if (filter_truncate_ratio < 0.0) {
    // Choose the filter domain window based on the "filter_truncate_threshold"
    //    filter_truncate_threshold = exp(-filter_truncate_ratio^m_exp);
    //    -> filter_truncate_ratio^m_exp = -log(filter_truncate_threshold)
    filter_truncate_ratio_A = pow(-log(filter_truncate_threshold), 1.0/m_exp);
    // Do the same thing for the other Gaussian:
    filter_truncate_ratio_B = pow(-log(filter_truncate_threshold), 1.0/n_exp);
  }
  Filter2D<RealNum, int> filterXY_A =
    GenFilterGenGauss2D(width_a, //"a_x", "a_y" gaussian width parameters
                        m_exp,   //"m" exponent parameter
                        filter_truncate_ratio_A);

  Filter2D<RealNum, int> filterXY_B =
    GenFilterGenGauss2D(width_b, //"b_x", "b_y" gaussian width parameters
                        n_exp,   //"n" exponent parameter
                        filter_truncate_ratio_B);

  return _GenFilterDogg2D(width_a,
                          width_b,
                          m_exp,
                          n_exp,
                          filterXY_A, filterXY_B,
                          pA,
                          pB,
                          pReportProgress);
} //GenFilterDogg2D(...,filter_truncate_ratio, filter_truncate_threshold,...)





template<class RealNum>
// Create a 3-D filter and fill it with a difference of (generalized) Gaussians:
Filter3D<RealNum, int> 
GenFilterDogg3D(RealNum width_a[3],  //"a" parameter in formula
                RealNum width_b[3],  //"b" parameter in formula
                RealNum m_exp,       //"m" parameter in formula
                RealNum n_exp,       //"n" parameter in formula
                RealNum filter_truncate_ratio,//how many sigma before cutoff?
                RealNum filter_truncate_threshold,//cutoff when decay below threshold
                RealNum *pA=NULL, //optional:report A,B coeffs to user
                RealNum *pB=NULL, //optional:report A,B coeffs to user
                ostream *pReportProgress = NULL)
{
  RealNum filter_truncate_ratio_A = filter_truncate_ratio;
  RealNum filter_truncate_ratio_B = filter_truncate_ratio;
  if (filter_truncate_ratio < 0.0) {
    // Choose the filter domain window based on the "filter_truncate_threshold"
    //    filter_truncate_threshold = exp(-filter_truncate_ratio^m_exp);
    //    -> filter_truncate_ratio^m_exp = -log(filter_truncate_threshold)
    filter_truncate_ratio_A = pow(-log(filter_truncate_threshold), 1.0/m_exp);
    // Do the same thing for the other Gaussian:
    filter_truncate_ratio_B = pow(-log(filter_truncate_threshold), 1.0/n_exp);
  }
  Filter3D<RealNum, int> filter_A =
    GenFilterGenGauss3D(width_a, //"a_x", "a_y", "a_z" gaussian width parameters
                        m_exp,   //"m" exponent parameter
                        filter_truncate_ratio_A);

  Filter3D<RealNum, int> filter_B =
    GenFilterGenGauss3D(width_b, //"b_x", "b_y", "b_z" gaussian width parameters
                        n_exp,   //"n" exponent parameter
                        filter_truncate_ratio_B);

  return _GenFilterDogg3D(width_a,
                          width_b,
                          m_exp,
                          n_exp,
                          filter_A, filter_B,
                          pA,
                          pB,
                          pReportProgress);

} //GenFilterDogg3D(...,filter_truncate_ratio, filter_truncate_threshold,...)










template<class RealNum>
RealNum
ApplyGauss3D(int const image_size[3], 
             RealNum ***aaafSource,
             RealNum ***aaafDest,
             RealNum ***aaafMask,
             RealNum const sigma[3],
             RealNum filter_truncate_ratio,
             RealNum filter_truncate_threshold,
             bool normalize = true,
             ostream *pReportProgress = NULL)
{

  if (filter_truncate_ratio <= 0) {
    assert(filter_truncate_threshold > 0.0);
    //    filter_truncate_threshold = exp(-(1/2)*filter_truncate_ratio^2);
    //    -> filter_truncate_ratio^2 = -2*log(filter_truncate_threshold)
    filter_truncate_ratio = sqrt(-2*log(filter_truncate_threshold));
  }
  return ApplyGauss3D(image_size, 
                      aaafSource,
                      aaafDest,
                      aaafMask,
                      sigma,
                      filter_truncate_ratio,
                      normalize,
                      pReportProgress);

} // ApplyGauss3D(..., filter_truncate_ratio, filter_truncate_threshold, ...)






template<class RealNum>
void
ApplyDog3D(int const image_size[3], //source image size
           RealNum ***aaafSource,   //source image
           RealNum ***aaafDest,     //save results here
           RealNum ***aaafMask,     //ignore voxels where mask==0
           RealNum sigma_a[3],
           RealNum sigma_b[3],
           RealNum filter_truncate_ratio,
           RealNum filter_truncate_threshold,
           RealNum *pA = NULL,
           RealNum *pB = NULL)
{
  RealNum ***aaafTemp; //temporary array to store partially processed tomogram
  RealNum *afTemp;     //temporary array to store partially processed tomogram

  Alloc3D(image_size,
          &afTemp,
          &aaafTemp);

  RealNum A, B;        // let the user know what A B coefficients were used

  A = ApplyGauss3D(image_size,
                   aaafSource,
                   aaafDest,         // <-- save result here
                   aaafMask,
                   sigma_a,
                   filter_truncate_ratio,
                   filter_truncate_threshold,
                   true,
                   &cerr);

  B = ApplyGauss3D(image_size,
                   aaafSource,
                   aaafTemp,         // <-- save result here
                   aaafMask,
                   sigma_b,
                   filter_truncate_ratio,
                   filter_truncate_threshold,
                   true,
                   &cerr);

  // Subtract the second convolved signal from the first
  for (int iz = 0; iz < image_size[2]; iz++)
    for (int iy = 0; iy < image_size[1]; iy++)
      for (int ix = 0; ix < image_size[0]; ix++)
        aaafDest[iz][iy][ix] -= aaafTemp[iz][iy][ix];

  // Deallocate the temporary array
  Dealloc3D(image_size,
            &afTemp,
            &aaafTemp);

  // Report the A and B normalization coefficients to the caller?
  if (pA)
    *pA = A;
  if (pB)
    *pB = B;

} // ApplyDog3D(...,filter_truncate_ratio,filter_truncate_threshold,...)





template<class RealNum>
void
ApplyDogScaleFree3D(int const image_size[3], //source image size
                    RealNum ***aaafSource,   //source image
                    RealNum ***aaafDest,     //save results here
                    RealNum ***aaafMask,     //ignore voxels where mask==0
                    RealNum const sigma[3],  //Gaussian width in x,y,z drections
                    RealNum delta_sigma_over_sigma, //difference in Gauss widths
                    RealNum filter_truncate_ratio,
                    RealNum filter_truncate_threshold,
                    RealNum *pA = NULL,
                    RealNum *pB = NULL,
                    ostream *pReportProgress = NULL)
{

  if (filter_truncate_ratio < 0)
    //    filter_truncate_threshold = exp(-(1/2)*filter_truncate_ratio^2);
    //    -> filter_truncate_ratio^2 = -2*log(filter_truncate_threshold)
    filter_truncate_ratio = sqrt(-2*log(filter_truncate_threshold));

  ApplyDogScaleFree3D(image_size,
                      aaafSource,
                      aaafDest,
                      aaafMask,
                      sigma,
                      delta_sigma_over_sigma,
                      filter_truncate_ratio,
                      pA,
                      pB,
                      pReportProgress);

} // ApplyDogScaleFree3D(..,filter_truncate_ratio,filter_truncate_threshold,...)






void
HandleGGauss(Settings settings,
             MrcSimple &tomo_in,
             MrcSimple &tomo_out,
             MrcSimple &mask,
             float voxel_width[3])
{
  Filter3D<float, int>
    filter = GenFilterGenGauss3D(settings.width_a,
                                 settings.m_exp,
                                 settings.filter_truncate_ratio,
                                 settings.filter_truncate_threshold,
                                 static_cast<float*>(NULL),
                                 &cerr);

  filter.Apply(tomo_in.header.nvoxels,
               tomo_in.aaafI,
               tomo_out.aaafI,
               mask.aaafI,
               false,
               &cerr);

} //HandleGGauss()




void
HandleGauss(Settings settings,
            MrcSimple &tomo_in,
            MrcSimple &tomo_out,
            MrcSimple &mask,
            float voxel_width[3])
{
  cerr << "filter_type = Gaussian\n";

  float A; // let the user know what A coefficient was used

  A = ApplyGauss3D(tomo_in.header.nvoxels,
                   tomo_in.aaafI,
                   tomo_out.aaafI,   // <-- save result here
                   mask.aaafI,
                   settings.width_a,
                   settings.filter_truncate_ratio,
                   settings.filter_truncate_threshold,
                   true,
                   &cerr);

      
  cerr << " Filter Used:\n"
    " h(x,y,z)   = A*exp(-0.5*((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^))\n"
    //" h(x,y,z)   = A*exp(-((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))\n"
    " ... where  A = " << A << "\n" 
    "   (a_x, a_y, a_z) = "
       << "(" << settings.width_a[0]
       << " " << settings.width_a[1]
       << " " << settings.width_a[2] << ")\n";
  cerr << " You can plot a slice of this function\n"
       << "     in the X direction using:\n"
       << " draw_filter_1D.py -gauss "
       << A << " " << settings.width_a[0] << endl;
  cerr << " and in the Y direction using:\n"
       << " draw_filter_1D.py -gauss "
       << A << " " << settings.width_a[1] << endl;
  cerr << " and in the Z direction using:\n"
       << " draw_filter_1D.py -gauss "
       << A << " " << settings.width_a[2] << endl;

} //HandleGauss()




void
HandleDogg(Settings settings,
           MrcSimple &tomo_in,
           MrcSimple &tomo_out,
           MrcSimple &mask,
           float voxel_width[3])
{
  cerr << "filter_type = Difference-of-Generalized-Gaussians (DOGG)\n";

  Filter3D<float, int> filter;
  float A, B;       // let the user know what A B coefficients were used

  filter = GenFilterDogg3D(settings.width_a,//"a" parameter in formula
                           settings.width_b,//"b" parameter in formula
                           settings.m_exp,  //"m" parameter in formula
                           settings.n_exp,  //"n" parameter in formula
                           settings.filter_truncate_ratio,
                           settings.filter_truncate_threshold,
                           &A,
                           &B,
                           &cerr);

  filter.Apply(tomo_in.header.nvoxels,
               tomo_in.aaafI,
               tomo_out.aaafI,
               mask.aaafI,
               false,
               &cerr);

} // HandleDogg()



void
HandleDoggXY(Settings settings,
             MrcSimple &tomo_in,
             MrcSimple &tomo_out,
             MrcSimple &mask,
             float voxel_width[3])
{

  cerr << "filter_type = Difference-of-Generalized-Gaussians in the XY plane\n";
  // Generate a filter
  //
  //   h(x,y,z) = h_xy(x,y) * h_z(z)
  //
  //Take advantage of the fact that the filter we
  //are using (ie. the function we are convolving with the source image),
  //is the product of a function of X,Y, with a function of Z.
  //This makes it a "seperable" filter:  We can perform the filters in the Z
  //direction, followed by filtering the result in the X & Y directions.
  //This reduces the computation by a factor of O(filter.halfwidth[2]))
  //(A 1-D convolution followed by a 2-D convolution is much faster per 
  // voxel than a full 3-D convolution.)

  // First, generate the filter in the Z direction:

  int filter_truncate_halfwidthZ = -1;
  Filter1D<float, int> filterZ;
  if (settings.filter_truncate_ratio > 0.0)
    filter_truncate_halfwidthZ = floor(settings.width_a[2] *
                                       settings.filter_truncate_ratio);
  else
    filter_truncate_halfwidthZ = floor(settings.width_a[2] *
                                       sqrt(-2*log(settings.filter_truncate_threshold)));

  filterZ = GenFilterGauss1D(settings.width_a[2],
                             filter_truncate_halfwidthZ);

  float C;       // let the user know what C coefficient was used
  C = filterZ.afH[filter_truncate_halfwidthZ]; //(C=peak height located at the
                                              //   middle of the array)

  // Then generate the filter in the XY directions

  Filter2D<float, int> filterXY;
  float A, B;       // let the user know what A B coefficients were used

  filterXY = GenFilterDogg2D(settings.width_a,//"a" parameter in formula
                             settings.width_b,//"b" parameter in formula
                             settings.m_exp,  //"m" parameter in formula
                             settings.n_exp,  //"n" parameter in formula
                             settings.filter_truncate_ratio,
                             settings.filter_truncate_threshold,
                             &A,
                             &B);

  // Precompute the effect of the filter in the Z direction.

  // Create temporary 1-D arrays to perform the filter in the Z-direction:
  float *afIZorig = new float [tomo_in.header.nvoxels[2]];
  float *afIZnew  = new float [tomo_in.header.nvoxels[2]];
  float *afMask         = NULL;
  if (mask.aaafI)
    afMask       = new float [tomo_in.header.nvoxels[2]];

  // Then apply the filter in the Z direction
  // (and store the filtered 3-D image in the original tomogram array)
  for (int ix = 0; ix < tomo_in.header.nvoxels[0]; ix++) {
    for (int iy = 0; iy < tomo_in.header.nvoxels[1]; iy++) {
      for (int iz = 0; iz < tomo_in.header.nvoxels[2]; iz++) {
        afIZorig[iz] = tomo_in.aaafI[iz][iy][ix];
        if (afMask)
          afMask[iz] = mask.aaafI[iz][iy][ix];
      }
      filterZ.Apply(tomo_in.header.nvoxels[2],
                    afIZorig,
                    afIZnew,
                    afMask,
                    true);

      //It would be wasteful to allocate a temporary array to store this
      //Instead store the result of the convolution in the original array:
      for (int iz = 0; iz < tomo_in.header.nvoxels[2]; iz++)
        tomo_in.aaafI[iz][iy][ix] = afIZnew[iz];
    } //for (int iy = 0; iy < tomo_in.header.nvoxels[1]; iy++) {
  } // for (int ix = 0; ix < tomo_in.header.nvoxels[0]; ix++) {
  delete [] afIZorig;
  delete [] afIZnew;
  if (afMask)
    delete [] afMask;

  // Now apply the filter in the X and Y directions:
  cerr << "  progress: processing plane#" << endl;
  for (int iz = 0; iz < tomo_in.header.nvoxels[2]; iz++) {
    cerr << "  " << iz+1 << " / " << tomo_in.header.nvoxels[2] << "\n";
    float **aafMaskXY = NULL;
    if (mask.aaafI)
      aafMaskXY = mask.aaafI[iz];
    filterXY.Apply(tomo_in.header.nvoxels,
                   tomo_in.aaafI[iz],
                   tomo_out.aaafI[iz],
                   aafMaskXY,
                   false);
  }

  cerr << " Filter Used:\n"
    //" h(x,y,z) = (h_a(x,y) - h_b(x,y)) * C * exp(-0.5*(z/s)^2)\n"
    " h(x,y,z) = (h_a(x,y) - h_b(x,y)) * C * exp(-(z/s)^2)\n"
    " h_a(x,y) = A*exp(-((x/a_x)^2 + (y/a_y)^2)^(m/2))\n"
    " h_b(x,y) = B*exp(-((x/b_x)^2 + (y/b_y)^2)^(n/2))\n"
    "  ... where  A = " << A << "\n"
    "             B = " << B << "\n" 
    "             C = " << C << "\n" 
    "             m = " << settings.m_exp << "\n"
    "             n = " << settings.n_exp << "\n" 
    "   (a_x, a_y) = "
       << "(" << settings.width_a[0]
       << " " << settings.width_a[1] << ")\n"
    "   (b_x, b_y) = "
       << "(" << settings.width_b[0]
       << " " << settings.width_b[1] << ")\n"
    "            s = " << settings.width_a[2] << endl;
  cerr << " You can plot a slice of this function\n"
       << "     in the X direction using:\n"
    " draw_filter_1D.py -dogg " << A*C << " " << B*C
       << " " << settings.width_a[0] << " " << settings.width_b[0]
       << " " << settings.m_exp << " " << settings.n_exp << endl;
  cerr << " and in the Y direction using:\n"
    " draw_filter_1D.py -dogg " << A*C << " " << B*C
       << " " << settings.width_a[1] << " " << settings.width_b[1]
       << " " << settings.m_exp << " " << settings.n_exp << endl;
  cerr << " and in the Z direction using:\n"
    " draw_filter_1D.py -gauss " << C   // * (A-B)   <--COMMENTING OUT (A-B)
       << " " << settings.width_a[2] << endl;
} //HandleDoggXY()












void
HandleDog(Settings settings,
          MrcSimple &tomo_in,
          MrcSimple &tomo_out,
          MrcSimple &mask,
          float voxel_width[3])
{
  cerr << "filter_type = Difference of Gaussians (DOG)\n";

  float A, B;

  ApplyDog3D(tomo_in.header.nvoxels,
             tomo_in.aaafI,
             tomo_out.aaafI,
             mask.aaafI,
             settings.width_a,
             settings.width_b,
             settings.filter_truncate_ratio,
             settings.filter_truncate_threshold,
             &A,
             &B);

  cerr << " Filter Used:\n"
    " h(x,y,z)   = h_a(x,y,z) - h_b(x,y,z)\n"
    //" h_a(x,y,z) = A*exp(-0.5*((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))\n"
    " h_a(x,y,z) = A*exp(-((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))\n"
    //" h_b(x,y,z) = B*exp(-0.5*((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))\n"
    " h_b(x,y,z) = B*exp(-((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))\n"
    "  ... where      A = " << A << "\n"
    "                 B = " << B << "\n" 
    "   (a_x, a_y, a_z) = "
       << "(" << settings.width_a[0]
       << " " << settings.width_a[1]
       << " " << settings.width_a[2] << ")\n"
    "   (b_x, b_y, b_z) = "
       << "(" << settings.width_b[0]
       << " " << settings.width_b[1]
       << " " << settings.width_b[2] << ")\n";
  cerr << " You can plot a slice of this function\n"
       << "     in the X direction using:\n"
    " draw_filter_1D.py -dog " << A << " " << B
       << " " << settings.width_a[0] << " " << settings.width_b[0] << endl;
  cerr << " and in the Y direction using:\n"
    " draw_filter_1D.py -dog " << A << " " << B
       << " " << settings.width_a[1] << " " << settings.width_b[1] << endl;
  cerr << " and in the Z direction using:\n"
    " draw_filter_1D.py -dog " << A << " " << B
       << " " << settings.width_a[2] << " " << settings.width_b[2] << endl;

} //HandleDog()





void
HandleDogScaleFree(Settings settings,
                   MrcSimple &tomo_in,
                   MrcSimple &tomo_out,
                   MrcSimple &mask,
                   float voxel_width[3])
{
  cerr << "filter_type = Fast Laplacian of Gaussians (LOG)\n"
       << "  (This will be approximated as a Difference of Gaussians,\n"
       << "   as explained below.)\n";
  //cerr << "filter_type = Difference of Gaussians Scale Free (DOGSF)\n";

  float A, B;

  ApplyDogScaleFree3D(tomo_in.header.nvoxels,
                      tomo_in.aaafI,
                      tomo_out.aaafI,
                      mask.aaafI,
                      settings.dogsf_width,
                      settings.delta_sigma_over_sigma,
                      settings.filter_truncate_ratio,
                      settings.filter_truncate_threshold,
                      &A,
                      &B,
                      &cerr);


  // REMOVE THIS CRUFT EVENTUALLY:
  //for (int iz = 0; iz < tomo_in.header.nvoxels[2]; iz++)
  //  for (int iy = 0; iy < tomo_in.header.nvoxels[1]; iy++)
  //    for (int ix = 0; ix < tomo_in.header.nvoxels[0]; ix++)
  //      // Note: the DOG filter we used above is (in this case)
  //      // approximation of the LOG filter, (Laplacian-of-Gaussian),
  //      // which is the second derivative of the image intensity
  //      // after convolution with a Gaussian.
  //      // The ApplyDogScaleFree3D() function is not aware of the physical
  //      // size of each voxel, so it measures distance in voxels.
  //      // To report the 2nd derivative (with respect to position) in physical
  //      // units (instead of voxels) we need to divide by voxel_width^2
  //
  //     tomo_out.aaafI[iz][iy][ix] /= SQR(voxel_width[0]);
  //
  //      // Here we assume the voxel width is the same in x,y,z directions
  //      // and equals voxel_width[0].


  cerr << " Filter Used:\n"
    " h(x,y,z)   = h_a(x,y,z) - h_b(x,y,z)\n"
    //" h_a(x,y,z) = A*exp(-0.5*((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))\n"
    " h_a(x,y,z) = A*exp(-((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))\n"
    //" h_b(x,y,z) = B*exp(-0.5*((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))\n"
    " h_b(x,y,z) = B*exp(-((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))\n"
    "  ... where      A = " << A << "\n"
    "                 B = " << B << "\n" 
    "   (a_x, a_y, a_z) = "
       << "(" << settings.width_a[0]
       << " " << settings.width_a[1]
       << " " << settings.width_a[2] << ")\n"
    "   (b_x, b_y, b_z) = "
       << "(" << settings.width_b[0]
       << " " << settings.width_b[1]
       << " " << settings.width_b[2] << ")\n";
  cerr << " You can plot a slice of this function\n"
       << "     in the X direction using:\n"
    " draw_filter_1D.py -dog " << A << " " << B
       << " " << settings.width_a[0] << " " << settings.width_b[0] << endl;
  cerr << " and in the Y direction using:\n"
    " draw_filter_1D.py -dog " << A << " " << B
       << " " << settings.width_a[1] << " " << settings.width_b[1] << endl;
  cerr << " and in the Z direction using:\n"
    " draw_filter_1D.py -dog " << A << " " << B
       << " " << settings.width_a[2] << " " << settings.width_b[2] << endl;

} // HandleDogUseDelta()




// Find all scale-invariant blobs in the image:
template<class RealNum>
void
BlobDog3D(int const image_size[3], //source image size
          RealNum ***aaafSource,   //source image
          RealNum ***aaafMask,     //ignore voxels where mask==0
          const vector<RealNum>& blob_widths, // list of radii to try (ordered)
          vector<array<int,3> >& minima_crds, // store minima x,y,z coords here
          vector<array<int,3> >& maxima_crds, // store maxima x,y,z coords here
          vector<RealNum>& minima_sigma, // corresponding radius for that minima
          vector<RealNum>& maxima_sigma, // corresponding radius for that maxima
          vector<RealNum>& minima_scores, // what was the blobs score?
          vector<RealNum>& maxima_scores, // (score = intensity after filtering)
          RealNum delta_sigma_over_sigma, //difference in Gauss widths parameter
          RealNum filter_truncate_ratio,     //how many sigma before truncating?
          RealNum filter_truncate_threshold, //decay in filter before truncating
          RealNum minima_threshold,       // discard unremarkable minima
          RealNum maxima_threshold,       // discard unremarkable maxima
          bool use_threshold_ratios=true, // threshold=ratio*best_score ?
          RealNum min_separation = 0.0,
          ostream *pReportProgress = NULL,
          RealNum ****aaaafI = NULL, //preallocated memory for filtered images
          RealNum **aafI = NULL)     //preallocated memory for filtered images
                                     //(sometimes the caller wants to access)
{
  
  if (filter_truncate_ratio <= 0) {
    assert(filter_truncate_threshold > 0.0);
    //    filter_truncate_threshold = exp(-(1/2)*filter_truncate_ratio^2);
    //    -> filter_truncate_ratio^2 = -2*log(filter_truncate_threshold)
    filter_truncate_ratio = sqrt(-2*log(filter_truncate_threshold));
  }

  BlobDog3D(image_size,
            aaafSource,
            aaafMask,
            blob_widths,
            minima_crds,
            maxima_crds,
            minima_sigma,
            maxima_sigma,
            minima_scores,
            maxima_scores,
            delta_sigma_over_sigma,
            filter_truncate_ratio,
            minima_threshold,
            maxima_threshold,
            use_threshold_ratios,
            min_separation,
            pReportProgress,
            aaaafI,
            aafI);

} //BlobDog3D(...,filter_truncate_ratio,filter_truncate_threshold,...)




void
HandleMinDistance(Settings settings,
                  MrcSimple &tomo_in,
                  MrcSimple &tomo_out,
                  MrcSimple &mask,
                  float voxel_width[3])
{
  fstream coords_file;
  coords_file.open(settings.in_coords_file_name.c_str(), ios::in);
  if (! coords_file)
    throw InputErr("Error: unable to open \""+
                   settings.in_coords_file_name +"\" for reading.\n");
  vector<array<int,3> > crds;
  while (coords_file) {
    float x, y, z;
    coords_file >> x;
    coords_file >> y;
    coords_file >> z;
    int ix, iy, iz;
    ix = static_cast<int>(x / voxel_width[0]);
    iy = static_cast<int>(y / voxel_width[1]);
    iz = static_cast<int>(z / voxel_width[2]);
    array<int, 3> ixiyiz;
    ixiyiz[0] = ix;
    ixiyiz[1] = iy;
    ixiyiz[2] = iz;
    crds.push_back(ixiyiz);
  }

  cerr << " ------ calculating distance to points in "
       << settings.in_coords_file_name << " ------\n"
       << endl;

  // At some point I was trying to be as general as possible and allowed
  // for the possibility that voxels need not be cubes (same width x,y,z)
  // Now, I realize that allowing for this possibility would slow
  // down the next step considerably, so I just assume cube-shaped voxels:
  float voxel_width_ = voxel_width[0];
  assert((voxel_width[0] == voxel_width[1]) &&
         (voxel_width[1] == voxel_width[2]));

  for (int iz=0; iz<tomo_out.header.nvoxels[2]; iz++) {
    cerr << "processing Z-plane " << iz+1 << " / "
         << tomo_out.header.nvoxels[2] << "\n";
    for (int iy=0; iy<tomo_out.header.nvoxels[1]; iy++) {
      for (int ix=0; ix<tomo_out.header.nvoxels[0]; ix++) {
        int rminsq_int = SQR(tomo_out.header.nvoxels[0] +
                             tomo_out.header.nvoxels[1] +
                             tomo_out.header.nvoxels[2]);

        if (mask.aaafI && (mask.aaafI[iz][iy][ix] == 0.0))
          continue;

        // Loop over all of the points and find the nearest one.
        // Calculate the distance to that point, and store that in the tomogram.
        // (Note: This is a terribly inefficient way to do this.
        //        A paint-bucket-like-tool using expanding spheres would
        //        be much faster and visit each voxel only once on average.)
        for (vector<array<int,3> >::const_iterator pxyz=crds.begin();
             pxyz != crds.end();
             pxyz++) {
          int rsqi=SQR(ix-(*pxyz)[0])+SQR(iy-(*pxyz)[1])+SQR(iz-(*pxyz)[2]);
          if (rsqi < rminsq_int)
            rminsq_int = rsqi;
        }
        float rmin = sqrt(rminsq_int * SQR(voxel_width_));
        tomo_out.aaafI[iz][iy][ix] = rmin;
      }
    }
  }
} //HandleMinDistance()







void
ApplySphereDecals(MrcSimple &tomo_in,
                  MrcSimple &tomo_out,
                  MrcSimple &mask,
                  vector<array<int,3> > &crds,
                  vector<float> &radii,
                  vector<float> &shell_thicknesses,
                  vector<float> &voxel_intensities_foreground,
                  float voxel_intensity_background = 0.0,
                  float voxel_intensity_background_rescale = 0.4,
                  bool voxel_intensity_foreground_normalize = false)
{

  tomo_out = tomo_in; //copy the voxels from the original image to tomo_out
  float tomo_ave  =  AverageArr(tomo_in.header.nvoxels,
                                tomo_in.aaafI,
                                mask.aaafI);
  float tomo_stddev = StdDevArr(tomo_in.header.nvoxels,
                                tomo_in.aaafI,
                                mask.aaafI);

  double score_sum_sq = 0.0;
  for (int i = 0; i < voxel_intensities_foreground.size(); i++)
    score_sum_sq += SQR(voxel_intensities_foreground[i]);

  double score_rms = sqrt(score_sum_sq / voxel_intensities_foreground.size());

  // The spherical shells will have brightnesses chosen according to their
  // scores. Rescale the tomogram intensities so that they are approximately
  // in the range of scores.  This way we can see both at the same time
  // (when viewing the tomogram using IMOD, for example)
  for (int iz = 0; iz < tomo_out.header.nvoxels[2]; iz++) {
    for (int iy = 0; iy < tomo_out.header.nvoxels[1]; iy++) {
      for (int ix = 0; ix < tomo_out.header.nvoxels[0]; ix++) {
        tomo_out.aaafI[iz][iy][ix] =
          ((tomo_out.aaafI[iz][iy][ix] - tomo_ave) / tomo_stddev)
          *
          score_rms
          *
          voxel_intensity_background_rescale;
      }
    }
  }


  for (int iz=0; iz<tomo_out.header.nvoxels[2]; iz++)
    for (int iy=0; iy<tomo_out.header.nvoxels[1]; iy++)
      for (int ix=0; ix<tomo_out.header.nvoxels[0]; ix++)
        tomo_out.aaafI[iz][iy][ix] += voxel_intensity_background;


  for (int i = 0; i < crds.size(); i++) {
    cerr << "processing coordinates " << i+1 << " / "
         << crds.size() << " : "
         << crds[i][0] << " " << crds[i][1] << " " << crds[i][2] << "\n";

    if (mask.aaafI &&
        (mask.aaafI[crds[i][2]][crds[i][1]][crds[i][0]] == 0.0))
      continue;

    int Rs = ceil(radii[i]-0.5);
    if (Rs < 0) Rs = 0;
    float Rssqr_max = SQR(radii[i]);
    float Rssqr_min = 0.0;
    if ((shell_thicknesses[i] > 0.0) && (radii[i] - shell_thicknesses[i] > 0.0))
      Rssqr_min = SQR(radii[i] - shell_thicknesses[i]);
    // Normalize the brightness of each sphere?
    // (ie by dividing the intensity by the number of voxels in the sphere)
    float imultiplier = 1.0;
    long nvoxelspersphere = 1;
    if (voxel_intensity_foreground_normalize) {
      nvoxelspersphere = 0;
      for (int jz = -Rs; jz <= Rs; jz++) {
        for (int jy = -Rs; jy <= Rs; jy++) {
          for (int jx = -Rs; jx <= Rs; jx++) {
            float rsqr = jx*jx + jy*jy + jz*jz;
            if ((Rssqr_min <= rsqr) && (rsqr <= Rssqr_max))
              nvoxelspersphere++;
          }
        }
      }
    }
    imultiplier = 1.0 / nvoxelspersphere;
    for (int jz = -Rs; jz <= Rs; jz++) {
      for (int jy = -Rs; jy <= Rs; jy++) {
        for (int jx = -Rs; jx <= Rs; jx++) {
          int rsqr = jx*jx + jy*jy + jz*jz;
          if (! ((Rssqr_min <= rsqr) && (rsqr <= Rssqr_max)))
            continue;
          else if ((crds[i][0] + jx < 0) ||
                   (crds[i][0] + jx >= tomo_out.header.nvoxels[0]) ||
                   (crds[i][1] + jy < 0) ||
                   (crds[i][1] + jy >= tomo_out.header.nvoxels[1]) ||
                   (crds[i][2] + jz < 0) ||
                   (crds[i][2] + jz >= tomo_out.header.nvoxels[2]))
            continue;
          else if (mask.aaafI
                   &&
                   (mask.aaafI[crds[i][2]+jz]
                              [crds[i][1]+jy]
                              [crds[i][0]+jx]
                    == 0.0))
            continue;
          else
            tomo_out.aaafI[crds[i][2]+jz]
                          [crds[i][1]+jy]
                          [crds[i][0]+jx] =
                               voxel_intensities_foreground[i] * imultiplier;
        }
      }
    }
  } //for (int i = 0; i < crds.size(); i++) {

} //ApplySphereDecals()





void
HandleSphereDecals(Settings settings,
                   MrcSimple &tomo_in,
                   MrcSimple &tomo_out,
                   MrcSimple &mask,
                   float voxel_width[3])
{
  // At some point I was trying to be as general as possible and allowed
  // for the possibility that voxels need not be cubes (same width x,y,z)
  // Now, I realize that allowing for this possibility would slow
  // down the next step considerably, so I just assume cube-shaped voxels:
  float voxel_width_ = voxel_width[0];
  assert((voxel_width[0] == voxel_width[1]) &&
         (voxel_width[1] == voxel_width[2]));

  vector<array<int,3> > crds;
  vector<float> radii;
  vector<float> shell_thicknesses;
  vector<float> scores;

  cerr << " ------ calculating distance to points in "
       << settings.in_coords_file_name << " ------\n"
       << endl;

  fstream coords_file;
  coords_file.open(settings.in_coords_file_name.c_str(), ios::in);
  if (! coords_file)
    throw InputErr("Error: unable to open \""+
                   settings.in_coords_file_name +"\" for reading.\n");

  bool custom_radii = false;
  while (coords_file) {
    string strLine;
    getline(coords_file, strLine);
    if (strLine.size() == 0)
      continue;
    stringstream ssLine(strLine);
    float x, y, z;
    ssLine >> x;
    ssLine >> y;
    ssLine >> z;
    int ix, iy, iz;
    ix = static_cast<int>(x / voxel_width[0]);
    iy = static_cast<int>(y / voxel_width[1]);
    iz = static_cast<int>(z / voxel_width[2]);
    array<int, 3> ixiyiz;
    ixiyiz[0] = ix;
    ixiyiz[1] = iy;
    ixiyiz[2] = iz;

    float radius = -1.0;
    float score = settings.sphere_decals_foreground;
    if (ssLine) { // Does the file contain a 4th column? (the radius)
      float _radius;
      ssLine >> _radius;
      // convert from physical distance to # of voxels:
      _radius /= voxel_width[0];
      if (ssLine)
        radius = _radius;
      custom_radii = true;
    }

    if (radius < 0) { //If file does not contain a 4th column
      radius = settings.sphere_decals_radius;
      if (settings.sphere_decals_radius < 0)
        radius = 0.5;   // (sphere will be 1 voxel wide by default)
    }
    if (settings.sphere_decals_radius >= 0) //override the radius ?
      radius = settings.sphere_decals_radius;
    else
      radius *= settings.sphere_decals_scale;

    float shell_thickness = settings.sphere_decals_shell_thickness;
    if (settings.sphere_decals_shell_thickness_is_ratio) {
      shell_thickness *= radius;
      if (shell_thickness < settings.sphere_decals_shell_thickness_min)
        shell_thickness = settings.sphere_decals_shell_thickness_min;
    }
    shell_thickness = settings.sphere_decals_shell_thickness;
    if (settings.sphere_decals_shell_thickness_is_ratio) {
      shell_thickness *= radius;
      // The spherical shells superimposed on the tomogram should be
      // at least 1 voxel wide in order to be visible to the user
      if (shell_thickness < settings.sphere_decals_shell_thickness_min)
        shell_thickness = 1.0;
    }

    if (ssLine) {
      float _score;
      ssLine >> _score;
      if (ssLine)
        score = _score;
    }

    if (settings.score_lower_bound > settings.score_upper_bound) {
      if (! ((settings.score_lower_bound <= score) ||
             (score <= settings.score_upper_bound)))
        // If the score is not sufficiently high (or low), skip this blob
        continue; 
    }

    if (! settings.sphere_decals_foreground_use_score)
      score = settings.sphere_decals_foreground;

    crds.push_back(ixiyiz);
    radii.push_back(radius);
    shell_thicknesses.push_back(shell_thickness);
    scores.push_back(score);

  } //while (coords_file) {...




  reverse(crds.begin(), crds.end());
  reverse(radii.begin(), radii.end());
  reverse(scores.begin(), scores.end());

  ApplySphereDecals(tomo_in,
                    tomo_out,
                    mask,
                    crds,
                    radii,
                    shell_thicknesses,
                    scores,
                    settings.sphere_decals_background,
                    settings.sphere_decals_background_scale,
                    settings.sphere_decals_foreground_norm);

} //HandleSphereDecals()





void
HandleBlobDetector(Settings settings,
                   MrcSimple &tomo_in,
                   MrcSimple &tomo_out,
                   MrcSimple &mask,
                   float voxel_width[3])
{
  vector<array<int,3> > minima_crds_int;
  vector<array<int,3> > maxima_crds_int;

  vector<float> minima_sigma;
  vector<float> maxima_sigma;

  vector<float> minima_scores;
  vector<float> maxima_scores;


  // Optional: Preallocate space for BlobDog3D()
  float*** aaaafI[3];
  float* aafI[3];

  Alloc3D(tomo_in.header.nvoxels,
          &(aafI[0]),
          &(aaaafI[0]));

  if (tomo_out.aaafI) {
    // Optional: Instead of allocating aaaafI[1] and aafI[1], borrow the memory
    // you've already allocated to tomo_out.aaafI and afI. (Goal: Save memory.)
    aafI[1] = tomo_out.afI;
    aaaafI[1] = tomo_out.aaafI;
  } else {
    Alloc3D(tomo_in.header.nvoxels,
            &(aafI[1]),
            &(aaaafI[1]));
  }

  Alloc3D(tomo_in.header.nvoxels,
          &(aafI[2]),
          &(aaaafI[2]));
  // This way we can save memory and also save the 
  // filtered image to a file which we can view using IMOD.

  BlobDog3D(tomo_in.header.nvoxels,
            tomo_in.aaafI,
            mask.aaafI,
            settings.blob_widths,  // a list of different sigma to try (ordered)
            minima_crds_int,  // store minima x,y,z coords here
            maxima_crds_int,  // store maxima x,y,z coords here
            minima_sigma, // corresponding radius for that minima
            maxima_sigma, // corresponding radius for that maxima
            minima_scores, // what was the blob's score?
            maxima_scores, // ("score" = intensity after filtering)
            settings.delta_sigma_over_sigma, //difference in Gauss widths parameter
            settings.filter_truncate_ratio,
            settings.filter_truncate_threshold,
            settings.score_upper_bound,
            settings.score_lower_bound,
            settings.score_bounds_are_ratios,
            settings.blob_min_separation,
            &cerr,
            aaaafI,
            aafI);


  long n_minima = minima_crds_int.size();
  long n_maxima = maxima_crds_int.size();
  vector<array<float,3> > minima_crds(n_minima);
  vector<array<float,3> > maxima_crds(n_maxima);

  // The user expects results in units of physical distance, not voxels.
  // Compensate for that now:
  for (int i = 0; i < n_minima; i++) {
    minima_crds[i][0] = minima_crds_int[i][0] * voxel_width[0];
    minima_crds[i][1] = minima_crds_int[i][1] * voxel_width[1];
    minima_crds[i][2] = minima_crds_int[i][2] * voxel_width[2];
    minima_sigma[i] *= voxel_width[0];      //Gaussian width has units of length
    // REMOVE THIS CRUFT EVENTUALLY:
    //minima_scores[i] /= SQR(voxel_width[0]);//"score"~=Laplacian_of_Gaussian and
  }
  for (int i = 0; i < n_maxima; i++) {
    maxima_crds[i][0] = maxima_crds_int[i][0] * voxel_width[0];
    maxima_crds[i][1] = maxima_crds_int[i][1] * voxel_width[1];
    maxima_crds[i][2] = maxima_crds_int[i][2] * voxel_width[2];
    maxima_sigma[i] *= voxel_width[0];      //Gaussian width has units of length
    // REMOVE THIS CRUFT EVENTUALLY:
    //maxima_scores[i] /= SQR(voxel_width[0]);//"score"~=Laplacian_of_Gaussian and
  }

  //string out_file_name_base = settings.out_file_name;
  //if ((EndsWith(settings.out_file_name, ".rec")) ||
  //    (EndsWith(settings.out_file_name, ".mrc")))
  //  out_file_name_base =
  //    settings.out_file_name.substr(0,
  //                                  settings.out_file_name.length()-4);

  if ((minima_crds_int.size() > 0) && (settings.blob_minima_file_name != "")) {
    fstream minima_file;
    minima_file.open(settings.blob_minima_file_name.c_str(), ios::out);
    if (! minima_file)
      throw InputErr("Error: unable to open \""+ settings.blob_minima_file_name +"\" for reading.\n");
    for (int i=0; i < minima_crds_int.size(); i++) {
      minima_file << minima_crds[i][0] << " "
                  << minima_crds[i][1] << " "
                  << minima_crds[i][2] << " "
                  << minima_sigma[i] / settings.blob_width_multiplier << " "
        // See comment at the end of this function:
                  << minima_scores[i] << "\n";
        // Here we assume the voxel width is the same in x,y,z directions
        // and equals voxel_width[0].
    }
  }

  if ((maxima_crds_int.size() > 0) && (settings.blob_maxima_file_name != "")) {
    fstream maxima_file;
    maxima_file.open(settings.blob_maxima_file_name.c_str(), ios::out);
    if (! maxima_file)
      throw InputErr("Error: unable to open \""+ settings.blob_maxima_file_name +"\" for reading.\n");
    for (int i=0; i < maxima_crds_int.size(); i++) {
      maxima_file << maxima_crds[i][0] << " "
                  << maxima_crds[i][1] << " "
                  << maxima_crds[i][2] << " "
                  << maxima_sigma[i] / settings.blob_width_multiplier << " "
        // See comment at the end of this function:
                  << maxima_scores[i] << "\n";
        // Here we assume the voxel width is the same in x,y,z directions
        // and equals voxel_width[0].
    }
  }

  // ---------------------------------------------------------------------
  // OPTIONAL: Now render the final image with all of the spherical shells
  //           superimposed on top of the original image:
  //
  // Does the user want us to create a new image showing where the blobs are?
  if (tomo_out.aaafI) {

    vector<array<int,3> > display_crds_int(minima_crds_int);
    display_crds_int.insert(display_crds_int.end(),
                            maxima_crds_int.rbegin(),
                            maxima_crds_int.rend());

    //concatinate "maxima_sigma" with "minima_sigma"
    vector<float> display_sigma(minima_sigma);
    display_sigma.insert(display_sigma.end(),
                         maxima_sigma.rbegin(),
                         maxima_sigma.rend());

    //concatinate "maxima_scores" with "minima_scores"
    vector<float> display_scores(minima_scores);
    display_scores.insert(display_scores.end(),
                          maxima_scores.rbegin(),
                          maxima_scores.rend());

    vector<float> display_shell_thicknesses(display_crds_int.size());
    vector<float> display_radii(display_crds_int.size());
    for (int i = 0; i < display_crds_int.size(); i++) {
      // Choose the size of the hollow spheres around each object so they are
      // large enough that the original objects underneath are still visible.
      display_radii[i] = sqrt(3.0) * display_sigma[i] / voxel_width[0];
      display_shell_thicknesses[i] = settings.sphere_decals_shell_thickness;
      if (settings.sphere_decals_shell_thickness_is_ratio)
        display_shell_thicknesses[i] *= display_radii[i];
      display_radii[i] *= settings.sphere_decals_scale;
      // The spherical shells superimposed on the tomogram should be
      // at least 1 voxel wide in order to be visible to the user
      if (display_shell_thicknesses[i] < settings.sphere_decals_shell_thickness_min)
        display_shell_thicknesses[i] = 1.0;
    }


    ApplySphereDecals(tomo_in,
                      tomo_out,
                      mask,
                      display_crds_int,
                      display_radii,
                      display_shell_thicknesses,
                      display_scores,
                      settings.sphere_decals_background,
                      settings.sphere_decals_background_scale,
                      false);

  } //if (tomo_out.aaafI)



  Dealloc3D(tomo_in.header.nvoxels,
            &aafI[0],
            &aaaafI[0]);

  // Optional: Instead of allocating aaaafI[1] and aafI[1], borrow the memory
  // you've already allocated to tomo_out.aaafI and afI. (Goal: Save memory.)
  if (! tomo_out.aaafI) {  // <-- Was tomo_out.aaafI available?
    // However, if we didn't allocate it this memory, then don't delete it:
    Dealloc3D(tomo_in.header.nvoxels,
              &aafI[1],
              &aaaafI[1]);
  }

  Dealloc3D(tomo_in.header.nvoxels,
            &aafI[2],
            &aaaafI[2]);


  // REMOVE THIS CRUFT EVENTUALLY:
  //for (int iz = 0; iz < tomo_in.header.nvoxels[2]; iz++)
  //  for (int iy = 0; iy < tomo_in.header.nvoxels[1]; iy++)
  //    for (int ix = 0; ix < tomo_in.header.nvoxels[0]; ix++)
  //      // Here we assume the voxel width is the same in x,y,z directions
  //      // and equals voxel_width[0].
  //      tomo_out.aaafI[iz][iy][ix] /= SQR(voxel_width[0]);

} //HandleBlobDetector()






void
HandleThresholds(Settings settings,
                 MrcSimple &tomo_in,
                 MrcSimple &tomo_out,
                 MrcSimple &mask,
                 float voxel_width[3])
{
  cerr << "Applying thresholds" << endl;

  if (settings.out_thresh2_use_clipping_sigma) {
    float stddev_intensity = StdDevArr(tomo_out.header.nvoxels,
                                       tomo_out.aaafI,
                                       mask.aaafI);
    float ave_intensity = AverageArr(tomo_out.header.nvoxels,
                                     tomo_out.aaafI,
                                     mask.aaafI);

    // REMOVE THIS CRUFT LATER:
    //size_t *histY;
    //float *histX;
    //int nbins = -1;
    //float hist_bin_width = 1.0;
    //IntensityHistogramArr(&histX,
    //                      &histY,
    //                      nbins,
    //                      hist_bin_width,
    //                      tomo_out.header.nvoxels,
    //                      tomo_out.aaafI,
    //                      mask.aaafI);
    //for (int i=0; i < nbins; i++)
    //  cerr << histX[i] << " " << histY[i] << endl;
    //delete [] histX;
    //delete [] histY;

    settings.out_threshold_01_a = ave_intensity +
      settings.out_threshold_01_a * stddev_intensity;
    settings.out_threshold_01_b = ave_intensity +
      settings.out_threshold_01_b * stddev_intensity;
    cerr << "ave="<< ave_intensity <<", stddev="<<stddev_intensity << endl;
    cerr << "  Clipping intensities between ["
         << settings.out_threshold_01_a << ", "
         << settings.out_threshold_01_b << "]" << endl;
  }
  for (int iz=0; iz<tomo_out.header.nvoxels[2]; iz++) {
    for (int iy=0; iy<tomo_out.header.nvoxels[1]; iy++) {
      for (int ix=0; ix<tomo_out.header.nvoxels[0]; ix++) {
        if (! settings.use_dual_thresholds) {
          tomo_out.aaafI[iz][iy][ix] =
            Threshold2(tomo_out.aaafI[iz][iy][ix],
                       settings.out_threshold_01_a,
                       settings.out_threshold_01_b,
                       (settings.out_thresh2_use_clipping
                        ? settings.out_threshold_01_a
                        : settings.out_thresh_a_value),
                       (settings.out_thresh2_use_clipping
                        ? settings.out_threshold_01_b
                        : settings.out_thresh_b_value));
        }
        else
          tomo_out.aaafI[iz][iy][ix] =
            Threshold4(tomo_out.aaafI[iz][iy][ix],
                       settings.out_threshold_01_a,
                       settings.out_threshold_01_b,
                       settings.out_threshold_10_a,
                       settings.out_threshold_10_b,
                       settings.out_thresh_a_value,
                       settings.out_thresh_b_value);
      }
    }
  }
} //HandleThresholds()



void
HandleExtrema(Settings settings,
              MrcSimple &tomo_in,
              MrcSimple &tomo_out,
              MrcSimple &mask,
              float voxel_width[3])
{
  //By now the user has stored some filtered image in tomo_out.  If not, it
  //should have the same contents as tomo_in.(we initially set tomo_out=tomo_in)
  tomo_out.FindMinMaxMean();

  //default values disable the thresholds:
  float local_minima_threshold = settings.find_minima_threshold;
  float local_maxima_threshold = settings.find_maxima_threshold;
  //float local_minima_threshold =tomo_out.header.dmax;//keep all minima
  //float local_maxima_threshold =tomo_out.header.dmin;//keep all maxima

  vector<array<int, 3> > minima_crds_int;
  vector<array<int, 3> > maxima_crds_int;

  cerr << "---- searching for local minima & maxima ----\n";
  for (int iz=0; iz<tomo_out.header.nvoxels[2]; iz++) {
    cerr << "  " << iz+1 << " / " << tomo_out.header.nvoxels[2] << "\n";
    for (int iy=0; iy<tomo_out.header.nvoxels[1]; iy++) {
      for (int ix=0; ix<tomo_out.header.nvoxels[0]; ix++) {
        // Search the 26 surrounding voxels to see if this voxel is
        // either a minima or a maxima
        bool is_minima = true;
        bool is_maxima = true;
        for (int jz = -1; jz <= 1; jz++) {
          for (int jy = -1; jy <= 1; jy++) {
            for (int jx = -1; jx <= 1; jx++) {
              if ((ix+jx <0) || (ix+jx >= tomo_out.header.nvoxels[0]) ||
                  (iy+jy <0) || (iy+jy >= tomo_out.header.nvoxels[1]) ||
                  (iz+jz <0) || (iz+jz >= tomo_out.header.nvoxels[2]))
                continue;
              if (jx==0 && jy==0 && jz==0)
                continue;
              if ((! mask.aaafI) && (mask.aaafI[iz+jz][iy+jy][ix+jx] == 0)){
                is_minima = false;
                is_maxima = false;
                continue;
              }
              if (tomo_out.aaafI[iz+jz][iy+jy][ix+jx] <=
                  tomo_out.aaafI[iz][iy][ix])
                is_minima = false;
              if (tomo_out.aaafI[iz+jz][iy+jy][ix+jx] >=
                  tomo_out.aaafI[iz][iy][ix])
                is_maxima = false;
            }
          }
        }
        if (is_minima && settings.find_minima) {
          if (((! mask.aaafI) || (mask.aaafI[iz][iy][ix] != 0)) &&
              (tomo_out.aaafI[iz][iy][ix] < local_minima_threshold)) {
            array<int, 3> ixiyiz;
            ixiyiz[0] = ix;
            ixiyiz[1] = iy;
            ixiyiz[2] = iz;
            minima_crds_int.push_back(ixiyiz);
          }
        }
        if (is_maxima && settings.find_maxima) {
          if (((! mask.aaafI) || (mask.aaafI[iz][iy][ix] != 0)) &&
              (tomo_out.aaafI[iz][iy][ix] > local_maxima_threshold)) {
            array<int, 3> ixiyiz;
            ixiyiz[0] = ix;
            ixiyiz[1] = iy;
            ixiyiz[2] = iz;
            maxima_crds_int.push_back(ixiyiz);
          }
        }
        assert(! (is_minima && is_maxima));
      } //for (int ix=0; ix<tomo_out.header.nvoxels[0]; ix++) {
    } //for (int iy=0; iy<tomo_out.header.nvoxels[1]; iy++) {
  } //for (int iz=0; iz<tomo_out.header.nvoxels[2]; iz++) {

  if (settings.blob_min_separation > 0.0) {
    if (settings.find_minima_occlusion_radius > 0.0) {
      vector<float> minima_sigma(minima_crds_int.size(),
                                 settings.find_minima_occlusion_radius/sqrt(2));
      vector<float> minima_scores(minima_crds_int.size());
      for (size_t i = 0; i < minima_crds_int.size(); i++) {
        int ix = minima_crds_int[i][0];
        int iy = minima_crds_int[i][1];
        int iz = minima_crds_int[i][2];
        minima_scores[i] = tomo_out.aaafI[iz][iy][ix];
      }
      DiscardOverlappingBlobs(minima_crds_int,
                              minima_sigma, 
                              minima_scores,
                              false,
                              settings.blob_min_separation,
                              tomo_in.header.nvoxels,
                              &cerr);
    }
    if (settings.find_maxima_occlusion_radius > 0.0) {
      vector<float> maxima_sigma(maxima_crds_int.size(),
                                 settings.find_maxima_occlusion_radius/sqrt(2));
      vector<float> maxima_scores(maxima_crds_int.size());
      for (size_t i = 0; i < maxima_crds_int.size(); i++) {
        int ix = maxima_crds_int[i][0];
        int iy = maxima_crds_int[i][1];
        int iz = maxima_crds_int[i][2];
        maxima_scores[i] = tomo_out.aaafI[iz][iy][ix];
      }
      DiscardOverlappingBlobs(maxima_crds_int,
                              maxima_sigma, 
                              maxima_scores,
                              true,
                              settings.blob_min_separation,
                              tomo_in.header.nvoxels,
                              &cerr);
    }

  } //if (settings.blob_min_separation > 0)

  //string out_file_name_base = settings.out_file_name;
  //if ((EndsWith(settings.out_file_name, ".rec")) ||
  //    (EndsWith(settings.out_file_name, ".mrc")))
  //  out_file_name_base =
  //    settings.out_file_name.substr(0,
  //                                  settings.out_file_name.length()-4);
  if ((minima_crds_int.size()) > 0 && settings.find_minima) {
    fstream minima_file;
    minima_file.open(settings.find_minima_file_name.c_str(), ios::out);
    if (! minima_file)
      throw InputErr("Error: unable to open \""+ settings.find_minima_file_name +"\" for reading.\n");
    for (int i=0; i < minima_crds_int.size(); i++)
      minima_file << minima_crds_int[i][0] * voxel_width[0] << " "
                  << minima_crds_int[i][1] * voxel_width[1] << " "
                  << minima_crds_int[i][2] * voxel_width[2] << "\n";
  }
  if ((maxima_crds_int.size() > 0) && settings.find_maxima) {
    fstream coords_file;
    coords_file.open(settings.find_maxima_file_name.c_str(), ios::out);
    if (! coords_file)
      throw InputErr("Error: unable to open \""+ settings.find_maxima_file_name +"\" for reading.\n");
    for (int i=0; i < maxima_crds_int.size(); i++)
      coords_file << maxima_crds_int[i][0] * voxel_width[0] << " "
                  << maxima_crds_int[i][1] * voxel_width[1] << " "
                  << maxima_crds_int[i][2] * voxel_width[2] << "\n";
  }
} //HandleExtrema()





#ifndef DISABLE_TEMPLATE_MATCHING


      // First, lets characterize the problem we are trying to solve
      // graphically and mathematically:
      //
      // 1-dimensional example:
      //
      //  I_i
      //            original image voxel intensities
      //  /|\                 __   
      //   |                 |  |   __
      //   |            __   |  |  |  |
      //   |           |  |  |  |__|  |
      //   |           |  |__|        |   __
      //   |           |              |  |  |
      //   |           |              |__|  |
      //   |      __   |                    |
      //   |     |  |__|                    |__
      //   |   __|                             |____
      //   |_____________________________________________\  i (voxel position)
      //        .. i0-2  i0-1  i0  i0+1  i0+2 ...        /
      //
      //
      // Notation:
      //  i = index (this integer corresponds to a voxel from the template)
      //      Don't worry about the fact that this is a 3-dimensional image.
      //      Pretend that images are 1-dimensional, and i is an integer.
      // I_ = the collection of intensities of voxels from the original image
      //    = (I_1, I_2, I_3, ... )  in vector notation
      //
      //
      // p_ = intensities of voxels from the image after translation by -i0
      //    = (  p_1,      p_2,      p_3,    ... )
      //    = ( I_(1-i0), I_(2-i0), I_(3-i0) ... )
      //
      // This will shift the image so that voxels that appear nearby i0 are
      // now located nearby 0.  (I.E, it will center the image at position i0)
      //
      //
      //  [Alternate notation:  p_ =  T(-i0) I_   
      //  where "T(-i0)" is the Translational Shift Operator acting on vector I_
      //   T(-i0) I_  =  (I_(1-i0), I_(2-i0), I_(3-i0), ... ) ]
      //
      //
      // Now, we want to compare this with the shape of another function, q_i
      // which stores the shape of template we are looking for in our image.
      // 
      //  q_i
      //
      //  /|\             template intensities  
      //   |                 __,--------.__              
      //   |              __/              \__
      //   |             /                    \
      //   |           _/                      \_
      //        ___,--'                          '--.____
      //                               
      //  /_________________________________________________\ i (voxel position)
      //  \            -4 -3 -2 -1  0  1  2  3  4           /
      //
      // q_ = the collection of intensities of voxels from the template
      //      which we are comparing with the image.
      //    = (..., q_-2, q_-1, q_0, q_1, q_2, q_3, ... )  in vector notation
      //
      // PROBLEM:
      //
      // Unfortunatley, the original image will be noisy, and the overall
      // brightness and contrast of the object will vary significantly from
      // image to image (and sometimes at different locations in the same image)
      // For example, we cannot say what the density (brightness) of, say,
      // a ribosome in a cell should be or even that it should be,
      // say 10% more dense than its surroundings.
      //
      // To see if the template fits the shape of the image at all,
      //   (...at this location in the image)
      // we slide and stretch the template function in the vertical
      // direction until it agrees with the original image as much as possible.
      // 
      // Mathematically, our goal is to solve for the scalars b,c which minimize
      // the sum-of-squares-difference between the intensities of the image and
      // the template, after scaling and intensity offsets have been applied
      //
      // Notation:
      // We can scale the intensity of voxels in the template
      // by multiplying the vector q_ by a scalar (which we call "c")
      // Likewise can offset intensity of voxels by adding a vector b*J
      // where "b" is a scalar
      //   and "J_" is a vector filled with 1s,  J_ = (1, 1, 1, ...)
      //
      // In this notation, we see that the goal is to choose b, c to minimize:
      //
      //   variance  =  |   p_   -   c*(q_ + b*J_) |^2
      //                                             
      //                   /|\      /|\     /|\
      //                    |        |       |
      //                    |        |       |
      //               translated  scaling intensity
      //                  image            offset
      //               (centered 
      //                 at i0)
      //
      // The "variance" is the mean squared difference between
      //   the intensity of the voxels in our source image,  and ...
      //   the intensity of our voxels in our template after shifting & scaling
      //
      // After we have solved for b and c, we will consider the template
      // to be a good match (at this location) if:
      //    1) the variance is small
      //    2) c is large
      //
      // Equivalently, we define the scoring function
      //
      //    c / sqrt(variance)
      //
      //    (A good match <--> lambda is small)
      //
      // First, we have to find an easy way to find the optimal values of b & c:
      //
      // --- Solution Method: ---
      //
      // First shift both the image and the template up and down, so that
      // have the same average value (considering voxels nearby)
      // Introduce :
      //   P_ = T(-i0) p_ - <p_> J  = ( p_1-<p_>, p_2-<p_>, p_3-<p_>, ... )
      //   Q_ =        q_ - <q_> J  = ( q_1-<q_>, q_2-<q_>, q_3-<q_>, ... )
      //
      // Now the goal is to find c to minimize
      //
      //    | P_ - c * Q_ |^2
      //
      // The optimal solution is to choose c so that p'_ - c*q'_
      // lies along the dotted line below
      //
      //              _.
      //        P_    /|
      //             / :
      //            /  :
      //           /   : P_ - c*Q_
      //          /    : (optimal c)
      //         /     : 
      //        /      :
      //       /       :
      //      /        :
      //     /         :
      //    /__________________\ Q_
      //                       /
      //
      //    __________\ 
      //              / c*Q_
      //                         
      //         where |c*Q_|^2 = <P_, P_> - <P_,Q_>^2/|Q_|^2
      //
      // SOLUTION:
      //
      //        c = <P_, Q_> / |Q_|^2
      //
      // variance = | P_  -  c Q_ |^2
      //          = < (P_ - c Q_), (P_ - c Q_) >
      //          = <P_, P_> - 2*c*<P_, Q_> + c^2*<Q_, Q_>
      //          = <P_, P_> - <P_,Q_>^2/|Q_|^2
      //
      // ...where we introduce the weighted inner product ("dot product")
      // between vectors P_ and Q_:
      //
      // <P_, Q_>  =  sum_i( w_i * P_i * Q_i )
      //
      // The norm (vector length) is also weighted by w_i:
      // | Q |^2 = <Q_, Q_> = sum_i( w_i * Q_i * Q_i )
      //
      //
      //
      // Note: Why do we use weighted inner products (w_i)?
      //
      // We only want to match the template to a region in the image locally,
      // ie., in the viscinity of i0 (see above)
      // To insure this, we weight all the terms in all the sumations by
      //       w_i  = a function of i which is large when i is near 0,
      //              and zero elsewhere
      //       w_i is usually a Gaussian whose width is at least as wide 
      //           as template we are fitting:
      //       w_i = A * exp(-(i/s)^2)
      //
      // ... and we only perform the sum in regions where w_i is non-zero



void
HandleTemplateGGauss(Settings settings,
                     MrcSimple &tomo_in,
                     MrcSimple &tomo_out,
                     MrcSimple &mask,
                     float voxel_width[3])
{

  // Filter weights, w_i:
  Filter3D<float, int>
    w = GenFilterGenGauss3D(settings.template_background_radius,
                            settings.n_exp,
                            //template_profile.halfwidth);
                            settings.filter_truncate_ratio,
                            settings.filter_truncate_threshold,
                            static_cast<float*>(NULL),
                            static_cast<ostream*>(NULL));
  // GenFilterGenGauss3D creates normalized gaussians with integral 1.
  // That's not what we want for the weights, w_i:
  // The weights w_i should be 1 in the viscinity we care about, and 0 outside
  // So, we can undo this normalization by dividing all the w_i weights
  // by their maximum value max(w_i)    (at the central peak, at halfwidth).
  // This will mean the maximum w_i is 1, and the others decay to 0, as we want
  float wpeak = w.aaafH[w.halfwidth[2]][w.halfwidth[1]][w.halfwidth[0]];
  w.MultiplyScalar(1.0 / wpeak);

  // Template function, Q:
  Filter3D<float, int>
    Q = GenFilterGenGauss3D(settings.width_a, //"a" parameter in formula
                            settings.m_exp,   //"m" parameter in formula
                            w.halfwidth,
                            static_cast<float*>(NULL),
                            static_cast<ostream*>(NULL));

  // Q_ = q_ - <q_>
  float qave = Q.Average(w.aaafH);
  Q.AddScalar(-qave);

  // It's also convenient to scale the voxel intensities in Q_ so that the 
  // range of values, |Q_|^2, (=the standard deviation) is also unity.
  // This makes it meaningful to compare "c" with "variance" directly.
  // (If c/sqrt(variance), then it means that the fit is

  // Now calculate Q_dot_Q
  float Q_dot_Q = Q.SumSqr(w.aaafH);
  Q.MultiplyScalar(1.0/sqrt(Q_dot_Q));
  Q_dot_Q = Q.SumSqr(w.aaafH); // = |Q_|^2  (same everywhere)
  assert(abs(Q_dot_Q - 1.0) < 0.001);  // (should equal 1 but may vary due to roundoff error of 32bit floats)


  cerr << " ------ Calculating the average of nearby voxels: ------\n";

  float template_background_sigma[3];
  for (int d = 0; d < 3; d++)
    template_background_sigma[d] = settings.template_background_radius[d] / sqrt(2);

  // P = original image (after subtracting average nearby intensities):

  MrcSimple P = tomo_in; //this will take care of allocating the array

  // First, let's calculate the weighted average voxel intensity in the 
  // source image
  if (settings.n_exp == 2.0) {

    // then do it the fast way with seperable (ordinary) Gaussian filters
    ApplyGauss3D(tomo_in.header.nvoxels,
                 tomo_in.aaafI,
                 P.aaafI,    // <-- save result here
                 mask.aaafI,
                 template_background_sigma,
                 settings.filter_truncate_ratio,
                 settings.filter_truncate_threshold,
                 true);
  }
  else {
    w.Apply(tomo_in.header.nvoxels,
            tomo_in.aaafI,
            P.aaafI,    // <-- save result here
            mask.aaafI,
            true,
            &cerr);
  }

  // Subtract the average value from the image intensity, and store in P:
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
        P.aaafI[iz][iy][ix] = tomo_in.aaafI[iz][iy][ix]-P.aaafI[iz][iy][ix];

  cerr << "\n"
    " ------ Calculating <P_, Q_> / |Q_|^2 ... -------\n" << endl;

  // Now calculate <P_, Q_>
  // Recall <P_, Q_> = sum_i (w_i * Q_i * P_i)
  // precompute w_i * Q_i and call it "Q_times_w_i"
  //        <P_, Q_> = sum_i (Q_times_w_i * P_i)
  // Then we can reuse the code we normally use to calculate sum_i(A_i B_i)
  Filter3D<float, int> Q_times_w = w;
  for (int jz=-w.halfwidth[2]; jz<=w.halfwidth[2]; jz++)
    for (int jy=-w.halfwidth[1]; jy<=w.halfwidth[1]; jy++)
      for (int jx=-w.halfwidth[0]; jx<=w.halfwidth[0]; jx++)
        Q_times_w.aaafH[jz+w.halfwidth[2]]
                       [jy+w.halfwidth[1]]
                       [jx+w.halfwidth[0]]
          = 
          (w.aaafH[jz+w.halfwidth[2]]
                  [jy+w.halfwidth[1]]
                  [jx+w.halfwidth[0]]
           *
           Q.aaafH[jz+w.halfwidth[2]]
                  [jy+w.halfwidth[1]]
                  [jx+w.halfwidth[0]]
           );
               

  // Calculate P_dot_Q = <P_, Q_> = sum_i (Q_times_w_i * P_i)

  MrcSimple P_dot_Q = tomo_in;  //(note: this will allocate P_dot_Q.aaafI)

  Q_times_w.Apply(tomo_in.header.nvoxels,
                  P.aaafI,
                  P_dot_Q.aaafI, // <-- store result here
                  mask.aaafI,
                  false,
                  &cerr);


  //cerr << " ----- DEBUGGING: Calculating RMSE of fit directly -----\n"<< endl;

  //Q.SlowDebugScanTemplateError(P.header.nvoxels,
  //                             P.aaafI,
  //                             tomo_out.aaafI,
  //                             P_dot_Q.aaafI,
  //                             w.aaafH,
  //                             2.0,
  //                             mask.aaafI,
  //                             &cerr);
  //for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
  //  for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
  //    for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
  //      tomo_out.aaafI[iz][iy][ix] = 
  //        P_dot_Q.aaafI[iz][iy][ix] / tomo_out.aaafI[iz][iy][ix];




  // Calculate <P_, P_>

  // Because memory is so precious, we must reuse arrays whenever we can.
  // So we will now use "P" (which used to denote voxel intensity)
  // to store the square of intensity instead
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
        P.aaafI[iz][iy][ix] *= P.aaafI[iz][iy][ix];

  cerr << "\n"
    " ------ Calculating <P_, P_> ... ------\n" << endl;

  //   Original code:  Commented out because it requires too much memory:
  //MrcSimple P_dot_P = tomo_in; //(this takes care of allocating the array)
  //   Memory efficient, ugly code:
  // We're done with tomo_in.  Re-use this array to store P_dot_P temporarilly,
  // but give it a new name, to make the code slightly less confusing to read:
  float ***P_dot_P_aaafI = tomo_in.aaafI;

  w.Apply(tomo_in.header.nvoxels,
          P.aaafI,
          //P_dot_P.aaafI,
          P_dot_P_aaafI,  // <-- store result here
          mask.aaafI,
          false,
          &cerr);

  // Now calculate "c", and "rmse" (sqrt(variance))
  // We can calculate this from the P_dot_Q and P_dot_P arrays we have now.
  // Save the result in tomo_out, and save to a file.

  // First, write out RMSE, root-mean-squared-error  (sqrt(variance))
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
        //     variance = <P_, P_> - <P_,Q_>^2 / |Q_|^2
        //              = <P_, P_> - <P_,Q_>^2 / <Q_,Q_>

        //float variance = (P_dot_P.aaafI[iz][iy][ix]
        float variance = (P_dot_P_aaafI[iz][iy][ix]
                          -
                          SQR(P_dot_Q.aaafI[iz][iy][ix]) / Q_dot_Q);

        //Optional:
        //Compensate for dividing by w.aaafH[][][] by "wpeak" earlier.
        //This enables us to interpret variance as RMSE, (a.k.a. root-
        //mean-squared-error.  The formula above only calculates the
        //"mean" if w_i are normalized, which they were before we divided
        //them all by wpeak.  (Also: This insures that this fast method is
        //equivalent to the "SlowDebug" method commented out above.)
        //(Note: It is important for other reasons that w_i (w.aaafH[][][])
        //       were not normalized until this point.)

        variance *= wpeak;

        if (variance < 0.0)
          variance = 0.0;

        tomo_out.aaafI[iz][iy][ix] = sqrt(variance);
      } //for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
    } //for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
  } //for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)

  tomo_out.FindMinMaxMean();
  string rmse_file_name;
  if ((EndsWith(settings.out_file_name, ".rec")) ||
      (EndsWith(settings.out_file_name, ".mrc")))
    rmse_file_name =
      (settings.out_file_name.substr(0,
                                     settings.out_file_name.length()-4)
       + string("_rmse.mrc"));
  else 
    rmse_file_name = settings.out_file_name + string("_rmse.mrc");
  tomo_out.Write(rmse_file_name);


  // Then write out the "c" parameters to a file
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
        //            c = <P_, Q_> / |Q_|^2
        //              = <P_, Q_> / <Q_,Q_>
        float c = P_dot_Q.aaafI[iz][iy][ix] / Q_dot_Q;
        tomo_out.aaafI[iz][iy][ix] = c;
      } //for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
    } //for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
  } //for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
  tomo_out.FindMinMaxMean();
  //No need to write this file out now. It happens later. Commenting out:
  //tomo_out.Write(settings.out_file_name + string("_c.rec"));
  

  // Optional: Create a (3D) image file containing the filter we are using
  //           multiplied by the weight function.  Q_times_w is the function
  //           that we are actually convolving with our original image.
  //           It's sometimes useful to visualize it, so I save it as a file
  MrcSimple Q_times_w_mrc(Q_times_w.array_size,
                          Q_times_w.aaafH);
  Q_times_w_mrc.Write("Q_times_w.mrc");
  // Do the same for the filter (Q) as well.
  MrcSimple Q_mrc(Q.array_size,
                  Q.aaafH);
  Q_mrc.Write("Q.mrc");


  //{ // FOR DEBUGGING ONLY:
  //  int Q_rand_halfwidth[3];
  //  Q_rand_halfwidth[0] = floor(2 * w.halfwidth[0]);
  //  Q_rand_halfwidth[1] = floor(2 * w.halfwidth[1]);
  //  Q_rand_halfwidth[2] = floor(2 * w.halfwidth[2]);
  //  Filter3D<float, int>
  //    Q_rand = GenFilterGenGauss3D(settings.width_a,//"a" parameter in formula
  //                                 settings.m_exp,  //"m" parameter in formula
  //                                 Q_rand_halfwidth,
  //                                 static_cast<float*>(NULL),
  //                                 &cerr);
  //  RANDOM_INIT();
  //  MrcSimple Q_rand_mrc(Q_rand.array_size,
  //                       Q_rand.aaafH);
  //  Q_rand_mrc.Rescale01();
  //  for (int jz=-Q_rand.halfwidth[2]; jz<=Q_rand.halfwidth[2]; jz++) {
  //    for (int jy=-Q_rand.halfwidth[1]; jy<=Q_rand.halfwidth[1]; jy++) {
  //      for (int jx=-Q_rand.halfwidth[0]; jx<=Q_rand.halfwidth[0]; jx++) {
  //        Q_rand_mrc.aaafI[jz+Q_rand.halfwidth[2]]
  //                        [jy+Q_rand.halfwidth[1]]
  //                        [jx+Q_rand.halfwidth[0]] *= 0.5;
  //        Q_rand_mrc.aaafI[jz+Q_rand.halfwidth[2]]
  //                        [jy+Q_rand.halfwidth[1]]
  //                        [jx+Q_rand.halfwidth[0]] += 0.2*RANDOM_GAUSSIAN();
  //      }
  //    }
  //  }
  //  Q_rand_mrc.FindMinMaxMean();
  //  Q_rand_mrc.Write("Q_rand.mrc");
  //}

} //HandleTemplateGGauss()






void
HandleTemplateGauss(Settings settings,
                    MrcSimple &tomo_in,
                    MrcSimple &tomo_out,
                    MrcSimple &mask,
                    float voxel_width[3])
{
  // Filter weights, w_i:
  Filter3D<float, int>
    w = GenFilterGenGauss3D(settings.template_background_radius,
                            static_cast<float>(2.0), //settings.template_background_exponent,
                            //template_profile.halfwidth);
                            settings.filter_truncate_ratio,
                            settings.filter_truncate_threshold,
                            static_cast<float*>(NULL),
                            static_cast<ostream*>(NULL));
  // GenFilterGenGauss3D() creates normalized gaussians with integral 1.
  // That's not what we want for the weights, w_i:
  // The weights w_i should be 1 in the viscinity we care about, and 0 outside
  // So, we can undo this normalization by dividing all the w_i weights
  // by their maximum value max(w_i)    (at the central peak, at halfwidth).
  // This will mean the maximum w_i is 1, and the others decay to 0, as we want
  float wpeak = w.aaafH[w.halfwidth[2]][w.halfwidth[1]][w.halfwidth[0]];
  w.MultiplyScalar(1.0 / wpeak);

  // Template function, Q:
  Filter3D<float, int>
    q = GenFilterGenGauss3D(settings.width_a, //"a" parameter in formula
                            static_cast<float>(2.0), //settings.m_exp,
                            w.halfwidth,
                            static_cast<float*>(NULL),
                            static_cast<ostream*>(NULL));

  // Q_ = q_ - <q_>
  Filter3D<float, int> Q = q;
  float qave = q.Average(w.aaafH);
  Q.AddScalar(-qave);

  // It's also convenient to scale the voxel intensities in Q_ so that the 
  // range of values, |Q_|^2, (=the standard deviation) is also unity.
  // This makes it meaningful to compare "c" with "variance" directly.
  // (If c/sqrt(variance), then it means that the fit is

  // Now calculate Q_dot_Q
  float Q_dot_Q = Q.SumSqr(w.aaafH);
  Q.MultiplyScalar(1.0/sqrt(Q_dot_Q));

  // We might as well do the same thing for q and qave:
  q.MultiplyScalar(1.0/sqrt(Q_dot_Q));
  qave = q.Average(w.aaafH);
  Q_dot_Q = Q.SumSqr(w.aaafH); // = |Q_|^2  (same everywhere)
  assert(abs(Q_dot_Q - 1.0) < 0.001);  // (should equal 1 but may vary due to roundoff error of 32bit floats)



  cerr << "\n"
    " ------ Calculating the average of nearby voxels: ------\n";
  // P = original image (after subtracting average nearby intensities):

  MrcSimple P = tomo_in; //this will take care of allocating the array

  // First, let's calculate the weighted average voxel intensity in the 
  // source image
  //
  // THIS WORKS, BUT THERE IS A FASTER WAY.  COMMENTING OUT:
  //
  //w.Apply(tomo_in.header.nvoxels,
  //        tomo_in.aaafI,
  //        P.aaafI,   // <-- save result here
  //        mask.aaafI,
  //        true,
  //        &cerr);
  //
  // TRY THIS INSTEAD (works only when w is an ordinary Gaussian):

  float template_background_sigma[3];
  for (int d = 0; d < 3; d++)
    template_background_sigma[d] = settings.template_background_radius[d] / sqrt(2);

  ApplyGauss3D(tomo_in.header.nvoxels,
               tomo_in.aaafI,
               P.aaafI,   // <-- save result here
               mask.aaafI,
               template_background_sigma, //width of Gaussian
               settings.filter_truncate_ratio,
               settings.filter_truncate_threshold,
               true,
               &cerr);

  // Subtract the average value from the image intensity, and store in P:
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
        P.aaafI[iz][iy][ix] = tomo_in.aaafI[iz][iy][ix]-P.aaafI[iz][iy][ix];



  cerr << "\n"
    " ------ Calculating <P_, Q_> / |Q_|^2 ... -------\n" << endl;




  //  REMOVE THIS CRUFT:
  //
  //// Now calculate <P_, Q_>
  //// Recall <P_, Q_> = sum_i (w_i * Q_i * P_i)
  //// precompute w_i * Q_i and call it "Q_times_w_i"
  ////        <P_, Q_> = sum_i (Q_times_w_i * P_i)
  //// Then we can reuse the code we normally use to calculate sum_i(A_i B_i)
  //Filter3D<float, int> Q_times_w = w;
  //for (int jz=-w.halfwidth[2]; jz<=w.halfwidth[2]; jz++)
  //  for (int jy=-w.halfwidth[1]; jy<=w.halfwidth[1]; jy++)
  //    for (int jx=-w.halfwidth[0]; jx<=w.halfwidth[0]; jx++)
  //      Q_times_w.aaafH[jz+w.halfwidth[2]]
  //                     [jy+w.halfwidth[1]]
  //                     [jx+w.halfwidth[0]]
  //        = 
  //        (w.aaafH[jz+w.halfwidth[2]]
  //                [jy+w.halfwidth[1]]
  //                [jx+w.halfwidth[0]]
  //         *
  //         Q.aaafH[jz+w.halfwidth[2]]
  //                [jy+w.halfwidth[1]]
  //                [jx+w.halfwidth[0]]
  //         );
               



  // Calculate P_dot_Q = <P_, Q_> = sum_i (Q_times_w_i * P_i)
  MrcSimple P_dot_Q = tomo_in; //(this takes care of allocating the array)



  // The original code works, but there's a faster way.  Commenting out:
  //Q_times_w.Apply(tomo_in.header.nvoxels,
  //                P.aaafI,
  //                P_dot_Q.aaafI,
  //                mask.aaafI,
  //                false,
  //                &cerr);
  // Instead, use "ApplyGauss3D()", a fast version which only works for
  // Gaussians.  First we have to calculate the width of that Gaussian...
  float radius_Q_times_w[3];
  float sigma_Q_times_w[3];
  for (int d=0; d < 3; d++) {
    // The product of two Gaussians (widths a, b) is a Gaussian
    // whose width is narrower than the width of either Gaussian
    // e^(-(r/a)^2) * e^(-(r/b)^2) = e^(-(r/R)^2),   R=sqrt(1/(1/a^2+1/b^2))
    radius_Q_times_w[d] =
      sqrt(1.0 / (1.0 / SQR(settings.width_a[d])
                  +
                  1.0 / SQR(settings.template_background_radius[d])));
    sigma_Q_times_w[d] = radius_Q_times_w[d] / sqrt(2);
  }
  // We can use the faster "ApplyGauss3D()" to perform the convolution
  ApplyGauss3D(P.header.nvoxels, // = tomo_in.header.nvoxels,
               P.aaafI,
               P_dot_Q.aaafI,        // <-- save result here
               mask.aaafI,
               sigma_Q_times_w,
               settings.filter_truncate_ratio,
               settings.filter_truncate_threshold,
               false, // Don't normalize. Gaussian will have peak height=1
               &cerr);
  // The function above convolves the image, P, with a Gaussian whose
  // central peak has a height of 1.  What we want to do instead is convolve
  // it with a Gaussian whose central peak has a height equal to the height
  // of the "Q" gaussian.  This number is in the middle of the Q.aaafH array
  float qpeak = q.aaafH[Q.halfwidth[2]]
                       [Q.halfwidth[1]]
                       [Q.halfwidth[0]];


  //  REMOVE THIS CRUFT:
  //float q_times_w_peak = Q_times_w.aaafH[Q_times_w.halfwidth[2]]
  //                                      [Q_times_w.halfwidth[1]]
  //                                      [Q_times_w.halfwidth[0]];
  ////should be the same, since w peaked at 1
  //assert((abs(qpeak-qave)-abs(q_times_w_peak)) < 0.001*abs(q_times_w_peak));


  // Now multiply all of the voxels in the colvoved image (P_dot_Q)
  // by "qpeak" to compensate for using the wrong normalization earlier:
  for(int iz=0; iz<P_dot_Q.header.nvoxels[2]; iz++)
    for(int iy=0; iy<P_dot_Q.header.nvoxels[1]; iy++)
      for(int ix=0; ix<P_dot_Q.header.nvoxels[0]; ix++)
        P_dot_Q.aaafI[iz][iy][ix] *= qpeak;
  // (The speed gains are worth the extra code complexity)

  // We forgot to take into account the fact that Q_ = q_ - <q_>
  // and Q_times_w is the product of Q with another Gaussian.
  // The result is:
  // Q_times_w is not a Gaussian, but the difference of two Gaussians,
  //    one with width "radius_Q_times_w" (which we computed already), and
  //    one with width "settings.template_background_radius"
  // We have convolved P with q_times_w (before subtracting <q_> = "qave")
  // and stored it in P_dot_Q.  Now convolve P with qave*w and subtract
  // it from the result we obtained for P_dot_Q.
  // Again, we can use the faster "ApplyGauss3D()" to perform the convolution:
  ApplyGauss3D(P.header.nvoxels, // = tomo_in.header.nvoxels,
               P.aaafI,
               tomo_in.aaafI,        // <-- save result here
               mask.aaafI,
               template_background_sigma,
               settings.filter_truncate_ratio,
               settings.filter_truncate_threshold,
               false,  // Don't normalize. Gaussian will have peak height=1
               &cerr);

  for(int iz=0; iz<P_dot_Q.header.nvoxels[2]; iz++)
    for(int iy=0; iy<P_dot_Q.header.nvoxels[1]; iy++)
      for(int ix=0; ix<P_dot_Q.header.nvoxels[0]; ix++)
        P_dot_Q.aaafI[iz][iy][ix] -= qave * tomo_in.aaafI[iz][iy][ix];



  // Calculate <P_, P_>
  // (This should be easier than calculating <P_, Q_>)

  // Because memory is so precious, we must reuse arrays whenever we can.
  // So we will now use "P" (which used to denote voxel intensity)
  // to store the square of intensity instead
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
        P.aaafI[iz][iy][ix] *= P.aaafI[iz][iy][ix];




  cerr << "\n"
    " ------ Calculating <P_, P_> ... ------\n" << endl;


  //   Original code:  Commented out because it requires too much memory:
  //MrcSimple P_dot_P = tomo_in; //(this takes care of allocating the array)
  //   Memory efficient, ugly code:
  // We're done with tomo_in.  Re-use this array to store P_dot_P temporarilly,
  // but give it a new name, to make the code slightly less confusing to read:
  float ***P_dot_P_aaafI = tomo_in.aaafI;

  // The code below works too, but there's a faster way.  COMMENTING OUT:
  //w.Apply(tomo_in.header.nvoxels,
  //        P.aaafI,
  //        //P_dot_P.aaafI,
  //        P_dot_P_aaafI,
  //        mask.aaafI,
  //        false,
  //        &cerr);
  // Instead, we can use "ApplyGauss3D()" to perform the convolution quickly

  ApplyGauss3D(P.header.nvoxels, // = tomo_in.header.nvoxels,
               P.aaafI,
               //P_dot_P.aaafI,
               P_dot_P_aaafI,        // <-- save result here
               mask.aaafI,
               template_background_sigma,  //width of "w" Gaussian
               settings.filter_truncate_ratio,
               settings.filter_truncate_threshold,
               false,  // don't normalize. The "w" weights have a peak of 1
               &cerr);

  // Now calculate "c", and "rmse" (sqrt(variance))
  // We can calculate this from the P_dot_Q and P_dot_P arrays we have now.
  // Save the result in tomo_out, and save to a file.

  // First, write out RMSE, root-mean-squared-error  (sqrt(variance))
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
        //     variance = <P_, P_> - <P_,Q_>^2 / |Q_|^2
        //              = <P_, P_> - <P_,Q_>^2 / <Q_,Q_>

        //float variance = (P_dot_P.aaafI[iz][iy][ix]
        float variance = (P_dot_P_aaafI[iz][iy][ix]
                          -
                          SQR(P_dot_Q.aaafI[iz][iy][ix]) / Q_dot_Q);

        //Optional:
        //Compensate for dividing by w.aaafH[][][] by "wpeak" earlier.
        //This enables us to interpret variance as RMSE, (a.k.a. root-
        //mean-squared-error.  The formula above only calculates the
        //"mean" if w_i are normalized, which they were before we divided
        //them all by wpeak.
        //(Note: It is important for other reasons that w_i (w.aaafH[][][])
        //       were not normalized until this point.)

        variance *= wpeak;

        if (variance < 0.0)
          variance = 0.0;

        tomo_out.aaafI[iz][iy][ix] = sqrt(variance);
      } //for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
    } //for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
  } //for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
  tomo_out.FindMinMaxMean();
  tomo_out.FindMinMaxMean();
  string rmse_file_name;
  if ((EndsWith(settings.out_file_name, ".rec")) ||
      (EndsWith(settings.out_file_name, ".mrc")))
    rmse_file_name =
      (settings.out_file_name.substr(0,
                                     settings.out_file_name.length()-4)
       + string("_rmse.mrc"));
  else 
    rmse_file_name = settings.out_file_name + string("_rmse.mrc");
  tomo_out.Write(rmse_file_name);

  // Then write out the "c" parameters to a file
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
        //            c = <P_, Q_> / |Q_|^2
        //              = <P_, Q_> / <Q_,Q_>
        float c = P_dot_Q.aaafI[iz][iy][ix] / Q_dot_Q;
        tomo_out.aaafI[iz][iy][ix] = c;
      } //for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
    } //for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
  } //for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
  tomo_out.FindMinMaxMean();
  //No need to write this file out now. It happens later. Commenting out:
  //tomo_out.Write(settings.out_file_name + string("_c.rec"));



  // Optional: Create a (3D) image file containing the filter we are using
  //           multiplied by the weight function.  Q_times_w is the function
  //           that we are actually convolving with our original image.
  //           It's sometimes useful to visualize it, so I save it as a file
  //MrcSimple Q_times_w_mrc(Q_times_w.array_size,
  //                        Q_times_w.aaafH);
  //Q_times_w_mrc.Write("Q_times_w.mrc");
  // Do the same for the filter (Q) as well.
  MrcSimple Q_mrc(Q.array_size,
                      Q.aaafH);
  Q_mrc.Write("Q.mrc");

} //HandleTemplateGauss()




void
HandleLocalFluctuations(Settings settings,
                        MrcSimple &tomo_in,
                        MrcSimple &tomo_out,
                        MrcSimple &mask,
                        float voxel_width[3])
// Calculate the fluctuations of nearby voxel intensities
{

  // Filter weights, w_i:
  Filter3D<float, int>
    w = GenFilterGenGauss3D(settings.template_background_radius,
                            settings.template_background_exponent,
                            //template_profile.halfwidth);
                            settings.filter_truncate_ratio,
                            settings.filter_truncate_threshold,
                            static_cast<float*>(NULL),
                            static_cast<ostream*>(NULL));
  // GenFilterGenGauss3D() creates normalized gaussians with integral 1.
  // That's not what we want for the weights, w_i:
  // The weights w_i should be 1 in the viscinity we care about, and 0 outside
  // So, we can undo this normalization by dividing all the w_i weights
  // by their maximum value max(w_i)    (at the central peak, at halfwidth).
  // This will mean the maximum w_i is 1, and the others decay to 0, as we want
  float wpeak = w.aaafH[w.halfwidth[2]][w.halfwidth[1]][w.halfwidth[0]];
  w.MultiplyScalar(1.0 / wpeak);


  cerr << " ------ Calculating the average of nearby voxels: ------\n";
  // P = original image (after subtracting average nearby intensities):

  MrcSimple P = tomo_in; //this will take care of allocating the array

  float template_background_sigma[3];
  for (int d = 0; d < 3; d++)
    template_background_sigma[d] = settings.template_background_radius[d] / sqrt(2);

  // First, let's calculate the weighted average voxel intensity in the 
  // source image
  if (settings.template_background_exponent == 2.0) {

    // then do it the fast way with seperable (ordinary) Gaussian filters
    ApplyGauss3D(tomo_in.header.nvoxels,
                 tomo_in.aaafI,
                 P.aaafI,    // <-- save result here
                 mask.aaafI,
                 template_background_sigma,
                 settings.filter_truncate_ratio,
                 settings.filter_truncate_threshold,
                 true,
                 &cerr);
  }
  else {
    w.Apply(tomo_in.header.nvoxels,
            tomo_in.aaafI,
            P.aaafI,    // <-- save result here
            mask.aaafI,
            true,
            &cerr);
  }

  // Subtract the average value from the image intensity, and store in P:
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
        P.aaafI[iz][iy][ix] = tomo_in.aaafI[iz][iy][ix]-P.aaafI[iz][iy][ix];

  // Calculate <P_, P_>

  // Because memory is so precious, we must reuse arrays whenever we can.
  // So we will now use "P" (which used to denote voxel intensity)
  // to store the square of intensity instead
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
        P.aaafI[iz][iy][ix] *= P.aaafI[iz][iy][ix];

  cerr << "\n"
    " ------ Calculating fluctuations around that average ------\n" << endl;

  //   Original code:  Commented out because it requires too much memory:
  //MrcSimple P_dot_P = tomo_in; //(this takes care of allocating the array)
  //   Memory efficient, ugly code:
  // We're done with tomo_in.  Re-use this array to store P_dot_P temporarilly,
  // but give it a new name, to make the code slightly less confusing to read:
  float ***P_dot_P_aaafI = tomo_in.aaafI;

  if (settings.template_background_exponent == 2.0) {
    // then do it the fast way with seperable (ordinary) Gaussian filters
    ApplyGauss3D(tomo_in.header.nvoxels,
                 P.aaafI,
                 //P_dot_P.aaafI,
                 P_dot_P_aaafI,    // <-- save result here
                 mask.aaafI,
                 template_background_sigma,
                 settings.filter_truncate_ratio,
                 settings.filter_truncate_threshold,
                 true,
                 &cerr);
  }
  else {
    w.Apply(tomo_in.header.nvoxels,
            P.aaafI,
            //P_dot_P.aaafI,
            P_dot_P_aaafI,  // <-- store result here
            mask.aaafI,
            false,
            &cerr);
  }

  // Now calculate "rms" (sqrt(variance))
  // Save the result in tomo_out

  // Rrite out RMS variance from the average of nearby voxels:
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
        //     variance = <P_, P_>

        //float variance = (P_dot_P.aaafI[iz][iy][ix]
        float variance = P_dot_P_aaafI[iz][iy][ix];

        //Optional:
        //Compensate for dividing by w.aaafH[][][] by "wpeak" earlier.
        //This enables us to interpret variance as RMSE, (a.k.a. root-
        //mean-squared-error.  The formula above only calculates the
        //"mean" if w_i are normalized, which they were before we divided
        //them all by wpeak.  (Also: This insures that this fast method is
        //equivalent to the "SlowDebug" method commented out above.)
        //(Note: It is important for other reasons that w_i (w.aaafH[][][])
        //       were not normalized until this point.)

        variance *= wpeak;

        if (variance < 0.0)
          variance = 0.0;

        tomo_out.aaafI[iz][iy][ix] = sqrt(variance);
      } //for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
    } //for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
  } //for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)

  tomo_out.FindMinMaxMean();
  
} //HandleLocalFluctuations()



#endif //#ifndef DISABLE_TEMPLATE_MATCHING







#ifndef DISABLE_BOOTSTRAPPING

template<class RealNum>
void
ScrambleImage3D(int image_size[3],
                RealNum const ***aaafSource,
                  //(ellipsoidal) scramble radius in x,y,z directions:
                int const scramble_radius[3], 
                RealNum ***aaafDest)
// Scramble the contents of an image by replacing the current voxel with a
// randmly chosen voxel nearby
// (which lies within an ellipsoid specified by "scramble_radius")
{
  for (int iz=0; iz < image_size[2]; iz++) {
    for (int iy=0; iy < image_size[1]; iy++) {
      for (int ix=0; ix < image_size[0]; ix++) {
        RealNum r = 2.0;
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





void
HandleBootstrapDogg(Settings settings,
                     MrcSimple &tomo_in,
                     MrcSimple &tomo_out,
                     MrcSimple &mask);
{
  cerr << "filter_type = Difference-of-Generalized-Gaussians with BOOTSTRAPPING\n"
       << "              (BOOTSTRAP_DOGG)\n";

  Filter3D<float, int> filter;
  float A, B;       // let the user know what A B coefficients were used

  filter = GenFilterDogg3D(settings.width_a,//"a" parameter in formula
                           settings.width_b,//"b" parameter in formula
                           settings.m_exp,  //"m" parameter in formula
                           settings.n_exp,  //"n" parameter in formula
                           settings.filter_truncate_ratio,
                           settings.filter_truncate_threshold,
                           &A,
                           &B,
                           &cerr);

  if (settings.bs_ntests > 0) {

    // IN OLDER VERSIONS OF THE CODE, I WOULD ONLY CONSIDER VOXELS WHOSE
    // FILTERED VALUE EXCEEDED SOME USER-SPECIFIED THRESHOLD.
    // NOW, WE USE THE "MASK" INSTEAD.  IN OTHER WORDS, THE USER IS
    // RESPONSIBLE FOR RUNNING A FILTER ON THE INPUT IMAGE (OR USING SOME
    // OTHER METHOD) TO SELECT OUT THE VOXELS THEY WANT US TO CONSIDER.
    // SAVE THESE VOXELS IN A DIFFERENT FILE, AND LOAD IT USING "-mask"
    //
    // (If the don't, this procedure will take a really, really long time.)
    //
    if (mask.aaafI == NULL) {
      cerr << "WARNING: THIS FILTER IS VERY SLOW UNLESS YOU USE THE \"-mask\" ARGUMENT\n"
           << "         TO SPECIFY THE VOXELS YOU WANT TO CONSIDER.\n"
           << "         (BY DEFAULT, THIS PROGRAM USES ALL THE VOXELS).\n";
      for (int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
        for (int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
          for (int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
            vvGoodVoxels.push_back(ixiyiz);
            vCountFP.push_back(0);
          }
        }
      }
    }
    else { 
      for (int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
        for (int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
          for (int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
            if (mask.aaafI[ix][iy][iz] != 0.0) {
              vvGoodVoxels.push_back(ixiyiz);
              vCountFP.push_back(0);
            }
          }
        }
      }
    }
          
    //Initialize the random number generator we will need later:
    RANDOM_INIT(settings.bs_random_seed);
        
    for (int i_sim = 0; i_sim < settings.bs_ntests; i_sim++) {
      cerr <<
        "\n"
        " ---------------------------------------------------\n"
        "   False positives bootstrap simulation# " << isim+1
           << " / " << settings.bs_ntests << "\n"
           << endl;

      ScrambleImage(tomo_in.header.nvoxels,
                    tomo_in.aaafI,
                    settings.bs_scramble_radius,
                    aaafScrambled);

      for (i_vox = 0; i_vox < vvGoodVoxels.size(); i_vox++) {
        int ix = vGoodVoxels[i_vox][0];
        int iy = vGoodVoxels[i_vox][1];
        int iz = vGoodVoxels[i_vox][2];
        float g =
          filter.ApplyToVoxel(ix, iy, iz,
                              tomo_in.header.nvoxels,
                              //tomo_in.aaafI,
                              aaafScrambled,
                              mask.aaafI,
                              false);

        //if (g * settings.bs_threshold_sign >= //filter applied to scramble
        //    tomo_out.aaafI[iz][iy][ix]) //>= applied to orig image
        //  vCountFP[i_vox] += 1;

        if (g * settings.bs_threshold_sign >=
            settings.bs_threshold)
          vCountFP[i_vox] += 1;

      }
    } // for (int i_sim = 0; i_sim < settings.bs_ntests; i_sim++)

    // Now copy the measured probabilities into the output image:
    for (i_vox = 0; i_vox < vvGoodVoxels.size(); i_vox++) {
      int ix = vGoodVoxels[i_vox][0];
      int iy = vGoodVoxels[i_vox][1];
      int iz = vGoodVoxels[i_vox][2];
      // The tomogram should store the 1 - (probability of a false positive)
      // at this location.  In the code above, this probability
      // has been calculated by counting the number of times that a
      // (locally) scrambled image generated a signal as strong as the
      // result of applying the same filter to the unscrambled image there.
      tomo_out.aaafI[iz][iy][ix] =
        (1.0
         -
         static_cast<float>(vCountFP[i_vox]) / settings.fp_nsims);
    }

    ////Dealloc3D(tomo_in.header.nvoxels,
    ////        &aiCount,
    ////        &aaaiCount);

  } //if (settings.bs_ntests > 0)

} //HandleBootstrapDogg()
#endif //#ifndef DISABLE_BOOTSTRAPPING


    







int main(int argc, char **argv) {
  try {
    Settings settings; // parse the command-line argument list from the shell
    settings.ParseArgs(argc, argv);

    // Read the input tomogram
    cerr << "Reading tomogram \""<<settings.in_file_name<<"\"" << endl;
    MrcSimple tomo_in;
    tomo_in.Read(settings.in_file_name, false);
    // (Note: You can also use "tomo_in.Read(cin);" or "cin >> tomo;")
    tomo_in.PrintStats(cerr);      //Optional (display the tomogram size & format)

    // ---- mask ----

    // Optional: if there is a "mask", read that too
    MrcSimple mask;
    if (settings.mask_file_name != "") {
      cerr << "Reading mask \""<<settings.mask_file_name<<"\"" << endl;
      mask.Read(settings.mask_file_name, false);
      if ((mask.header.nvoxels[0] != tomo_in.header.nvoxels[0]) ||
          (mask.header.nvoxels[1] != tomo_in.header.nvoxels[1]) ||
          (mask.header.nvoxels[2] != tomo_in.header.nvoxels[2]))
        throw InputErr("Error: The size of the mask does not match the size of the tomogram.\n");
      // The mask should be 1 everywhere we want to consider, and 0 elsewhere.
      if (settings.use_mask_select) {
        for (int iz=0; iz<mask.header.nvoxels[2]; iz++)
          for (int iy=0; iy<mask.header.nvoxels[1]; iy++)
            for (int ix=0; ix<mask.header.nvoxels[0]; ix++)
              if (mask.aaafI[iz][iy][ix] == settings.mask_select)
                mask.aaafI[iz][iy][ix] = 1.0;
              else
                mask.aaafI[iz][iy][ix] = 0.0;
      }
    }

    if (settings.rescale01_in)
      tomo_in.Rescale01(mask.aaafI);

    // ---- make an array that will store the new tomogram we will create ----

    cerr << "allocating space for new tomogram..." << endl;
    MrcSimple tomo_out = tomo_in; //this will take care of allocating the array

    float voxel_width[3] = {1.0, 1.0, 1.0};

    // ---- Voxel size? ----

    if (settings.voxel_width > 0.0) {
      // Did the user manually specify the width of each voxel?
      voxel_width[0] = settings.voxel_width;
      voxel_width[1] = settings.voxel_width;
      voxel_width[2] = settings.voxel_width;
    }
    else {
      // Otherwise, infer it from the header of the MRC file
      voxel_width[0] = tomo_in.header.cellA[0]/tomo_in.header.nvoxels[0];
      voxel_width[1] = tomo_in.header.cellA[1]/tomo_in.header.nvoxels[1];
      voxel_width[2] = tomo_in.header.cellA[2]/tomo_in.header.nvoxels[2];
      if (settings.voxel_width_divide_by_10) {
        voxel_width[0] *= 0.1;
        voxel_width[1] *= 0.1;
        voxel_width[2] *= 0.1;
      }
      cerr << "voxel width in physical units = ("
           << voxel_width[0] << ", "
           << voxel_width[1] << ", "
           << voxel_width[2] << ")\n";
    }

    if ((voxel_width[0] <= 0.0) ||
        (voxel_width[1] <= 0.0) ||
        (voxel_width[2] <= 0.0))
      throw InputErr("Error in tomogram header: Invalid voxel width(s).\n"
                     "Use the -w argument to specify the voxel width.");

    if (abs((voxel_width[0] - voxel_width[1])
            /
            (0.5*(voxel_width[0] + voxel_width[1]))) > 0.0001)
      throw InputErr("Error in tomogram header: Unequal voxel widths in the x and y directions.\n"
                     "Use the -w argument to specify the voxel width.");
    for (int d=0; d<3; d++) {
      settings.width_a[d] /= voxel_width[d];
      settings.width_b[d] /= voxel_width[d];
      settings.dogsf_width[d] /= voxel_width[d];
      #ifndef DISABLE_BOOTSTRAPPING
      settings.f_bs_scramble_radius[d] /= voxel_width[d];
      settings.bs_scramble_radius[d] = ceil(settings.bs_scramble_radius[d]);
      #endif
      #ifndef DISABLE_TEMPLATE_MATCHING
      settings.template_background_radius[d] /= voxel_width[d];
      #endif
      //Commenting out:  The user cannot set filter_truncate_halfwidth directly anymore:
      //settings.filter_truncate_halfwidth[d] = floor(settings.filter_truncate_ratio * 
      //                                              MAX(settings.width_a[d],
      //                                                  settings.width_b[d]));
    }

    // At some point I was trying to be as general as possible and allowed
    // for the possibility that voxels need not be cubes (same width x,y,z)
    // Now, I realize that allowing for this possibility would slow
    // down the next step considerably, so I just assume cube-shaped voxels:
    assert((voxel_width[0] == voxel_width[1]) &&
           (voxel_width[1] == voxel_width[2]));
    for (int ir = 0; ir < settings.blob_widths.size(); ir++)
      settings.blob_widths[ir] /= voxel_width[0];

    settings.sphere_decals_radius /= voxel_width[0];
    settings.sphere_decals_shell_thickness /= voxel_width[0];



    // ---- filtering ----

    //cerr << "applying filter (window size in voxels: "
    //     << 1 + 2*settings.filter_truncate_halfwidth[0] << ","
    //     << 1 + 2*settings.filter_truncate_halfwidth[1] << ","
    //     << 1 + 2*settings.filter_truncate_halfwidth[2] << ")"
    //     << " ..." << endl;

    if (settings.filter_type == Settings::NONE) {
      cerr << "filter_type = Intensity Map <No convolution filter specified>\n";
      // Not needed:
      //for (int iz = 0; iz < size[2]; iz++)
      //  for (int iy = 0; iy < size[1]; iy++)
      //    for (int ix = 0; ix < size[0]; ix++)
      //      tomo_out.aaafI[iz][iy][ix]=tomo_in.aaafI[iz][iy][ix];
      // (We have copied the contents from tomo_in into tomo_out already.)
    } 

    else if (settings.filter_type == Settings::GAUSS) {

      HandleGauss(settings, tomo_in, tomo_out, mask, voxel_width);
        
    } // if (settings.filter_type == Settings::GAUSS)


    else if (settings.filter_type == Settings::GGAUSS) {

      HandleGGauss(settings, tomo_in, tomo_out, mask, voxel_width);

    } //else if (settings.filter_type == Settings::GGAUSS)




    else if (settings.filter_type == Settings::DOG) {

      HandleDog(settings, tomo_in, tomo_out, mask, voxel_width);

    } //if (settings.filter_type == Settings::DOG)



    else if (settings.filter_type == Settings::DOGG) {

      HandleDogg(settings, tomo_in, tomo_out, mask, voxel_width);

    } //if (settings.filter_type == Settings::DOGG)



    else if (settings.filter_type == Settings::DOGGXY) {

      HandleDoggXY(settings, tomo_in, tomo_out, mask, voxel_width);

    } //else if (settings.filter_type = Settings::DOGGXY)


    else if (settings.filter_type == Settings::DOG_SCALE_FREE) {

      HandleDogScaleFree(settings, tomo_in, tomo_out, mask, voxel_width);

    } //if (settings.filter_type == Settings::DOG_SCALE_FREE)



    else if (settings.filter_type == Settings::BLOB) {

      HandleBlobDetector(settings, tomo_in, tomo_out, mask, voxel_width);

    } //if (settings.filter_type == Settings::DOG)


    else if (settings.filter_type == Settings::TEMPLATE_GGAUSS) {

      #ifndef DISABLE_TEMPLATE_MATCHING
      HandleTemplateGGauss(settings, tomo_in, tomo_out, mask, voxel_width);
      #endif //#ifndef DISABLE_TEMPLATE_MATCHING

    } //else if (settings.filter_type == Settings::TEMPLATE_GGAUSS)


    else if (settings.filter_type == Settings::LOCAL_FLUCTUATIONS) {

      HandleLocalFluctuations(settings, tomo_in, tomo_out, mask, voxel_width);

    } //else if (settings.filter_type == Settings::TEMPLATE_GGAUSS)




    else if (settings.filter_type == Settings::TEMPLATE_GAUSS) {

      #ifndef DISABLE_TEMPLATE_MATCHING
      HandleTemplateGauss(settings, tomo_in, tomo_out, mask, voxel_width);
      #endif //#ifndef DISABLE_TEMPLATE_MATCHING

    } //else if (settings.filter_type == Settings::TEMPLATE_GAUSS)




    else if (settings.filter_type == Settings::BOOTSTRAP_DOGG) {

      #ifdef DISABLE_BOOSTRAPPING
      HandleBootstrappDogg(settings, tomo_in, tomo_out, mask, voxel_width);
      #endif //#ifndef DISABLE_BOOTSTRAPPING

    } //else if (settings.filter_type == Settings::BOOTSTRAP_DOGG) {


    // ----- distance_filter -----

    else if (settings.filter_type == settings.MIN_DISTANCE) {

      HandleMinDistance(settings, tomo_in, tomo_out, mask, voxel_width);

    }

    // ----- sphere_decals_filter -----

    else if (settings.filter_type == settings.SPHERE_DECALS) {

      HandleSphereDecals(settings, tomo_in, tomo_out, mask, voxel_width);

    }


    else {
      assert(false);  //should be one of the choices above
    }





    // --- Exchange light voxels for dark voxels ? ---

    if (settings.invert_output)
      tomo_out.Invert(mask.aaafI);




    // ----- thresholding: -----
    
    if (settings.use_thresholds) {

      HandleThresholds(settings, tomo_in, tomo_out, mask, voxel_width);

    }

    // ----- find local minima or maxima -----

    if (settings.find_minima || settings.find_maxima) {

      HandleExtrema(settings, tomo_in, tomo_out, mask, voxel_width);

    }





    // ----- masking: -----
    // Also, after thresholding, check again to see if the mask is zero
    // at this location.  If so, make sure that aaafI is 0 there
    // (or whatever the user has requested us to put there).
    // Do this after thresholding and inverting.
    if ((mask.aaafI) && (settings.use_mask_out))
      for (int iz=0; iz<mask.header.nvoxels[2]; iz++)
        for (int iy=0; iy<mask.header.nvoxels[1]; iy++)
          for (int ix=0; ix<mask.header.nvoxels[0]; ix++)
            if (mask.aaafI[iz][iy][ix] == 0.0)
              tomo_out.aaafI[iz][iy][ix] = settings.mask_out;

    tomo_out.FindMinMaxMean();


    // --- Rescale so that the lowest, highest voxels have density 0 and 1? ---

    if (settings.rescale01_out)
      tomo_out.Rescale01(mask.aaafI);








    // ------ Write the file ------
    if (settings.out_file_name != "") {
      cerr << "writing tomogram (in 32-bit float mode)" << endl;
      tomo_out.Write(settings.out_file_name);
      //(You can also use "tomo_out.Write(cout);" or "cout<<tomo_out;")
    }

  } // try {

  catch (InputErr& e) {
    cerr << "\n" << e.what() << endl;
    exit(1);
  }
} //main()


