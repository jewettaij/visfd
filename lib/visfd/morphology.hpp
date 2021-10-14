///   @file morphology.hpp
///   @brief  Functions relevant to image morphology.
///   @author Andrew Jewett
///   @date 2021-7-07

#ifndef _MORPHOLOGY_HPP
#define _MORPHOLOGY_HPP


#include <cstring>
#include <cassert>
#include <cmath>
#include <limits>
#include <ostream>
#include <vector>
using namespace std;
#include <err_visfd.hpp> // defines the "VisfdErr" exception type
#include <eigen3_simple.hpp>  // defines matrix diagonalizer (DiagonalizeSym3())
#include <visfd_utils.hpp>    // defines invert_permutation(), FindSpheres()
#include <alloc2d.hpp>    // defines Alloc2D() and Dealloc2D()
#include <alloc3d.hpp>    // defines Alloc3D() and Dealloc3D()
#include <filter1d.hpp>   // defines "Filter1D" (used in ApplySeparable())
#include <filter3d.hpp>   // defines common 3D image filters
#include <morphology_implementation.hpp>




namespace visfd {



/// @brief   Find all of the local minima in an image (aaafI).
///          The locations of minima are stored in the
///          *pv_minima_crds, and sorted in increasing order.
///          (IE they are sorted so that the most significant local minima
///          will appear first in the list.)
///          The caller can automatically discard minima which
///          are not sufficiently low, by supplying a "threshold" argument.
///          The optional aaafMask array (if not nullptr) can be used to ignore
///          certain voxels in the image (whose aaafMask entries are zero).
/// @note    If the caller wants to create an image showing where the
///          local minima are located, pass a non-NULL aaaiDest[][][] argument.
/// @note    Local minima on the boundaries of the image
///          (or near the edge of the mask)
///          are not as trustworthy since some of the neighboring voxels
///          will not be available for comparison.  These minima
///          can be ignored by setting allow_borders=false.  The number of
///          neighbors around every voxel which are considered (eg, 6, 18, 26)
///          can be controlled using the "connectivity" argument.
/// @return  Returns the number of minima found.

template<typename Scalar, typename Coordinate, typename IntegerIndex, typename Label>

size_t
FindMinima(int const image_size[3],            //!< size of input image array
           Scalar const *const *const *aaafI,   //!< input image array
           Scalar const *const *const *aaafMask, //!< if not nullptr then zero entries indicate which voxels to ignore
           vector<array<Coordinate, 3> > &minima_crds, //!< store the location of each minima
           vector<Scalar> &minima_scores, //!< store the brightness of each minima
           vector<IntegerIndex> &minima_nvoxels, //!< store number of voxels in each minima (usually 1)
           Scalar threshold=std::numeric_limits<Scalar>::infinity(), // OPTIONAL: Ignore minima which are not sufficiently low or high
           int connectivity=3,       //!< OPTIONAL: square root of search radius around each voxel (1=nearest_neighbors, 2=diagonal2D, 3=diagonal3D)
           bool allow_borders=true,  //!< OPTIONAL: if true, plateaus that touch the image border (or mask boundary) are valid minima
           Label ***aaaiDest=nullptr,  //!< OPTIONAL: create an image showing where the minima are?
           ostream *pReportProgress=nullptr)  //!< OPTIONAL: print progress to the user?
{
  return _FindExtrema(image_size[3],
                      aaafI,
                      aaafMask,
                      minima_crds,
                      minima_scores,
                      minima_nvoxels,
                      true,    // search for minima
                      threshold,
                      connectivity,
                      allow_borders,
                      aaaiDest,
                      pReportProgress);
} // FindMinima()





/// @brief   This function finds local maxima in the image,
///          but is otherwise identical to FindMinima().
///          See the documentation for FindMinima() for more information.
/// @return  Returns the number of maxima found.

template<typename Scalar, typename Coordinate, typename IntegerIndex, typename Label>

size_t
FindMaxima(int const image_size[3],            //!< size of input image array
           Scalar const *const *const *aaafI,   //!< input image array
           Scalar const *const *const *aaafMask, //!< if not nullptr then zero entries indicate which voxels to ignore
           vector<array<Coordinate, 3> > &maxima_crds, //!< store the location of each maxima
           vector<Scalar> &maxima_scores, //!< store the brightness of each maxima
           vector<IntegerIndex> &maxima_nvoxels, //!< store number of voxels in each maximama (usually 1)
           Scalar threshold=std::numeric_limits<Scalar>::infinity(), // OPTIONAL: Ignore maxima which are not sufficiently low or high
           int connectivity=3,       //!< OPTIONAL: square root of search radius around each voxel (1=nearest_neighbors, 2=diagonal2D, 3=diagonal3D)
           bool allow_borders=true,  //!< OPTIONAL: if true, plateaus that touch the image border (or mask boundary) are valid maxima
           Label ***aaaiDest=nullptr,  //!< OPTIONAL: create an image showing where the maxima are?
           ostream *pReportProgress=nullptr)  //!< OPTIONAL: print progress to the user?
{
  return _FindExtrema(image_size[3],
                      aaafI,
                      aaafMask,
                      maxima_crds,
                      maxima_scores,
                      maxima_nvoxels,
                      false,    // search for maxima
                      threshold,
                      connectivity,
                      allow_borders,
                      aaaiDest,
                      pReportProgress);
} // FindMaxima()




/// @brief Computes the grayscale dilation of an image with an arbitrary
///        structure factor, as defined here:
///       https://en.wikipedia.org/wiki/Dilation_(morphology)#Grayscale_dilation
///        The structure factor is implemented as a vector of tuples.
///        Each tuple contains the voxel coordinates and brightness value ("b")
///        for the structuer factor at that location.
///        (A version of this function using a spherical structure factors
///        is defined elsewhere.)
/// @note  This function has not been optimized for speed.
template<typename Scalar>
void
Dilate(vector<tuple<int,int,int,Scalar> > structure_factor, // a list of (ix,iy,iz,b) values, one for each voxel
       const int image_size[3],                //!< size of the image in x,y,z directions
       Scalar const *const *const *aaafSource, //!< source image array
       Scalar ***aaafDest,                     //!< filter results stored here
       Scalar const *const *const *aaafMask = nullptr,   //!< OPTIONAL: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
       ostream *pReportProgress = nullptr      //!< OPTIONAL: print progress to the user?
       )
{
  for (int iz=0; iz < image_size[2]; iz++) {
    if (pReportProgress)
      *pReportProgress << "  z = " << iz+1 << " of " << image_size[2] << endl;
    #pragma omp parallel for collapse(2)
    for (int iy=0; iy < image_size[1]; iy++) {
      for (int ix=0; ix < image_size[0]; ix++) {
        if ((aaafMask) && (aaafMask[iz][iy][ix] == 0.0))
          continue;
        Scalar max_f_plus_b = -std::numeric_limits<Scalar>::infinity();
        for (auto pVoxel = structure_factor.begin();
             pVoxel != structure_factor.end();
             pVoxel++) {
          int jx = std::get<0>(*pVoxel);
          int jy = std::get<1>(*pVoxel);
          int jz = std::get<2>(*pVoxel);
          int Ix = ix + jx;
          int Iy = iy + jy;
          int Iz = iz + jz;
          if ((Ix < 0) || (Ix >= image_size[0]) ||
              (Iy < 0) || (Iy >= image_size[1]) ||
              (Iz < 0) || (Iz >= image_size[2]))
            continue;
          if ((aaafMask) && (aaafMask[Iz][Iy][Ix] == 0.0))
            continue;
          Scalar b = std::get<3>(*pVoxel);
          Scalar f = aaafSource[Iz][Iy][Ix];
          max_f_plus_b = std::max(max_f_plus_b, f + b);
        }
        aaafDest[iz][iy][ix] = max_f_plus_b;
      }
    }
  }
} // Dilate(vector<tuple<int,int,int,Scalar> > structure_factor,...)





/// @brief Computes the grayscale erosion of an image with an arbitrary
///        structure factor, as defined here:
///        https://en.wikipedia.org/wiki/Erosion_(morphology)#Grayscale_erosion
///        The structure factor is implemented as a vector of tuples.
///        Each tuple contains the voxel coordinates and brightness value ("b")
///        for the structuer factor at that location.
///        (A version of this function using a spherical structure factors
///        is defined elsewhere.)
/// @note  This function has not been optimized for speed.
template<typename Scalar>
void
Erode(vector<tuple<int,int,int,Scalar> > structure_factor, // a list of (ix,iy,iz,b) values, one for each voxel
      const int image_size[3],                //!< size of the image in x,y,z directions
      Scalar const *const *const *aaafSource, //!< source image array
      Scalar ***aaafDest,                     //!< filter results stored here
      Scalar const *const *const *aaafMask = nullptr,   //!< OPTIONAL: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
      ostream *pReportProgress = nullptr      //!< OPTIONAL: print progress to the user?
      )
{
  for (int iz=0; iz < image_size[2]; iz++) {
    if (pReportProgress)
      *pReportProgress << "  z = " << iz+1 << " of " << image_size[2] << endl;
    #pragma omp parallel for collapse(2)
    for (int iy=0; iy < image_size[1]; iy++) {
      for (int ix=0; ix < image_size[0]; ix++) {
        if ((aaafMask) && (aaafMask[iz][iy][ix] == 0.0))
          continue;
        Scalar min_f_minus_b = std::numeric_limits<Scalar>::infinity();
        for (auto pVoxel = structure_factor.begin();
             pVoxel != structure_factor.end();
             pVoxel++) {
          int jx = std::get<0>(*pVoxel);
          int jy = std::get<1>(*pVoxel);
          int jz = std::get<2>(*pVoxel);
          int Ix = ix + jx;
          int Iy = iy + jy;
          int Iz = iz + jz;
          if ((Ix < 0) || (Ix >= image_size[0]) ||
              (Iy < 0) || (Iy >= image_size[1]) ||
              (Iz < 0) || (Iz >= image_size[2]))
            continue;
          if ((aaafMask) && (aaafMask[Iz][Iy][Ix] == 0.0))
            continue;
          Scalar b = std::get<3>(*pVoxel);
          Scalar f = aaafSource[Iz][Iy][Ix];
          min_f_minus_b = std::min(min_f_minus_b, f - b);
        }
        aaafDest[iz][iy][ix] = min_f_minus_b;
      }
    }
  }
} // Erode(vector<tuple<int,int,int,Scalar> > structure_factor,...)





/// @brief Computes the grayscale dilation of an image with a flat sphere.
///       https://en.wikipedia.org/wiki/Dilation_(morphology)#Grayscale_dilation
template<typename Scalar>
void
DilateSphere(Scalar radius,           //!< radius of the sphere (in voxels)
             int const image_size[3], //!< size of the image in x,y,z directions
             Scalar const *const *const *aaafSource, //!< source image array
             Scalar ***aaafDest,         //!< filter results stored here
             Scalar const *const *const *aaafMask = nullptr,   //!< OPTIONAL: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
             Scalar radius_max = 0.0,  //!< OPTIONAL: smoothly vary between two radii?
             Scalar bmax = 0.0, //!< OPTIONAL: structure factor b varies between 0 and this value
             ostream *pReportProgress = nullptr //!< OPTIONAL: report progress to the user
             )
{
  // structure_factor = a vector of (ix,iy,iz,b) values, one for each voxel
  vector<tuple<int,int,int,Scalar> > structure_factor;
  int Ri = ceil(std::max(radius, radius_max));
  for (int iz = -Ri; iz <= Ri; iz++) {
    for (int iy = -Ri; iy <= Ri; iy++) {
      for (int ix = -Ri; ix <= Ri; ix++) {
        bool add_this_voxel = false;
        Scalar b = 0.0;
        if (bmax == 0.0) {
          Scalar r = sqrt(SQR(ix) + SQR(iy) + SQR(iz));
          if (r <= radius)
            add_this_voxel = true;
        }
        else if (radius_max > radius) {
          Scalar r = sqrt(SQR(ix) + SQR(iy) + SQR(iz));
          if (r <= radius)
            add_this_voxel = true;
          else if (r <= radius_max) {
            add_this_voxel = true;
            b = -(r - radius) / (radius_max - radius);
            b *= bmax;
          }
        }
        else {
          // Calculate the distance of all of the corners to the origin.  The
          // corners of the cube are located at the corners of voxel ix,iy,iz
          //   (ix +/- 0.5, iy +/- 0.5, iz +/- 0.5)
          // If all 8 corners are within the sphere, 
          //   then include the voxel and set b = 0.
          // If all 8 corners are outside the sphere, ignore the voxel.
          // If any of the 8 corners lie within the sphere, then include
          // this voxel in the structure factor, but assign it a b value
          // between -1 and 0.  (-1 if the voxel does not overlap significantly
          // with the sphere. 0 if the voxel mostly overlaps with the sphere.)
          // (There are multiple ways of doing that.
          //  Here, I am doing it a sloppy way, but it should be good enough.)
          Scalar r_max = -std::numeric_limits<Scalar>::infinity();
          Scalar r_min =  std::numeric_limits<Scalar>::infinity();
          for (int jz=0; jz <= 1; jz++) {
            for (int jy=0; jy <= 1; jy++) {
              for (int jx=0; jx <= 1 ; jx++) {
                Scalar r = sqrt(SQR(ix+jx-0.5)+SQR(iy+jy-0.5)+SQR(iz+jz-0.5));
                if (r < r_min)
                  r_min = r;
                if (r > r_max)
                  r_max = r;
              }
            }
          }
          if (r_max < radius)
            add_this_voxel = true;
          else if (r_min > radius)
            add_this_voxel = false;
          else {
            add_this_voxel = true;
            b = -(r_max - radius) / (r_max - r_min);
            b *= bmax;
          }
        } // if ((bmax != 0) && (radius_max > radius))
        if (add_this_voxel)
          structure_factor.push_back(tuple<int,int,int,Scalar>(ix,iy,iz,b));
      } // for (int ix=0; ix < image_size[0]; ix++)
    } // for (int iy=0; iy < image_size[1]; iy++)
  } // for (int iz=0; iz < image_size[2]; iz++)

  if (pReportProgress)
    *pReportProgress << "DilateSphere(r="<<radius
                     <<" voxels, rmax="<<radius_max<<") progress:" <<endl;

  Dilate(structure_factor,
         image_size,
         aaafSource,
         aaafDest,
         aaafMask,
         pReportProgress);

} // DilateSphere()




/// @brief Computes the grayscale dilation of an image with a flat sphere.
///       https://en.wikipedia.org/wiki/Dilation_(morphology)#Grayscale_dilation
template<typename Scalar>
void
ErodeSphere(Scalar radius,            //!< radius of the sphere (in voxels)
            int const image_size[3],  //!< size of the image in x,y,z directions
            Scalar const *const *const *aaafSource, //!< source image array
            Scalar ***aaafDest,         //!< filter results stored here
            Scalar const *const *const *aaafMask = nullptr,   //!< OPTIONAL: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
            Scalar radius_max = 0.0,  //!< OPTIONAL: smoothly vary between two radii?
            Scalar bmax = 0.0, //!< OPTIONAL: structure factor b varies between 0 and this value
            ostream *pReportProgress = nullptr //!< OPTIONAL: report progress to the user
            )
{
  // structure_factor = a vector of (ix,iy,iz,b) values, one for each voxel
  vector<tuple<int,int,int,Scalar> > structure_factor;
  int Ri = ceil(std::max(radius, radius_max));
  for (int iz = -Ri; iz <= Ri; iz++) {
    for (int iy = -Ri; iy <= Ri; iy++) {
      for (int ix = -Ri; ix <= Ri; ix++) {
        bool add_this_voxel = false;
        Scalar b = 0.0;
        if (bmax == 0.0) {
          Scalar r = sqrt(SQR(ix) + SQR(iy) + SQR(iz));
          if (r <= radius)
            add_this_voxel = true;
        }
        else if (radius_max > radius) {
          Scalar r = sqrt(SQR(ix) + SQR(iy) + SQR(iz));
          if (r <= radius)
            add_this_voxel = true;
          else if (r <= radius_max) {
            add_this_voxel = true;
            b = -(r - radius) / (radius_max - radius);
            b *= bmax;
          }
        }
        else {
          // Calculate the distance of all of the corners to the origin.  The
          // corners of the cube are located at the corners of voxel ix,iy,iz
          //   (ix +/- 0.5, iy +/- 0.5, iz +/- 0.5)
          // If all 8 corners are within the sphere, 
          //   then include the voxel and set b = 0.
          // If all 8 corners are outside the sphere, ignore the voxel.
          // If any of the 8 corners lie within the sphere, then include
          // this voxel in the structure factor, but assign it a b value
          // between -1 and 0.  (-1 if the voxel does not overlap significantly
          // with the sphere. 0 if the voxel mostly overlaps with the sphere.)
          // (There are multiple ways of doing that.
          //  Here, I am doing it a sloppy way, but it should be good enough.)
          Scalar r_max = -std::numeric_limits<Scalar>::infinity();
          Scalar r_min =  std::numeric_limits<Scalar>::infinity();
          for (int jz=0; jz <= 1; jz++) {
            for (int jy=0; jy <= 1; jy++) {
              for (int jx=0; jx <= 1 ; jx++) {
                Scalar r = sqrt(SQR(ix+jx-0.5)+SQR(iy+jy-0.5)+SQR(iz+jz-0.5));
                if (r < r_min)
                  r_min = r;
                if (r > r_max)
                  r_max = r;
              }
            }
          }
          if (r_max < radius)
            add_this_voxel = true;
          else if (r_min > radius)
            add_this_voxel = false;
          else {
            add_this_voxel = true;
            b = -(r_max - radius) / (r_max - r_min);
            b *= bmax;
          }
        } // if ((bmax != 0) && (radius_max > radius))
        if (add_this_voxel)
          structure_factor.push_back(tuple<int,int,int,Scalar>(ix,iy,iz,b));
      } // for (int ix=0; ix < image_size[0]; ix++)
    } // for (int iy=0; iy < image_size[1]; iy++)
  } // for (int iz=0; iz < image_size[2]; iz++)

  if (pReportProgress)
    *pReportProgress << "ErodeSphere(r="<<radius
                     <<" voxels, rmax="<<radius_max<<") progress:" <<endl;

  Erode(structure_factor,
        image_size,
        aaafSource,
        aaafDest,
        aaafMask,
        pReportProgress);

} // ErodeSphere()




/// @brief Computes the grayscale opening of an image with a flat sphere.
///        https://en.wikipedia.org/wiki/Opening_(morphology)
template<typename Scalar>
void
OpenSphere(Scalar radius,             //!< radius of the sphere (in voxels)
           const int image_size[3],   //!< size of the image in x,y,z directions
           Scalar const *const *const *aaafSource, //!< source image array
           Scalar ***aaafDest,        //!< filter results stored here
           Scalar const *const *const *aaafMask = nullptr,   //!< OPTIONAL: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
           Scalar radius_max = 0.0,  //!< OPTIONAL: smoothly vary between two radii?
           Scalar bmax = 0.0, //!< OPTIONAL: structure factor b varies between 0 and this value
           ostream *pReportProgress = nullptr //!< OPTIONAL: report progress to the user
           )
{
  float ***aaafTmp;  //temporary space to store image after erosion

  aaafTmp = Alloc3D<Scalar>(image_size);

  ErodeSphere(radius,
              image_size,
              aaafSource,
              aaafTmp,  //<--save the resulting image in aaafTmp
              aaafMask,
              radius_max,
              bmax,
              pReportProgress);

  DilateSphere(radius,
               image_size,
               aaafTmp,  //<--use aaafTmp as the source image
               aaafDest, //<--save the resulting image in aaafDest
               aaafMask,
               radius_max,
               bmax,
               pReportProgress);

  Dealloc3D(aaafTmp);

} // OpenSphere()




/// @brief Computes the grayscale closing of an image with a flat sphere.
///        https://en.wikipedia.org/wiki/Closing_(morphology)
template<typename Scalar>
void
CloseSphere(Scalar radius,            //!< radius of the sphere (in voxels)
            const int image_size[3],  //!< size of the image in x,y,z directions
            Scalar const *const *const *aaafSource, //!< source image array
            Scalar ***aaafDest,         //!< filter results stored here
            Scalar const *const *const *aaafMask = nullptr,   //!< OPTIONAL: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
            Scalar radius_max = 0.0,  //!< OPTIONAL: smoothly vary between two radii?
            Scalar bmax = 0.0, //!< OPTIONAL: structure factor b varies between 0 and this value
            ostream *pReportProgress = nullptr //!< OPTIONAL: report progress to the user
            )
{
  float ***aaafTmp;  //temporary space to store image after dilation

  aaafTmp = Alloc3D<Scalar>(image_size);

  DilateSphere(radius,
               image_size,
               aaafSource,
               aaafTmp,  //<--save the resulting image in aaafTmp
               aaafMask,
               radius_max,
               bmax,
               pReportProgress);

  ErodeSphere(radius,
              image_size,
              aaafTmp,  //<--use aaafTmp as the source image
              aaafDest, //<--save the resulting image in aaafDest
              aaafMask,
              radius_max,
              bmax,
              pReportProgress);

  Dealloc3D(aaafTmp);

} // CloseSphere()



/// @brief Computes the grayscale white top-hat transform of an image
///        using a flat spherical structure factor.
///        https://en.wikipedia.org/wiki/Top-hat_transform
template<typename Scalar>
void
WhiteTopHatSphere(Scalar radius,             //!< radius of the sphere (in voxels)
                  const int image_size[3],   //!< size of the image in x,y,z directions
                  Scalar const *const *const *aaafSource, //!<source image array
                  Scalar ***aaafDest,         //!< filter results stored here
                  Scalar const *const *const *aaafMask = nullptr,   //!< OPTIONAL: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
                  Scalar radius_max = 0.0,  //!< OPTIONAL: smoothly vary between two radii?
                  Scalar bmax = 0.0, //!< OPTIONAL: structure factor b varies between 0 and this value
                  ostream *pReportProgress = nullptr //!< OPTIONAL: report progress to the user
                  )
{
  float ***aaafTmp;  //temporary space to store image after erosion

  aaafTmp = Alloc3D<Scalar>(image_size);

  // Compute the image opening
  OpenSphere(radius,
             image_size,
             aaafSource,
             aaafTmp,  //<--save the resulting image in aaafTmp
             aaafMask,
             radius_max,
             bmax,
             pReportProgress);

  // Subtract it from the original image brightness
  for (int iz = 0; iz < image_size[2]; iz++)
    for (int iy = 0; iy < image_size[1]; iy++)
      for (int ix = 0; ix < image_size[0]; ix++)
        aaafDest[iz][iy][ix] -= aaafTmp[iz][iy][ix];

  Dealloc3D(aaafTmp);

} // WhiteTopHatSphere()



/// @brief Computes the grayscale black top-hat transform of an image
///        using a flat spherical structure factor.
///        https://en.wikipedia.org/wiki/Top-hat_transform
template<typename Scalar>
void
BlackTopHatSphere(Scalar radius,             //!< radius of the sphere (in voxels)
                  const int image_size[3],   //!< size of the image in x,y,z directions
                  Scalar const *const *const *aaafSource,//!< source image array
                  Scalar ***aaafDest,         //!< filter results stored here
                  Scalar const *const *const *aaafMask = nullptr,   //!< OPTIONAL: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
                  Scalar radius_max = 0.0,  //!< OPTIONAL: smoothly vary between two radii?
                  Scalar bmax = 0.0, //!< OPTIONAL: structure factor b varies between 0 and this value
                  ostream *pReportProgress = nullptr //!< OPTIONAL: report progress to the user
                  )
{
  float ***aaafTmp;  //temporary space to store image after erosion

  aaafTmp = Alloc3D<Scalar>(image_size);

  // Compute the image closing
  CloseSphere(radius,
              image_size,
              aaafSource,
              aaafTmp,  //<--save the resulting image in aaafTmp
              aaafMask,
              radius_max,
              bmax,
              pReportProgress);

  // Subtract the original image brightness
  for (int iz = 0; iz < image_size[2]; iz++)
    for (int iy = 0; iy < image_size[1]; iy++)
      for (int ix = 0; ix < image_size[0]; ix++)
        aaafDest[iz][iy][ix] = aaafTmp[iz][iy][ix] - aaafDest[iz][iy][ix];

  Dealloc3D(aaafTmp);

} // BlackTopHatSphere()




/// @brief Computes the grayscale dilation of an image with a flat sphere.
///       https://en.wikipedia.org/wiki/Dilation_(morphology)#Grayscale_dilation
template<typename Scalar>
void
DilateSphere(Scalar radius,           //!< radius of the sphere (in voxels)
             int const image_size[3], //!< size of the image in x,y,z directions
             Scalar const *const *const *aaafSource, //!< source image array
             Scalar ***aaafDest,         //!< filter results stored here
             Scalar const *const *const *aaafMask = nullptr,   //!< OPTIONAL: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
             Scalar radius_max = 0.0,  //!< OPTIONAL: smoothly vary between two radii?
             Scalar bmax = 0.0, //!< OPTIONAL: structure factor b varies between 0 and this value
             ostream *pReportProgress = nullptr //!< OPTIONAL: report progress to the user
             )

/// @brief Computes either the grayscale erosion or dilation of an image with an
///        ellipsoidal-shaped structure factor.  For more details, see:
///       https://en.wikipedia.org/wiki/Dilation_(morphology)#Grayscale_dilation
///       https://en.wikipedia.org/wiki/Erosion_(morphology)#Grayscale_erosion
///        The size of the ellipsoid to search over (principal axis lengths) 
///        is given by the "search_ellipsoid" argument.
///        The directions to search over are stored in "aaaafDirections"
///        aaaafTensor[iz][iy][ix][] is a 1-dimensional array which
///        stores the direction of the ellipsoid at voxel position ix,iy,iz.
///        This 1-dimensional array stores 6 numbers
///        aaaafTensor[iz][iy][ix][0..5]
///        which are the 6 non-redundant entries of a symmetric tensor whose
///        eigenvectors correspond the directions of the ellipsoid at this
///        location. (They should be sorted by eigenvalue in descending order.)
///        Since these directions are not algined with the X,Y,Z axis,
///        and they are changing throughout the image, the distance that we
///        travel in each of these principal directions is determined by the
///        search_iters[] array argument, and is typically less than one voxel.
///        This ensures that we visit every real nearby voxel at least once,
///        even though we may sometimes visit the same voxel multiple times.
///        
/// @note  This function has not been optimized for speed.
template<typename Scalar>
void
DilateOrErodeAniso(bool dilate,                            //!< dilate or erode?
                   Scalar search_radius[3],                //!< size of the eillipsoidal structure factor in voxels
                   const int image_size[3],                //!< size of the image in x,y,z directions
                   Scalar const *const *const *aaafSource, //!< source image array
                   Scalar ***aaafDest,                     //!< filter results stored here
                   TensorContainer const *const *const *aaaafTensor, //!< ellipsoid orientations stored here (see above)
                   Scalar const *const *const *aaafMask = nullptr,   //!< OPTIONAL: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
                   const int *search_radius_i=nullptr,     //!< iterate this many times in each direction during the search
                   ostream *pReportProgress=nullptr        //!< OPTIONAL: print progress to the user?
       )
{
  assert(aaaafDirections);
  // figure out the size of the nearby region we will search over
  // (ie the size of the "structure factor")
  int niters[3];
  Scalar incr[3];
  if (search_radius_i == nullptr) {
    for (int d=0; d < 3; d++) {
      assert(search_radius[d] >= 0);
      niters[d] = std::floor(search_radius[d]+0.5);
      incr[d] = 1.0;
    }
  }
  else {
    for (int d=0; d < 3; d++) {
      assert(search_radius_i[d] >= 0);
      niters[d] = search_radius_i[d];
      nincr[d] = 0.0;
      if (search_radius_i[d] > 1)
        incr[d] = search_radius_i[d] / (search_radius_i[d] - 1);
    }
  }

  // Now, figure out which voxels lie within the ellipsoid.
  // We will save the relative coordinates of these voxels in an array.
  // We will lookup the voxel indices from this array.
  CONTINUEHERE;

  // Now figure out the size of the aaaaaaDirections array (in the Z direction)
  // It should equal the maximum width of the ellipsoid (from any direction).
  CONTINUEHERE;



  // now loop over the voxels in the source image
  for (int iz=0; iz < image_size[2]; iz++) {
    if (pReportProgress)
      *pReportProgress << "  z = " << iz+1 << " of " << image_size[2] << endl;
    vector<array<int,2> > voxels_to_consider;
    for (int iy=0; iy < image_size[1]; iy++) {
      for (int ix=0; ix < image_size[0]; ix++) {
        if ((aaafMask) && (aaafMask[iz][iy][ix] == 0.0))
          continue;
        else
          voxels_to_consider.push_back(array<int,2>(ix,iy));
      }
    }
    // Now precompute the aaaaaaDirections array for voxels near this z value
    CONTINUEHERE;




    #pragma omp parallel for collapse(1)
    for (size_t ip = 0; ip < voxels_to_consider.size(); ip++)
    {
      ix = voxels_to_consider[ip][0];
      iy = voxels_to_consider[ip][1];
      Scalar max_f_plus_b = -std::numeric_limits<Scalar>::infinity();

      //CONTINUEHERE
      Scalar dx, dy, dz;
      Scalar du = 0.0;
      Scalar dv = 0.0;
      Scalar dw = 0.0;
      for (int jw = 0; jw <= niters[2]; jw++) {
        for (int jv = 0; jv <= niters[1]; jv++) {
          for (int ju = 0; ju <= niters[0]; ju++) {
            int Ix = std::floor(x+0.5);
            int Iy = std::floor(y+0.5);
            int Iz = std::floor(z+0.5);
            if ((Ix < 0) || (Ix >= image_size[0]) ||
                (Iy < 0) || (Iy >= image_size[1]) ||
                (Iz < 0) || (Iz >= image_size[2]))
              continue;
            if ((aaafMask) && (aaafMask[Iz][Iy][Ix] == 0.0))
              continue;

            Scalar f = aaafSource[Iz][Iy][Ix];

            Scalar b = 0.0; // In this implementation, we assume b=0.  See:
            //https://en.wikipedia.org/wiki/Erosion_(morphology)#Grayscale_erosion
            //https://en.wikipedia.org/wiki/Dilation_(morphology)#Grayscale_dilation
            max_f_plus_b = std::max(max_f_plus_b, f + b);
            min_f_minus_b = std::max(min_f_minus_b, f - b);

            dx += incr[0] * directions[0][0];
            dy += incr[0] * directions[0][1];
            dz += incr[0] * directions[0][2];
          }
          dx += incr[1] * directions[1][0];
          dy += incr[1] * directions[1][1];
          dz += incr[1] * directions[1][2];
        }
        dx += incr[2] * directions[2][0];
        dy += incr[2] * directions[2][1];
        dz += incr[2] * directions[2][2];
      }
      if (dilate)
        aaafDest[iz][iy][ix] = max_f_plus_b;
      else
        aaafDest[iz][iy][ix] = min_f_minus_b;
    }
  }
} // DilateOrErodeAniso()




} //namespace visfd



#endif //#ifndef _MORPHOLOGY_HPP
