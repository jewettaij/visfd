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

  // This function is defined in visfd_utils_implementation.hpp
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


        


/// @brief Computes the grayscale dilation of an image with an arbitrary
///        structure factor, as defined here:
///       https://en.wikipedia.org/wiki/Dilation_(morphology)#Grayscale_dilation
///        The structure factor is implemented as a vector of tuples.
///        Each tuple contains the voxel coordinates and brightness value ("b")
///        for the structuer factor at that location.
///        (A version of this function using a spherical structure factors
///        is defined elsewhere.)
template<typename Scalar>
void
Dilate(vector<tuple<int,int,int,Scalar> > structure_factor, // a list of (ix,iy,iz,b) values, one for each voxel
       const int image_size[3],                //!< size of the image in x,y,z directions
       Scalar const *const *const *aaafSource, //!< image array aaafI[iz][iy][ix]
       Scalar ***aaafDest,                     //!< image array aaafI[iz][iy][ix]       
       Scalar const *const *const *aaafMask = nullptr,   //!< optional: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0

       ostream *pReportProgress = nullptr)
{
  for (int iz=0; iz < image_size[2]; iz++) {
    if (pReportProgress)
      *pReportProgress << "  z = " << iz+1 << " of " << image_size[2] << endl;
    #pragma omp parallel for collapse(2)
    for (int iy=0; iy < image_size[1]; iy++) {
      for (int ix=0; ix < image_size[0]; ix++) {
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
template<typename Scalar>
void
Erode(vector<tuple<int,int,int,Scalar> > structure_factor, // a list of (ix,iy,iz,b) values, one for each voxel
      const int image_size[3],                //!< size of the image in x,y,z directions
      Scalar const *const *const *aaafSource, //!< image array aaafI[iz][iy][ix]
      Scalar ***aaafDest,                     //!< image array aaafI[iz][iy][ix]       
      Scalar const *const *const *aaafMask = nullptr,   //!< optional: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0

      ostream *pReportProgress = nullptr)
{
  for (int iz=0; iz < image_size[2]; iz++) {
    if (pReportProgress)
      *pReportProgress << "  z = " << iz+1 << " of " << image_size[2] << endl;
    #pragma omp parallel for collapse(2)
    for (int iy=0; iy < image_size[1]; iy++) {
      for (int ix=0; ix < image_size[0]; ix++) {
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
///        The radius paramater can be a floating point number, but its
///        units are in voxels.
template<typename Scalar>
void
DilateSphere(Scalar radius,              //!< radius of the sphere
             int const image_size[3],    //!< size of the image in x,y,z directions
             Scalar const *const *const *aaafSource, //!< image array aaafI[iz][iy][ix]
             Scalar ***aaafDest,         //!< image array aaafI[iz][iy][ix]
             Scalar const *const *const *aaafMask = nullptr,   //!< optional: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
             bool smooth_boundary = false, //!< attempt to produce a rounder smoother structure factor that varies between 0 and -1 at its boundary voxels
             ostream *pReportProgress = nullptr //!< report progress to the user
             )
{
  // structure_factor = a vector of (ix,iy,iz,b) values, one for each voxel
  vector<tuple<int,int,int,Scalar> > structure_factor;

  int Ri = ceil(radius);
  for (int iz = -Ri; iz <= Ri; iz++) {
    for (int iy = -Ri; iy <= Ri; iy++) {
      for (int ix = -Ri; ix <= Ri; ix++) {
        bool add_this_voxel = false;
        Scalar b = 0.0;
        if (! smooth_boundary) {
          Scalar r = sqrt(SQR(ix) + SQR(iy) + SQR(iz));
          if (r <= radius) {
            add_this_voxel = true;
            b = 1.0;
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
          if (r_max < radius) {
            add_this_voxel = true;
            b = 1.0;
          }
          else if (r_min > radius) {
            add_this_voxel = false;
            b = 0.0;
          }
          else {
            add_this_voxel = true;
            b = (r_max - radius) / (r_max - r_min);
            b = -b;
          }
        } // if (smooth_boundary)
        if (add_this_voxel)
          structure_factor.push_back(tuple<int,int,int,Scalar>(ix,iy,iz,b));
      } // for (int ix=0; ix < image_size[0]; ix++)
    } // for (int iy=0; iy < image_size[1]; iy++)
  } // for (int iz=0; iz < image_size[2]; iz++)

  if (pReportProgress)
    *pReportProgress << "DilateSphere() progress:" << endl;

  Dilate(structure_factor,
         image_size,
         aaafSource,
         aaafDest,
         aaafMask,
         pReportProgress);

} // DilateSphere()




/// @brief Computes the grayscale dilation of an image with a flat sphere.
///       https://en.wikipedia.org/wiki/Dilation_(morphology)#Grayscale_dilation
///        The radius paramater can be a floating point number, but its
///        units are in voxels.
template<typename Scalar>
void
ErodeSphere(Scalar radius,              //!< radius of the sphere
            int const image_size[3],    //!< size of the image in x,y,z directions
            Scalar const *const *const *aaafSource, //!< image array aaafI[iz][iy][ix]
            Scalar ***aaafDest,         //!< image array aaafI[iz][iy][ix]
            Scalar const *const *const *aaafMask = nullptr,   //!< optional: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
            bool smooth_boundary = false, //!< attempt to produce a rounder smoother structure factor that varies between 0 and -1 at its boundary voxels
            ostream *pReportProgress = nullptr //!< report progress to the user
            )
{
  // structure_factor = a vector of (ix,iy,iz,b) values, one for each voxel
  vector<tuple<int,int,int,Scalar> > structure_factor;

  int Ri = ceil(radius);
  for (int iz = -Ri; iz <= Ri; iz++) {
    for (int iy = -Ri; iy <= Ri; iy++) {
      for (int ix = -Ri; ix <= Ri; ix++) {
        bool add_this_voxel = false;
        Scalar b = 0.0;
        if (! smooth_boundary) {
          Scalar r = sqrt(SQR(ix) + SQR(iy) + SQR(iz));
          if (r <= radius) {
            add_this_voxel = true;
            b = 1.0;
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
          if (r_max < radius) {
            add_this_voxel = true;
            b = 1.0;
          }
          else if (r_min > radius) {
            add_this_voxel = false;
            b = 0.0;
          }
          else {
            add_this_voxel = true;
            b = (r_max - radius) / (r_max - r_min);
          }
        } // if (smooth_boundary)
        if (add_this_voxel)
          structure_factor.push_back(tuple<int,int,int,Scalar>(ix,iy,iz,b));
      } // for (int ix=0; ix < image_size[0]; ix++)
    } // for (int iy=0; iy < image_size[1]; iy++)
  } // for (int iz=0; iz < image_size[2]; iz++)

  if (pReportProgress)
    *pReportProgress << "ErodeSphere() progress:" << endl;

  Erode(structure_factor,
        image_size,
        aaafSource,
        aaafDest,
        aaafMask,
        pReportProgress);

} // ErodeSphere()




/// @brief Computes the grayscale opening of an image with a flat sphere.
///        https://en.wikipedia.org/wiki/Opening_(morphology)
///        The radius paramater can be a floating point number, but its
///        units are in voxels.
template<typename Scalar>
void
OpenSphere(Scalar radius,             //!< radius of the sphere
           const int image_size[3],   //!< size of the image in x,y,z directions
           Scalar const *const *const *aaafSource, //!< image array aaafI[iz][iy][ix]
           Scalar ***aaafDest,        //!< image array aaafI[iz][iy][ix]
           Scalar const *const *const *aaafMask = nullptr,   //!< optional: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
           bool smooth_boundary = false, //!< attempt to produce a round and smooth structure factor that varies between 0 and -1 at its boundary voxels
           ostream *pReportProgress = nullptr //!< report progress to the user
           )
{
  float *afTmp;      //temporary space to store image after erosion
  float ***aaafTmp;  //temporary space to store image after erosion

  Alloc3D(image_size,
          &afTmp,
          &aaafTmp);

  ErodeSphere(radius,
              image_size,
              aaafSource,
              aaafTmp,  //<--save the resulting image in aaafTmp
              aaafMask,
              smooth_boundary,
              pReportProgress);

  DilateSphere(radius,
               image_size,
               aaafTmp,  //<--use aaafTmp as the source image
               aaafDest, //<--save the resulting image in aaafDest
               aaafMask,
               smooth_boundary,
               pReportProgress);

  Dealloc3D(image_size,
            &afTmp,
            &aaafTmp);

} // OpenSphere()




/// @brief Computes the grayscale closing of an image with a flat sphere.
///        https://en.wikipedia.org/wiki/Closing_(morphology)
///        The radius paramater can be a floating point number, but its
///        units are in voxels.
template<typename Scalar>
void
CloseSphere(Scalar radius,            //!< radius of the sphere
            const int image_size[3],  //!< size of the image in x,y,z directions
            Scalar const *const *const *aaafSource, //!< image array aaafI[iz][iy][ix]
            Scalar ***aaafDest,       //!< image array aaafI[iz][iy][ix]
            Scalar const *const *const *aaafMask = nullptr,   //!< optional: ignore voxel ix,iy,iz if aaafMask!=nullptr and aaafMask[iz][iy][ix]==0
            bool smooth_boundary = false, //!< attempt to produce a round and smooth structure factor that varies between 0 and -1 at its boundary voxels
            ostream *pReportProgress = nullptr //!< report progress to the user
            )
{
  float *afTmp;      //temporary space to store image after dilation
  float ***aaafTmp;  //temporary space to store image after dilation

  Alloc3D(image_size,
          &afTmp,
          &aaafTmp);

  DilateSphere(radius,
               image_size,
               aaafSource,
               aaafTmp,  //<--save the resulting image in aaafTmp
               aaafMask,
               smooth_boundary,
               pReportProgress);

  ErodeSphere(radius,
              image_size,
              aaafTmp,  //<--use aaafTmp as the source image
              aaafDest, //<--save the resulting image in aaafDest
              aaafMask,
              smooth_boundary,
              pReportProgress);

  Dealloc3D(image_size,
            &afTmp,
            &aaafTmp);

} // CloseSphere()



} //namespace visfd



#endif //#ifndef _MORPHOLOGY_HPP
