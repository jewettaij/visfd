///   @file draw.hpp
///   @brief a collection of functions for drawing objects
///          within an image for the purpose of annotation.
///   @author Andrew Jewett
///   @date 2019-4-15

#ifndef _DRAW_HPP
#define _DRAW_HPP

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
#include <visfd_utils.hpp>    // defines invert_permutation(), AveArray(), ...
#include <filter1d.hpp>   // defines "Filter1D" (used in ApplySeparable())
#include <filter3d.hpp>   // defines common 3D image filters



namespace visfd {



/// @brief SimpleRegion is a class used for storing
/// It stores coordinates describing the size and shape of the elementary
/// region being described, as well as a data member named "value".
/// The "value" is a number.  It is typically 1 or -1, but can be any number.
/// Typically the number is interpeted this way:
/// If the region you are describing includes
/// the values of the object (eg. rect or sphere),
/// then set value = 1.  If the region EXCLUDES
/// the values of this object, set value = -1.
/// You can build complex regions using multiple
/// "SimpleRegion" objects.
/// For example: You can described a rectangular region
/// with a spherical cavity by setting value=1 for
/// the rectangle, and value=-1 for the sphere.
template<typename Scalar>
struct SimpleRegion {
  struct Rect {
    Scalar xmin;
    Scalar xmax;
    Scalar ymin;
    Scalar ymax;
    Scalar zmin;
    Scalar zmax;
    //Rect() { xmin=0; xmax=-1; ymin=0; ymax=-1; zmin=0; zmax=-1; }
  };
  struct Sphere {
    Scalar x0;
    Scalar y0;
    Scalar z0;
    Scalar r;
    //Sphere() { x0=0; y0=0; z0=0; r=-1; }
  };
  enum RegionType {RECT, SPHERE};
  RegionType type;
  union {
    Rect rect;
    Sphere sphere;
  } data;
  Scalar value;        // A number associated with this region (see above).

  SimpleRegion():type(RECT), value(1) {
    //initialize with non-sensical values
    data.rect.xmin=0;
    data.rect.xmax=-1;
    data.rect.ymin=0;
    data.rect.ymax=-1;
    data.rect.zmin=0;
    data.rect.zmax=-1;
  }
}; // struct SimpleRegion






template<typename Scalar>
void
DrawRegions(int const image_size[3], //!< image size
            Scalar ***aaafDest,  //!< array where we should write new image
            Scalar const *const *const *aaafMask, //!< Optional: ignore voxels where mask==0
            vector<SimpleRegion<Scalar> > vRegions,        //!< a list of primitive regions and their brightnesses
            bool negative_means_subtract = false  //!< should we interpret negative brightness as set subtraction?
            )
{

  // Special case:  Should the background voxels be initialized with 1s?
  // If the first region has a negative value and negative_means_subtract==true,
  // then it means that the caller wants to remove (subtract) those voxels from
  // the set of non-zero background voxels, by setting their brightness to zero.
  // That doesn't make sense if all of the voxels in the image
  // already have brightness zero.  So I interpret this special case as a
  // request to subtract these voxels from a binary image which is filled with
  // non-zero values.  Since a binary image only contains 1s or 0s, I will
  // initialize the image full of 1s.  Then in the next step, we will remove
  // the voxels in the first region (vRegions[0]) from this background of 1s.
  if (negative_means_subtract &&
      (vRegions.size() > 0) && (vRegions[0].value < 0)) {
    bool all_zero = true;
    for (int iz = 0; iz < image_size[2] && all_zero; iz++) {
      for (int iy = 0; iy < image_size[1] && all_zero; iy++) {
        for (int ix = 0; ix < image_size[0] && all_zero; ix++) {
          if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
            continue;
          if (aaafDest[iz][iy][ix] != 0.0)
            all_zero = false;
        }
      }
    }
    if (all_zero) {
      for (int iz = 0; iz < image_size[2] && all_zero; iz++) {
        for (int iy = 0; iy < image_size[1] && all_zero; iy++) {
          for (int ix = 0; ix < image_size[0] && all_zero; ix++) {
            if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
              continue;
            aaafDest[iz][iy][ix] = 1.0;
          }
        }
      }
    }
  } // if (negative_means_subtract && (vRegions[0].value < 0))


  for (int i=0; i < vRegions.size(); i++) {

    switch (vRegions[i].type) {
    case SimpleRegion<Scalar>::SPHERE:
      {
        Scalar R = vRegions[i].data.sphere.r;
        int Ri = ceil(R-0.5);
        Scalar value = vRegions[i].value;
        int ix  = floor(vRegions[i].data.sphere.x0 + 0.5);
        int iy  = floor(vRegions[i].data.sphere.y0 + 0.5);
        int iz  = floor(vRegions[i].data.sphere.z0 + 0.5);
        for (int jz = -Ri; jz <= Ri; jz++) {
          for (int jy = -Ri; jy <= Ri; jy++) {
            Scalar descr = R*R - (jy*jy + jz*jz);
            if (descr < 0.0)
              continue;
            int xrange = floor(sqrt(descr));
            for (int jx = -xrange; jx <= xrange; jx++) {
              assert(jx*jx + jy*jy + jz*jz <= R*R);
              if ((ix+jx < 0) ||
                  (ix+jx >= image_size[0]) ||
                  (iy+jy < 0) ||
                  (iy+jy >= image_size[1]) ||
                  (iz+jz < 0) ||
                  (iz+jz >= image_size[2]))
                continue;
              else if (aaafMask && (aaafMask[iz+jz][iy+jy][ix+jx] == 0.0))
                continue;

              // now decide what to do with this voxel
              if (value < 0) {
                if (negative_means_subtract &&
                    (aaafDest[iz+jz][iy+jy][ix+jx] > 0))
                  // then exclude this voxel from the existing set of voxels
                  aaafDest[iz+jz][iy+jy][ix+jx] = 0.0; // 0=excluded, 1=included
              }
              else
                // then set this voxel's brightness directly
                aaafDest[iz+jz][iy+jy][ix+jx] = value;
            } //for (int jx = -xrange; jx <= xrange; jx++)
          } //for (int jy = -Ri; jy <= Ri; jy++)
        } //for (int jz = -Ri; jz <= Ri; jz++)
      } //case SimpleRegion<Scalar>::SPHERE:
      break;

    case SimpleRegion<Scalar>::RECT:
      {
        Scalar value = vRegions[i].value;
        Scalar ixmin = floor(vRegions[i].data.rect.xmin + 0.5);
        Scalar ixmax = floor(vRegions[i].data.rect.xmax + 0.5);
        Scalar iymin = floor(vRegions[i].data.rect.ymin + 0.5);
        Scalar iymax = floor(vRegions[i].data.rect.ymax + 0.5);
        Scalar izmin = floor(vRegions[i].data.rect.zmin + 0.5);
        Scalar izmax = floor(vRegions[i].data.rect.zmax + 0.5);

        for (int iz=std::max<float>(izmin, 0);
             iz <= std::min<float>(izmax, image_size[2]-1);
             iz++) {
          for (int iy=std::max<float>(iymin, 0);
               iy <= std::min<float>(iymax, image_size[1]-1);
               iy++) {
            for (int ix=std::max<float>(ixmin, 0);
                 ix <= std::min<float>(ixmax, image_size[0]-1);
                 ix++) {
              if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
                continue;
              // now decide what to do with this voxel
              if (value < 0) {
                if (negative_means_subtract &&
                    (aaafDest[iz][iy][ix] > 0))
                  // then exclude this voxel from the existing set of voxels
                  aaafDest[iz][iy][ix] = 0.0; // 0=excluded, 1=included
              }
              else
                // then set this voxel's brightness directly
                aaafDest[iz][iy][ix] = value;
            } // loop over ix
          } // loop over iy
        } // loop over iz
      } //case SimpleRegion<Scalar>::RECT:
      break;

    default:
      assert(false); //this line should not be reached
      break;

    } //switch (vRegions[i].type)

  } //for (int i=0; i < vRegions.size(); i++)
} //DrawRegions()





/// @brief  Create a 3D image containing multiple spheres (or spherical shells)
///         various sizes, thicknesses, and locations, specified by the caller.
///         The resulting spheres can (optionally) be superimposed with the
///         existing image (if the aaafBackground image array is != nullptr).

template<typename Scalar>

void
DrawSpheres(int const image_size[3], //!< image size
            Scalar ***aaafDest,  //!< array where we should write new image
            Scalar const *const *const *aaafMask,   //!< Optional: ignore voxels where mask==0
            const vector<array<Scalar,3> > &centers, //!< Coordinates for the center of each sphere (blob)
            const vector<Scalar> *pDiameters=nullptr,        //!< Optional: Diameter of each spherical shell (in voxels)
            const vector<Scalar> *pShellThicknesses=nullptr, //!< Optional: Thickness of each spherical shell (in voxels)
            const vector<Scalar> *pVoxelIntensitiesForeground=nullptr, //!< Optional: Assign voxels in spherical shell to this brightness (the vector should contain a separate entry for every sphere).
            Scalar const *const *const *aaafBackground = nullptr, //!< Optional: Specify the background voxels?
            Scalar voxel_intensity_background_offset = 0.0, //!< Optional: Add this number to the brightness of the background voxels, if specified.  (Otherwise set them to this number.)
            Scalar voxel_intensity_background_rescale = 1.0, //!< Optional: Multiply the background voxel brightnesses (if specified) by this number.
            bool voxel_intensity_background_normalize = false, //!< Optional: Shift and rescale background brightness automatically?
            bool voxel_intensity_foreground_normalize = false, //!< Optional: divide brightnesses by number of voxels in sphere (shell)? (rarely useful)
            ostream *pReportProgress = nullptr //!<Optional: report progress to the user?
            )
{
  assert(image_size);
  assert(aaafDest);
  // handle the edge cases first:

  if (centers.size() == 0)
    return;

  vector<Scalar> diameters;
  vector<Scalar> shell_thicknesses;
  vector<Scalar> voxel_intensities_foreground;
  
  if (pDiameters == nullptr) {
    // if empty, then fill the vector with the default value
    diameters.resize(centers.size(), 0.0);
    pDiameters = &diameters;
  }
  if (pShellThicknesses == nullptr) {
    // if empty, then fill the vector with the default value
    shell_thicknesses.resize(pDiameters->size());
    for (size_t i = 0; i < pDiameters->size(); i++)
      // By default, make each sphere a solid sphere by
      // setting the shell's thickness equal to the sphere's radius.
      shell_thicknesses[i] = (*pDiameters)[i] / 2;
    pShellThicknesses = &shell_thicknesses;
  }
  if (pVoxelIntensitiesForeground == nullptr) {
    // if empty, then fill the vector with the default foreground voxel value
    voxel_intensities_foreground.resize(centers.size(), 1.0);
    pVoxelIntensitiesForeground = &voxel_intensities_foreground;
  }
  assert(centers.size() == pDiameters->size());
  assert(centers.size() == pShellThicknesses->size());
  assert(centers.size() == pVoxelIntensitiesForeground->size());

  // rescale the brightnesses in the tomogram according to their average
  // value and standard deviation:

  Scalar background_ave = 0.0;
  Scalar background_stddev = 1.0;

  if (aaafBackground) {
    background_ave  =  AverageArr(image_size,
                                  aaafBackground,
                                  aaafMask);
    background_stddev  =  StdDevArr(image_size,
                                    aaafBackground,
                                    aaafMask);
  }
  double score_ave = 0.0;
  double score_rms = 0.0;
  double voxel_intensity_foreground_rms = 1.0;
  if (voxel_intensity_background_normalize) {
    for (int i = 0; i < (*pVoxelIntensitiesForeground).size(); i++)
      score_ave += (*pVoxelIntensitiesForeground)[i];
    if ((*pVoxelIntensitiesForeground).size() > 0)
      score_ave /= (*pVoxelIntensitiesForeground).size();
    for (int i = 0; i < (*pVoxelIntensitiesForeground).size(); i++)
      score_rms += SQR((*pVoxelIntensitiesForeground)[i]);
    if ((*pVoxelIntensitiesForeground).size() > 0)
      score_rms = sqrt(score_rms / (*pVoxelIntensitiesForeground).size());
    voxel_intensity_foreground_rms = score_rms;
  }
  // The spheres (or spherical shells) will have brightnesses chosen according
  // to their "scores".
  // Now determine the brightness of the background voxels.
  for (int iz = 0; iz < image_size[2]; iz++) {
    for (int iy = 0; iy < image_size[1]; iy++) {
      for (int ix = 0; ix < image_size[0]; ix++) {
        if (! voxel_intensity_background_normalize) {
          // The the background voxel brightnesses are a rescaled version of
          // the original image (and the scaling factor is 1 by default).
          aaafDest[iz][iy][ix] =
            aaafBackground[iz][iy][ix] * voxel_intensity_background_rescale;
        }
        else {
          // Then we want to automatically shift and rescale the brightess of
          // the background voxels so that the features of that image remain
          // easy to see even when (very bright?) spheres are superimposed
          // upon them.  (Why? If the spheres are many times brighter
          // or darker than the voxels in the original image, the it will
          // probably be impossible to see differenes in the brightnesses
          // of the spheres or features in the background image when viewed
          // together, because many of the voxels will appear as either
          // all black or all white when displayed on a screen.)
          if (background_stddev > 0.0)
            aaafDest[iz][iy][ix] =
              ((aaafBackground[iz][iy][ix]-background_ave) / background_stddev)
              *
              voxel_intensity_foreground_rms
              *
              voxel_intensity_background_rescale;
          else
            aaafDest[iz][iy][ix] = 0.0;
        }
      }
    }
  }

  // Now add the voxel_intensity_background brightness offset.
  for (int iz=0; iz<image_size[2]; iz++)
    for (int iy=0; iy<image_size[1]; iy++)
      for (int ix=0; ix<image_size[0]; ix++)
        aaafDest[iz][iy][ix] += voxel_intensity_background_offset;
  // (Again, all we have done so far only effects the background voxels,
  //  since we are about to overwrite the foreground voxels.)


  bool warn_points_outside_image = false;

  for (int i = 0; i < centers.size(); i++) {
    if (pReportProgress)
      *pReportProgress << "processing coordinates " << i+1 << " / "
           << centers.size() << ": x,y,z(in_voxels)="
           << centers[i][0] << "," << centers[i][1] << "," << centers[i][2];

    int ix = centers[i][0];
    int iy = centers[i][1];
    int iz = centers[i][2];

    if ((ix < 0) || (iy < 0) || (iz < 0) ||
        (ix >= image_size[0]) ||
        (iy >= image_size[1]) ||
        (iz >= image_size[2]))
      warn_points_outside_image = true;

    int Rs = ceil((*pDiameters)[i]/2-0.5);
    if (Rs < 0) Rs = 0;
    Scalar Rssqr_max = SQR((*pDiameters)[i]/2);
    Scalar Rssqr_min = 0.0;
    if (((*pShellThicknesses)[i] > 0.0) &&
        ((*pDiameters)[i]/2 - (*pShellThicknesses)[i] > 0.0))
      Rssqr_min = SQR((*pDiameters)[i]/2 - (*pShellThicknesses)[i]);

    if (pReportProgress)
      *pReportProgress << ", diameter=" << (*pDiameters)[i]
                       << ", th=" << (*pShellThicknesses)[i]
                       <<"\n";

    // Normalize the brightness of each sphere?
    // (ie by dividing the intensity by the number of voxels in the sphere)
    Scalar imultiplier = 1.0;
    long nvoxelspersphere = 1;
    if (voxel_intensity_foreground_normalize) {
      nvoxelspersphere = 0;
      for (int jz = -Rs; jz <= Rs; jz++) {
        for (int jy = -Rs; jy <= Rs; jy++) {
          for (int jx = -Rs; jx <= Rs; jx++) {
            if ((ix+jx < 0) ||
                (ix+jx >= image_size[0]) ||
                (iy+jy < 0) ||
                (iy+jy >= image_size[1]) ||
                (iz+jz < 0) ||
                (iz+jz >= image_size[2]))
              continue;
            else if (aaafMask
                     &&
                     (aaafMask[static_cast<int>(iz)+jz]
                              [static_cast<int>(iy)+jy]
                              [static_cast<int>(ix)+jx]
                      == 0.0))
              continue;
            int rsqr = jx*jx + jy*jy + jz*jz;
            if ((Rssqr_min <= rsqr) && (rsqr <= Rssqr_max))
              nvoxelspersphere++;
          }
        }
      }
    }
    imultiplier = 1.0;
    if (nvoxelspersphere > 0)
      imultiplier = 1.0 / nvoxelspersphere;
    for (int jz = -Rs; jz <= Rs; jz++) {
      for (int jy = -Rs; jy <= Rs; jy++) {
        for (int jx = -Rs; jx <= Rs; jx++) {
          int rsqr = jx*jx + jy*jy + jz*jz;
          if (! ((Rssqr_min <= rsqr) && (rsqr <= Rssqr_max)))
            continue;
          else if ((ix+jx < 0) ||
                   (ix+jx >= image_size[0]) ||
                   (iy+jy < 0) ||
                   (iy+jy >= image_size[1]) ||
                   (iz+jz < 0) ||
                   (iz+jz >= image_size[2]))
            continue;
          else if (aaafMask
                   &&
                   (aaafMask[iz+jz][iy+jy][ix+jx] == 0.0))
            continue;
          else
            aaafDest[iz+jz][iy+jy][ix+jx] =
              (*pVoxelIntensitiesForeground)[i] * imultiplier;
        }
      }
    }
  } //for (int i = 0; i < centers.size(); i++) {


  if (pReportProgress && warn_points_outside_image)
    *pReportProgress <<
      "------------------------------------------------------------------------------\n"
      "WARNING:\n"
      "    Some coordinates in the text file lie outside the boundaries of the image.\n"
      "    Did you remember to set the voxel width correctly?\n"
      "    (Did you use the \"-w\" argument?)\n"
      "--------------------------------------------------------------------------=---\n"
                     << endl;

} //DrawSpheres()



} //namespace visfd



#endif //#ifndef _DRAW_HPP
