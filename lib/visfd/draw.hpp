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
#include <alloc2d.hpp>    // defines Alloc2D() and Dealloc2D()
#include <alloc3d.hpp>    // defines Alloc3D() and Dealloc3D()
#include <filter1d.hpp>   // defines "Filter1D" (used in ApplySeparable())
#include <filter3d.hpp>   // defines common 3D image filters







/// @brief  Create a 3D image containing multiple spheres (or spherical shells)
///         various sizes, thicknesses, and locations, specified by the caller.
///         The resulting spheres can (optionally) be superimposed with the
///         existing image (if the aaafBackground image array is != NULL).

template<class Scalar>
void
VisualizeBlobs(int const image_size[3], //!< image size
               Scalar ***aaafDest,  //!< array where we should write new image
               Scalar const *const *const *aaafMask,   //!< Optional: ignore voxels where mask==0
               const vector<array<Scalar,3> > &centers, //!< coordinates for the center of each sphere (blob)
               const vector<Scalar> *pDiameters=NULL,         //!< Optional: diameter of each spherical shell (in voxels)
               const vector<Scalar> *pShellThicknesses=NULL, //!< Optional: thickness of each spherical shell (in voxels)
               const vector<Scalar> *pVoxelIntensitiesForeground=NULL, //!< Optional: assign voxels in spherical shell to this value (the vector should contain a separate entry for every sphere)
               Scalar voxel_intensity_background = 0.0, //!< Optional: assign background voxels to this value
               Scalar const *const *const *aaafBackground = NULL,   //!< Optional: superimpose background image?
               Scalar voxel_intensity_background_rescale = 0.333, //!< Optional: superimpose with old image? This is the ratio of the fluctuations in voxel intensities of the newly created background image relative to the average foreground voxel intensity.
               bool voxel_intensity_foreground_normalize = false, //!< Optional: divide brightnesses by number of voxels in spherical shell? (rarely useful)
               ostream *pReportProgress = NULL //!<Optional: report progress to the user?
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
  
  if (pDiameters == NULL) {
    // if empty, then fill the vector with the default value
    diameters.resize(centers.size(), 0.0);
    pDiameters = &diameters;
  }
  if (pShellThicknesses == NULL) {
    // if empty, then fill the vector with the default value
    shell_thicknesses.resize(centers.size(), 1.0);
    pShellThicknesses = &shell_thicknesses;
  }
  if (pVoxelIntensitiesForeground == NULL) {
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
  for (int i = 0; i < (*pVoxelIntensitiesForeground).size(); i++)
    score_ave += (*pVoxelIntensitiesForeground)[i];
  if ((*pVoxelIntensitiesForeground).size() > 0)
    score_ave /= (*pVoxelIntensitiesForeground).size();
  double voxel_intensity_foreground_ave = score_ave;

  // The spherical shells will have brightnesses chosen according to their
  // scores. Rescale the tomogram intensities so that they are approximately
  // in the range of scores.  This way we can see both at the same time
  // (when viewing the tomogram using IMOD, for example)
  for (int iz = 0; iz < image_size[2]; iz++) {
    for (int iy = 0; iy < image_size[1]; iy++) {
      for (int ix = 0; ix < image_size[0]; ix++) {
        if (background_stddev > 0.0)
          aaafDest[iz][iy][ix] =
            ((aaafBackground[iz][iy][ix] - background_ave) / background_stddev)
            *
            voxel_intensity_foreground_ave
            *
            voxel_intensity_background_rescale;
        else
          aaafDest[iz][iy][ix] = 0.0;
      }
    }
  }


  for (int iz=0; iz<image_size[2]; iz++)
    for (int iy=0; iy<image_size[1]; iy++)
      for (int ix=0; ix<image_size[0]; ix++)
        aaafDest[iz][iy][ix] += voxel_intensity_background;


  for (int i = 0; i < centers.size(); i++) {
    if (pReportProgress)
      *pReportProgress << "processing coordinates " << i+1 << " / "
           << centers.size() << ": x,y,z(in_voxels)="
           << centers[i][0] << "," << centers[i][1] << "," << centers[i][2];

    if ((centers[i][0] < 0) || (centers[i][1] < 0) || (centers[i][2] < 0) ||
        (centers[i][0] >= image_size[0]) ||
        (centers[i][1] >= image_size[1]) ||
        (centers[i][2] >= image_size[2]))
        throw VisfdErr("Error: Coordinates in the text file lie outside the boundaries of the image.\n"
                       "       Did you set the voxel width correctly?\n"
                       "       (Did you use the \"-w\" argument?)\n");
      
    if (aaafMask &&
        (aaafMask[static_cast<int>(centers[i][2])]
                 [static_cast<int>(centers[i][1])]
                 [static_cast<int>(centers[i][0])] == 0.0))
      continue;

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
            Scalar rsqr = jx*jx + jy*jy + jz*jz;
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
          else if ((centers[i][0] + jx < 0) ||
                   (centers[i][0] + jx >= image_size[0]) ||
                   (centers[i][1] + jy < 0) ||
                   (centers[i][1] + jy >= image_size[1]) ||
                   (centers[i][2] + jz < 0) ||
                   (centers[i][2] + jz >= image_size[2]))
            continue;
          else if (aaafMask
                   &&
                   (aaafMask[static_cast<int>(centers[i][2])+jz]
                            [static_cast<int>(centers[i][1])+jy]
                            [static_cast<int>(centers[i][0])+jx]
                    == 0.0))
            continue;
          else
            aaafDest[static_cast<int>(centers[i][2])+jz]
                    [static_cast<int>(centers[i][1])+jy]
                    [static_cast<int>(centers[i][0])+jx] =
                         (*pVoxelIntensitiesForeground)[i] * imultiplier;
        }
      }
    }
  } //for (int i = 0; i < centers.size(); i++) {

} //VisualizeBlobs()







#endif //#ifndef _DRAW_HPP
