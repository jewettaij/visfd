#ifndef _UNSUPPORTED_HPP
#define _UNSUPPORTED_HPP

/// UNSUPPORTED CODE
/// DEPRECIATION WARNING:
/// THIS FILE IMPLEMENTS FEATURES WHICH ARE NO LONGER MAINTAINED AND TESTED.
/// THIS CODE WILL LIKELY BE DELETED IN THE FUTURE.


#include "settings.hpp"



/// @brief
///  Perform template matching with a Gaussian.
///  DEPRECIATION WARNING: This function may be removed in the future.
///  A comment describing of this function is provided in "handlers_unsupported.cpp"

void
HandleTemplateGauss(const Settings &settings,
                    MrcSimple &tomo_in,
                    MrcSimple &tomo_out,
                    MrcSimple &mask,
                    float voxel_width[3]);



#ifndef DISABLE_TEMPLATE_MATCHING
/// @brief
///  Perform template matching with a (spherically symmetric)
///  generalized Gaussian.
///  DEPRECIATION WARNING: This function may be removed in the future.
///  A comment describing of this function is provided in "handlers_unsupported.cpp"

void
HandleTemplateGGauss(const Settings &settings,
                     MrcSimple &tomo_in,
                     MrcSimple &tomo_out,
                     MrcSimple &mask,
                     float voxel_width[3]);

#endif //#ifndef DISABLE_TEMPLATE_MATCHING



#ifndef DISABLE_BOOTSTRAPPING
/// @brief
/// DEPRECIATION WARNING: This function will probably be removed in the future.
void
HandleBootstrapDogg(const Settings &settings,
                    MrcSimple &tomo_in,
                    MrcSimple &tomo_out,
                    MrcSimple &mask);

#endif //#ifndef DISABLE_BOOTSTRAPPING



#ifndef DISABLE_DOGGXY
void
HandleDoggXY(const Settings &settings,
             MrcSimple &tomo_in,
             MrcSimple &tomo_out,
             MrcSimple &mask,
             float voxel_width[3]);
#endif //#ifndef DISABLE_DOGGXY



#ifndef DISABLE_INTENSITY_PROFILES
void
HandleBlobIntensityProfiles(const Settings &settings,
                            MrcSimple &tomo_in,
                            MrcSimple &tomo_out,
                            MrcSimple &mask);
#endif //#ifndef DISABLE_INTENSITY_PROFILES



#ifndef DISABLE_INTENSITY_PROFILES
void
HandleBlobRadialIntensity(const Settings &settings,
                          MrcSimple &tomo_in,
                          MrcSimple &tomo_out,
                          MrcSimple &mask,
                          float voxel_width[3]);
#endif


void
HandleDistanceToPoints(const Settings &settings,
                       MrcSimple &tomo_in,
                       MrcSimple &tomo_out,
                       MrcSimple &mask,
                       float voxel_width[3]);


#endif //#ifndef _UNSUPPORTED_HPP
