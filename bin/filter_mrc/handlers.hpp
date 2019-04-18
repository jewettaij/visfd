#ifndef _HANDLERS_HPP
#define _HANDLERS_HPP

#include "settings.hpp"


void
HandleGGauss(Settings settings,
             MrcSimple &tomo_in,
             MrcSimple &tomo_out,
             MrcSimple &mask,
             float voxel_width[3]);




void
HandleGauss(Settings settings,
            MrcSimple &tomo_in,
            MrcSimple &tomo_out,
            MrcSimple &mask,
            float voxel_width[3]);




void
HandleDogg(Settings settings,
           MrcSimple &tomo_in,
           MrcSimple &tomo_out,
           MrcSimple &mask,
           float voxel_width[3]);




void
HandleDog(Settings settings,
          MrcSimple &tomo_in,
          MrcSimple &tomo_out,
          MrcSimple &mask,
          float voxel_width[3]);




void
HandleDogScaleFree(Settings settings,
                   MrcSimple &tomo_in,
                   MrcSimple &tomo_out,
                   MrcSimple &mask,
                   float voxel_width[3]);



void
HandleBlobsNonmaxSuppression(Settings settings,
                             MrcSimple &mask,
                             float voxel_width[3],
                             vector<array<float,3> >& crds,
                             vector<float>& diameters,
                             vector<float>& scores);



void
HandleDrawSpheres(Settings settings,
                  MrcSimple &tomo_in,
                  MrcSimple &tomo_out,
                  MrcSimple &mask,
                  float voxel_width[3]);



void
HandleBlobDetector(Settings settings,
                   MrcSimple &tomo_in,
                   MrcSimple &tomo_out,
                   MrcSimple &mask,
                   float voxel_width[3]);



void
HandleThresholds(Settings settings,
                 MrcSimple &tomo_in,
                 MrcSimple &tomo_out,
                 MrcSimple &mask,
                 float voxel_width[3]);



void
HandleExtrema(Settings settings,
              MrcSimple &tomo_in,
              MrcSimple &tomo_out,
              MrcSimple &mask,
              float voxel_width[3]);



void
HandleLocalFluctuations(Settings settings,
                        MrcSimple &tomo_in,
                        MrcSimple &tomo_out,
                        MrcSimple &mask,
                        float voxel_width[3]);


void
HandleWatershed(Settings settings,
                MrcSimple &tomo_in,
                MrcSimple &tomo_out,
                MrcSimple &mask,
                float voxel_width[3]);


void
HandleClusterConnected(Settings settings,
                       MrcSimple &tomo_in,
                       MrcSimple &tomo_out,
                       MrcSimple &mask,
                       float voxel_width[3]);


void
HandleRidgeDetector(Settings settings,
                    MrcSimple &tomo_in,
                    MrcSimple &tomo_out,
                    MrcSimple &mask,
                    float voxel_width[3]);



#endif //#ifndef _HANDLERS_H
