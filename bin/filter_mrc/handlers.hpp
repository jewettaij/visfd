#ifndef _HANDLERS_HPP
#define _HANDLERS_HPP

#include "settings.hpp"

void
HandleDilation(const Settings &settings,
               MrcSimple &tomo_in,
               MrcSimple &tomo_out,
               MrcSimple &mask,
               float voxel_width[3]);


void
HandleErosion(const Settings &settings,
              MrcSimple &tomo_in,
              MrcSimple &tomo_out,
              MrcSimple &mask,
              float voxel_width[3]);


void
HandleOpening(const Settings &settings,
              MrcSimple &tomo_in,
              MrcSimple &tomo_out,
              MrcSimple &mask,
              float voxel_width[3]);


void
HandleClosing(const Settings &settings,
              MrcSimple &tomo_in,
              MrcSimple &tomo_out,
              MrcSimple &mask,
              float voxel_width[3]);


void
HandleGGauss(const Settings &settings,
             MrcSimple &tomo_in,
             MrcSimple &tomo_out,
             MrcSimple &mask,
             float voxel_width[3]);


void
HandleGauss(const Settings &settings,
            MrcSimple &tomo_in,
            MrcSimple &tomo_out,
            MrcSimple &mask,
            float voxel_width[3]);


void
HandleDogg(const Settings &settings,
           MrcSimple &tomo_in,
           MrcSimple &tomo_out,
           MrcSimple &mask,
           float voxel_width[3]);


void
HandleDog(const Settings &settings,
          MrcSimple &tomo_in,
          MrcSimple &tomo_out,
          MrcSimple &mask,
          float voxel_width[3]);


void
HandleLoGDoG(const Settings &settings,
             MrcSimple &tomo_in,
             MrcSimple &tomo_out,
             MrcSimple &mask,
             float voxel_width[3]);


void
HandleBlobsNonmaxSuppression(const Settings &settings,
                             MrcSimple &mask,
                             float voxel_width[3],
                             vector<array<float,3> >& crds,
                             vector<float>& diameters,
                             vector<float>& scores);


void
HandleBlobScoreSupervisedMulti(const Settings &settings,
                               float voxel_width[3]);



void
HandleDrawSpheres(const Settings &settings,
                  MrcSimple &tomo_in,
                  MrcSimple &tomo_out,
                  MrcSimple &mask,
                  float voxel_width[3]);


void
HandleBlobDetector(const Settings &settings,
                   MrcSimple &tomo_in,
                   MrcSimple &tomo_out,
                   MrcSimple &mask,
                   float voxel_width[3]);


void
HandleThresholds(const Settings &settings,
                 MrcSimple &tomo_in,
                 MrcSimple &tomo_out,
                 MrcSimple &mask,
                 float voxel_width[3]);



void
HandleExtrema(const Settings &settings,
              MrcSimple &tomo_in,
              MrcSimple &tomo_out,
              MrcSimple &mask,
              float voxel_width[3]);


void
HandleLocalFluctuations(const Settings &settings,
                        MrcSimple &tomo_in,
                        MrcSimple &tomo_out,
                        MrcSimple &mask,
                        float voxel_width[3]);


void
HandleWatershed(const Settings &settings,
                MrcSimple &tomo_in,
                MrcSimple &tomo_out,
                MrcSimple &mask,
                float voxel_width[3]);


void
HandleClusterConnected(const Settings &settings,
                       MrcSimple &tomo_in,
                       MrcSimple &tomo_out,
                       MrcSimple &mask,
                       float voxel_width[3]);


void
HandleTV(const Settings &settings,
         MrcSimple &tomo_in,
         MrcSimple &tomo_out,
         MrcSimple &mask,
         float voxel_width[3]);


#endif //#ifndef _HANDLERS_H
