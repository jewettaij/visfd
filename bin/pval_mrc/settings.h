#ifndef _SETTINGS_H
#define _SETTINGS_H

#define DISABLE_BOOTSTRAPPING
//#define DISABLE_TEMPLATE_MATCHING


#include<vector>
#include<string>
using namespace std;

class Settings {
 public:
  string in_file_name;
  string in_coords_file_name;
  bool rescale01_in;
  string mask_file_name;
  int mask_select;
  bool use_mask_select;


  float voxel_width = 0.0;  //How many Angstroms per voxel? (if 0 then read from file)
  bool voxel_width_divide_by_10 = false; //Use nm instead of Angstroms?

  //float sigma;              //width of Gaussian blur to apply to image
  vector<float>  vfSigma;   //a list of widths to use for Gaussian blurring

  float filter_truncate_threshold;    //Filter intensity decay value before giving up
                            //When the filter strength is less than this value
                            //we ignore it. For difference-of-gaussian filters
                            //we choose the gaussian with the wider width. This
                            //parameter overrides other window-width settings.
                            //(Setting it to a number < 0 disables it.)

  float filter_truncate_ratio; //When averaging/filtering consider nearby
                            //voxels up to a distance of this many sigma away
                            //(Setting this to a number < 0 disables it.)

  bool precomputed_gaussian_blur; //did the user already use a Gaussian blur?
  // Did the user already precompute the volume and number of particles?
  float compartment_volume; //the volume of the cell (if known in advance)
  float num_particles;      //number of particles in this volume (if known)

  bool use_min_density;     //Are we interested in density minima or maxima?

  
  Settings();
  void ParseArgs(int argc, char **argv);

 private:
  void ParseArgs(vector<string>& vArgs);
  void ConvertArgvToVectorOfStrings(int argc, 
                                    char **argv, 
                                    vector<string>& dest);

}; //class Settings


#endif //#ifndef _SETTINGS_H
