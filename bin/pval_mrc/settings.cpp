#include <cmath>
#include <sstream>
#include <iostream>
#include <cassert>
using namespace std;

#ifndef DISABLE_OPENMP
#include <omp.h>       // (OpenMP-specific)
#endif

#include <err_report.h>
#include "settings.h"

//  Example usage:
// pval_mrc -w 19.2 -i image.rec -coords crds.txt -bin-width 300.0 
// pval_mrc -w 19.2 -i image.rec -mask mask.rec -gauss 545.0 -pg -n 7145 -v 4.89109e+07
// pval_mrc -w 19.2 -i image.rec -mask mask.rec -scan 300.0 700.0 1.02


Settings::Settings() {
  // Default settings
  in_file_name = "";
  image_size[0] = 0;
  image_size[1] = 0;
  image_size[2] = 0;
  in_coords_file_name = "";
  rescale01_in = false;
  mask_file_name = "";
  mask_select = 1;
  use_mask_select = false;
  use_min_density = true;

  voxel_width = 0.0;  //How many Angstroms per voxel? (if 0 then read from file)
  voxel_width_divide_by_10 = false; //Use nm instead of Angstroms?

  filter_truncate_threshold=0.02;    //Filter intensity decay value before giving up
                            //When the filter strength is less than this value
                            //we ignore it. For difference-of-gaussian filters
                            //we choose the gaussian with the wider width. This
                            //parameter overrides other window-width settings.
                            //(Setting it to a number < 0 disables it.)

  filter_truncate_ratio = -1.0; //When averaging/filtering consider nearby
                            //voxels up to a distance of this many sigma away
                            //(Setting this to a number < 0 disables it.)

  compartment_volume = -1.0; //impossible value
  precomputed_gaussian_blur=false; //user already created a blurred density-map?
  num_particles = -1.0;      //impossible value
  use_min_density = true;
  randomize_input_image = false;
  random_seed = 0;
}



void Settings::ParseArgs(int argc, char **argv) {
  vector<string> vArgs;
  ConvertArgvToVectorOfStrings(argc, argv, vArgs);
  ParseArgs(vArgs);
}



void Settings::ConvertArgvToVectorOfStrings(int argc,
                                            char **argv,
                                            vector<string>& dest)
{
  dest.resize(argc);
  for (int i=0; i < argc; ++i)
    dest[i].assign(argv[i]);
}




void
Settings::ParseArgs(vector<string>& vArgs)
{
  for (int i=1; i < vArgs.size(); ++i)
  {

    int num_arguments_deleted = 0;

    if ((vArgs[i] == "-in") || (vArgs[i] == "-i"))
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a file name.\n");
      in_file_name = vArgs[i+1];

      num_arguments_deleted = 2;

    } // if ((vArgs[i] == "-in") || (vArgs[i] == "-i"))

    if (vArgs[i] == "-image-size")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-'))
          throw invalid_argument("");
        image_size[0] = stoi(vArgs[i+1]);
        image_size[1] = stoi(vArgs[i+2]);
        image_size[2] = stoi(vArgs[i+3]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 3 positive integers.\n");
      }

      num_arguments_deleted = 2;

    } // if ((vArgs[i] == "-in") || (vArgs[i] == "-i"))


    else if ((vArgs[i] == "-out") || (vArgs[i] == "-o"))
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a file name.\n");
      out_file_name = vArgs[i+1];

      num_arguments_deleted = 2;

    } // if ((vArgs[i] == "-out") || (vArgs[i] == "-o"))


    else if ((vArgs[i] == "-coords") || (vArgs[i] == "-crds")) {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        in_coords_file_name = vArgs[i+1];
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a file name\n");
      }
      num_arguments_deleted = 2;
    }


    else if (vArgs[i] == "-w") {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        voxel_width = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by voxel width.\n");
      }
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-w")


    else if ((vArgs[i] == "-a2nm") || (vArgs[i] == "-ang-to-nm"))
    {
      voxel_width_divide_by_10 = true;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-a2nm")

    
    else if ((vArgs[i] == "-randomize") ||
             (vArgs[i] == "-rand") ||
             (vArgs[i] == "-rand"))
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a positive integer.\n");
      randomize_input_image = true;
      random_seed = stoi(vArgs[i+1]);
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-a2nm")


    else if (vArgs[i] == "-mask")
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a file name.\n");
      mask_file_name = vArgs[i+1];
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-mask")


    else if (vArgs[i] == "-mask-select")
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        use_mask_select = true;
        mask_select = stoi(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by an integer.\n");
      }
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-mask-select")


    else if (vArgs[i] == "-gauss")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        vfSigma.resize(1);
        vfSigma[0] = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a positive number (\"s\"),\n"
                       " the Gaussian width\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-gauss")


    else if (vArgs[i] == "-bin-width")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        vfSigma.resize(1);
        // See the comment in "extreme_prob.cpp" regarding the relationship
        // between "bin volume" and Gaussian width (sigma).  Short version:
        // bin_vol = (2*pi*sigma^2)^(3/2)
        // bin_width = bin_vol ^ 1/3
        //  <==> sigma = bin_width / sqrt(2*pi);
        vfSigma[0] = stof(vArgs[i+1]) / sqrt(2*M_PI);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a positive number (\"s\"),\n"
                       " the bin width\n");
      }
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-bin-width")


    else if ((vArgs[i] == "-scan") ||
             (vArgs[i] == "-scan-bin") ||
             (vArgs[i] == "-scan-bin-width") ||
             (vArgs[i] == "-scan-gauss"))
    {
      try {
        if ((i+3 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
            (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
            (vArgs[i+3] == "") || (vArgs[i+3][0] == '-'))
          throw invalid_argument("");

        float bin_width_min = stof(vArgs[i+1]);
        float bin_width_max = stof(vArgs[i+2]);
        if ((bin_width_min <= 0.0) ||
            (bin_width_max <= 0.0) ||
            (bin_width_min >= bin_width_max))
          throw invalid_argument("");

        float growth_ratio = stof(vArgs[i+3]);
        if (growth_ratio <= 1.0)
          throw invalid_argument("");
        int N;
        if (bin_width_min == bin_width_max)
          N = 1;
        else {
          N = 1 + ceil( log(bin_width_max/bin_width_min) / log(growth_ratio) );
          growth_ratio = pow(bin_width_max/bin_width_min, 1.0/N);
        }

        vfSigma.resize(N);

        // See the comment in "extreme_prob.cpp" regarding the relationship
        // between "bin volume" and Gaussian width (sigma).  Short version:
        // bin_vol = (2*pi*sigma^2)^(3/2)
        // bin_width = bin_vol ^ 1/3
        //  <==> sigma = bin_width / sqrt(2*pi);
        float first_bin_width = bin_width_min;
        float first_sigma = first_bin_width;
        if (vArgs[i] != "-scan-gauss")
          first_sigma /= sqrt(2*M_PI);
        vfSigma[0] = first_sigma;
        for (int n = 1; n < N; n++) {
          vfSigma[n] = vfSigma[n-1] * growth_ratio;
        }
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 3 positive numbers:\n"
                       "      bin_width_min  bin_width_max  growth_ratio N\n"
                       " The first two numbers should be in increasing order\n"
                       " and the last number should be > 1.0\n");
      }
      num_arguments_deleted = 4;
    } //if (vArgs[i] == "-scan")

    
    else if (vArgs[i] == "-window-ratio")
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        filter_truncate_ratio = stof(vArgs[i+1]);
        filter_truncate_threshold=-1.0; //(disables)override any filter_truncate_threshold settings
        //filter_truncate_ratio_exp = exp(-pow(filter_truncate_ratio,
        //                                   n_exp));
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number.\n");
      }
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-window-ratio")


    else if ((vArgs[i] == "-cutoff") || (vArgs[i] == "-window-cutoff"))
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        filter_truncate_threshold = stof(vArgs[i+1]);
        filter_truncate_ratio = -1.0; //(disables) override any filter_truncate_ratio settings
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number.\n");
      }      
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-cutoff")


    else if ((vArgs[i] == "-vol") ||
             (vArgs[i] == "-volume") ||
             (vArgs[i] == "-v"))
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") ||
            (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        compartment_volume = stof(vArgs[i+1]);
        if (compartment_volume < 0)
          throw invalid_argument("");
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a nonnegativenumber.\n");
      }      
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-vol")


    else if (vArgs[i] == "-n")
    {
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") ||
            (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        num_particles = stof(vArgs[i+1]);
        if (num_particles < 0)
          throw invalid_argument("");
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a nonnegativenumber.\n");
      }      
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-n")


    else if ((vArgs[i] == "-pre-gauss") ||
             (vArgs[i] == "-pg"))
    {
      precomputed_gaussian_blur = true;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-pre-gauss")


    else if ((vArgs[i] == "-min") ||
             (vArgs[i] == "-minima"))
    {
      use_min_density = true;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-min")


    else if ((vArgs[i] == "-max") ||
             (vArgs[i] == "-maxima"))
    {
      use_min_density = false;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-min")



    else if (vArgs[i] == "-np") {
      #ifdef DISABLE_OPENMP
      throw InputErr("Error: The " + vArgs[i] + 
                     " argument is only available if program was compiled\n"
                     " with support for OpenMP (multiprocessor support).\n");
      #else
      try {
        if ((i+1 >= vArgs.size()) ||
            (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        int num_threads = stoi(vArgs[i+1]);
        omp_set_num_threads(num_threads);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a positive integer, the\n"
                       "       number of threads (processors) requested.\n");
      }
      #endif //#ifdef DISABLE_OPENMP
      num_arguments_deleted = 2;
    }

    //Delete all the arguments we have processed in this iteration (if any)
    if (num_arguments_deleted > 0)
    {
      vArgs.erase(vArgs.begin()+i, vArgs.begin()+i+num_arguments_deleted);
      --i; //rewind by one so that when we invoke ++i
      //at the end of this loop, i's value will not change.
      //This will point us to the next unread argument.
    }

  } // loop over arguments "for (int i=1; i < vArgs.size(); ++i)"


  if (vArgs.size() > 1) {
    stringstream err_msg;
    err_msg << "Error: Unrecognized argument: \""
            << vArgs[1] << "\"\n";
    throw InputErr(err_msg.str());
  }


  if ((in_file_name == "") &&
      ((image_size[0] == 0) || (image_size[1]==0) || (image_size[2]==0)))
    throw InputErr("Error: You need to specify the size of the image that will be used to\n"
                   "       created to contain the density cloud.  You can EITHER specify this\n"
                   "       using the \"-in filename\" OR \"-image-size Nx Ny Nz\" arguments.\n");



  if (vfSigma.size() == 0) {
    throw InputErr("Error: You need to specify the effective bin width using \"-gauss sigma\".\n"
                   "       ...where \"sigma\" is the width of the Gaussian used for bluring.\n"
                   "       OR using \"-scan min_bin max_bin growth_ratio\"\n");
  }


} //Settings::ParseArgs(vector<string>& vArgs)
