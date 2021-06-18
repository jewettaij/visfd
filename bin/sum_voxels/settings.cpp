#include <cmath>
#include <sstream>
#include <iostream>
#include <cassert>
using namespace std;

#include "err.hpp"
#include "settings.hpp"


Settings::Settings() {
  // Default settings
  in_file_name = "";
  rescale01_in = false;
  mask_file_name = "";
  mask_select = 1;
  use_mask_select = false;

  calc_ave = false;
  calc_stddev = false;

  multiply_by_voxel_volume = false;
  voxel_width = 0.0;  //How many Angstroms per voxel? (if 0 then read from file)
  voxel_width_divide_by_10 = false; //Use nm instead of Angstroms?

  use_thresholds = false;
  use_dual_thresholds = false;
  in_threshold_01_a = 1.0;
  in_threshold_01_b = 1.0;
  in_threshold_10_a = 1.0;
  in_threshold_10_b = 1.0;
  in_thresh2_use_clipping = false;
  in_thresh_a_value = 0.0;
  in_thresh_b_value = 1.0;
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


    if ((vArgs[i] == "-volume") ||
        (vArgs[i] == "-vol") ||
        (vArgs[i] == "-integral") ||
        (vArgs[i] == "-integrate")) {
      multiply_by_voxel_volume = true;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-volume")



    else if ((vArgs[i] == "-ave") ||
             (vArgs[i] == "-average")) {
      calc_ave = true;
      num_arguments_deleted = 1;
    }


    else if (vArgs[i] == "-stddev") {
      calc_stddev = true;
      num_arguments_deleted = 1;
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


    else if ((vArgs[i] == "-no-rescale") || (vArgs[i] == "-norescale"))
    {
      rescale01_in = false;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-norescale")


    else if ((vArgs[i] == "-rescale") || (vArgs[i] == "-rescale-in"))
    {
      rescale01_in = true;
      in_threshold_01_a = 0.0;
      in_threshold_01_b = 1.0;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-rescale")


    else if (vArgs[i] == "-clip") {
      try {
        if (i+2 >= vArgs.size())
          throw invalid_argument("");
        use_thresholds = true;
        use_dual_thresholds = false;
        in_threshold_01_a = stof(vArgs[i+1]);
        in_threshold_01_b = stof(vArgs[i+2]);
        in_thresh2_use_clipping = true;
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 2 numbers.\n");
      }
      num_arguments_deleted = 3;
    }


    else if (vArgs[i] == "-thresh4") {
      try {
        if (i+4 >= vArgs.size())
          throw invalid_argument("");
        use_thresholds = true;
        use_dual_thresholds = true;
        in_threshold_01_a = stof(vArgs[i+1]);
        in_threshold_01_b = stof(vArgs[i+2]);
        in_threshold_10_a = stof(vArgs[i+3]);
        in_threshold_10_b = stof(vArgs[i+4]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 4 numbers.\n");
      }
      num_arguments_deleted = 5;
    }


    else if (vArgs[i] == "-thresh-range") {
      try {
        if (i+2 >= vArgs.size())
          throw invalid_argument("");
        use_dual_thresholds = true;
        in_threshold_01_a = stof(vArgs[i+1]);
        in_threshold_01_b = stof(vArgs[i+1]);
        in_threshold_10_a = stof(vArgs[i+2]);
        in_threshold_10_b = stof(vArgs[i+2]);
      }
      catch (invalid_argument& exc) {
        throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by 4 numbers.\n");
      }
      num_arguments_deleted = 3;
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

  if (vArgs.size() < 2) {
    throw InputErr("Error: Expected a file name.\n");
  }
  if (vArgs.size() > 2) {
    stringstream err_msg;
    err_msg << "Error: Too many arguments or mistake in argument syntax. Unrecognized args:\n";
    for (int i = 2; i < vArgs.size(); i++)
      err_msg << "       \"" << vArgs[i] << "\"\n";
    throw InputErr(err_msg.str());
  }

  if (calc_ave && calc_stddev)
    throw InputErr("Error: You can request \"-ave\" or \"-stddev\", but not both.\n");

  if (multiply_by_voxel_volume && (calc_ave || calc_stddev)) {
    throw InputErr("Error: The \"-ave\" or \"-stddev\" arguments do not make sense when used\n"
                   "       simulatneously with the \"-volume\" argument.\n"
                   "       (Omit the \"-volume\" argument.)\n");
  }
    
  in_file_name = vArgs[1];

} //Settings::ParseArgs(vector<string>& vArgs)
