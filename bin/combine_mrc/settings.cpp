#include <sstream>
using namespace std;
#include <err_report.h>
#include "settings.h"


string g_program_name = "merge_mrc";


vector<string> &split(const string &s, char delim, vector<string> &elems) {
  stringstream ss(s);
  string item;
  while (getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}


vector<string> split(const string &s, char delim) {
  vector<string> elems;
  split(s, delim, elems);
  return elems;
}


Settings::Settings() {
  // Default settings
  in1_file_name = "";
  in2_file_name = "";
  in_rescale01 = false;
  out_file_name = "";
  in1_use_thresholds = false;
  in1_threshold_01_a = 1.0;
  in1_threshold_01_b = 1.0;
  in1_threshold_10_a = 1.0;
  in1_threshold_10_b = 1.0;
  in2_use_thresholds = false;
  in2_threshold_01_a = 1.0;
  in2_threshold_01_b = 1.0;
  in2_threshold_10_a = 1.0;
  in2_threshold_10_b = 1.0;
  char operator_char = '+';
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


    if (vArgs[i] == "-rescale")
    {
      in_rescale01 = true;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-norescale")

    if (vArgs[i] == "-norescale")
    {
      in_rescale01 = false;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-norescale")

    // As of 2015-5-04, I don't consider "masks".
    // But I still parse the optional arguments which set them.
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

    else if (vArgs[i] == "-mask-out")
    {
      try {
        if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
          throw invalid_argument("");
        use_mask_out = true;
        mask_out = stof(vArgs[i+1]);
      }
      catch (invalid_argument& exc) {
          throw InputErr("Error: The " + vArgs[i] + 
                       " argument must be followed by a number.\n");
      }
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-mask-out")

    //Delete all the arguments we have processed in this iteration (if any)
    if (num_arguments_deleted > 0)
    {
      vArgs.erase(vArgs.begin()+i, vArgs.begin()+i+num_arguments_deleted);
      --i; //rewind by one so that when we invoke ++i
      //at the end of this loop, i's value will not change.
      //This will point us to the next unread argument.
    }
    
  } // loop over arguments "for (int i=1; i < vArgs.size(); ++i)"
  if (vArgs.size() != 5)
    throw InputErr(
      string("Error: Expected at least 4 arguments: file1.mrc operator file2.mrc out_file.mrc\n") +
      "  Usage examples:\n" +
      g_program_name + " file1.mrc  +  file2.mrc out_file.mrc\n" +
      g_program_name + " file1.mrc \"*\" file2.mrc out_file.mrc\n" +
      "  You do not have to threshold the images beforehand.  You can do it here.\n" +
      "  You can apply a discontinuous threshold (or 2-4 continuous thresholds):\n" +
      g_program_name + " file1.mrc,0.5 + file2.mrc,0.5 out_file.mrc\n" +
      g_program_name + " file1.mrc,0.48,0.49 + file2.mrc,0.48,0.49 out_file.mrc\n" +
      g_program_name + " file1.mrc,0.48,0.49,0.51,0.7 + file2.mrc,0.53,0.0.532,0.538,0.542 out_file.mrc\n" +
      "\n" +
      " Please supply two input files, an operator (+ or *) and an output file name.\n");


  vector<string> in1_file_args = split(vArgs[1], ',');

  operator_char = vArgs[2][0];

  vector<string> in2_file_args = split(vArgs[3], ',');

  out_file_name = vArgs[4];


  try {
    in1_file_name = in1_file_args[0];
    if (in1_file_args.size() > 1) {
      in1_use_thresholds = true;
      in1_threshold_01_a = stof(in1_file_args[1]);
      in1_threshold_01_b = in1_threshold_01_a;
      in1_threshold_10_a = in1_threshold_01_b;
      in1_threshold_10_b = in1_threshold_10_a;
    }
    if (in1_file_args.size() > 2) {
      in1_threshold_01_b = stof(in1_file_args[2]);
      in1_threshold_10_a = in1_threshold_01_b;
      in1_threshold_10_b = in1_threshold_10_a;
    }
    if (in1_file_args.size() > 3) {
      in1_threshold_10_a = stof(in1_file_args[3]);
      in1_threshold_10_b = in1_threshold_10_a;
    }
    if (in1_file_args.size() > 4)
      in1_threshold_10_b = stof(in1_file_args[4]);


    in2_file_name = in2_file_args[0];
    if (in2_file_args.size() > 1) {
      in2_use_thresholds = true;
      in2_threshold_01_a = stof(in2_file_args[1]);
      in2_threshold_01_b = in2_threshold_01_a;
      in2_threshold_10_a = in2_threshold_01_b;
      in2_threshold_10_b = in2_threshold_10_a;
    }
    if (in2_file_args.size() > 2) {
      in2_threshold_01_b = stof(in2_file_args[2]);
      in2_threshold_10_a = in2_threshold_01_b;
      in2_threshold_10_b = in2_threshold_10_a;
    }
    if (in2_file_args.size() > 3) {
      in2_threshold_10_a = stof(in2_file_args[3]);
      in2_threshold_10_b = in2_threshold_10_a;
    }
    if (in2_file_args.size() > 4)
      in2_threshold_10_b = stof(in2_file_args[4]);
  }
  catch (invalid_argument& exc) {
    throw InputErr("Error: The comma-separated threshold arguments following \n"
                   "       each file name must be numbers.\n");
  }

} // Settings::ParseArgs()
