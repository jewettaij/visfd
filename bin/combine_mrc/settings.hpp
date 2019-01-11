#ifndef _SETTINGS_HPP
#define _SETTINGS_HPP



#include<vector>
#include<string>
using namespace std;



class Settings {
 public:
  string in1_file_name;
  string in2_file_name;
  bool in_rescale01;
  string out_file_name;
  char operator_char = '+';
  bool  in1_use_thresholds;
  float in1_threshold_01_a;
  float in1_threshold_01_b;
  float in1_threshold_10_a;
  float in1_threshold_10_b;
  bool  in2_use_thresholds;
  float in2_threshold_01_a;
  float in2_threshold_01_b;
  float in2_threshold_10_a;
  float in2_threshold_10_b;

  // As of 2015-5-04, I don't consider "masks".  Please ignore next 5 arguments:
  string mask_file_name;
  int mask_select;
  bool use_mask_select;
  int mask_out;
  bool use_mask_out;

  Settings();
  void ParseArgs(int argc, char **argv);
 private:
  void ParseArgs(vector<string>& vArgs);
  void ConvertArgvToVectorOfStrings(int argc, 
                                    char **argv, 
                                    vector<string>& dest);
};


#endif //#ifndef _SETTINGS_HPP
