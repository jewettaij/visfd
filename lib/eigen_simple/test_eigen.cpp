#include <iostream>
using namespace std;

#include "eigen3_simple.h"
using namespace selfadjoint_eigen3;

int main(int argc, char **argv) {

  //double M[3][3] = {{2.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,3.0}};
  double M[3][3] = {{-1.3,2,4},{2,-2.3,3},{4,3,-5.3}};
  double eivals[3];
  double eivects[3][3];

  Diagonalize3(M, eivals, eivects);

  cerr << "eivals: "
       << eivals[0] << " " << eivals[1] << " " << eivals[2] << endl;
  cerr << "eivects[0]: "
       << eivects[0][0] << " " << eivects[0][1] << " " << eivects[0][2] << endl;
  cerr << "eivects[1]: "
       << eivects[1][0] << " " << eivects[1][1] << " " << eivects[1][2] << endl;
  cerr << "eivects[2]: "
       << eivects[2][0] << " " << eivects[2][1] << " " << eivects[2][2] << endl;
}

