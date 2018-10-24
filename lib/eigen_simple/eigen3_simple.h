#ifndef _EIGEN3_SIMPLE_H
#define _EIGEN3_SIMPLE_H

// This file contains code from "Eigen", a lightweight C++ template library
// for linear algebra.  (It was later simplified to remove some dependencies.)
//
// Copyright (C) 2008-2010 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2010 Jitse Niesen <jitse@maths.leeds.ac.uk>
// Copyright (C) 2018 Andrew Jewett <jewett@scripps.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <cmath>
#include <limits>
#include <cstring>
#include <cassert>
using namespace std;


#include "lin3_utils.h"



namespace selfadjoint_eigen3
{

  template <class Scalar>
  static inline void computeRoots3(const Scalar m[3][3],
                                   Scalar roots[3])
  {
    const Scalar s_inv3 = 1.0/3.0;
    const Scalar s_sqrt3 = std::sqrt(3.0);
 
    // The characteristic equation is x^3 - c2*x^2 + c1*x - c0 = 0.  The
    // eigenvalues are the roots to this equation, all guaranteed to be
    // real-valued, because the matrix is symmetric.
    Scalar c0 = m[0][0]*m[1][1]*m[2][2] + 2.0*m[1][0]*m[2][0]*m[2][1] - m[0][0]*m[2][1]*m[2][1] - m[1][1]*m[2][0]*m[2][0] - m[2][2]*m[1][0]*m[1][0];
    Scalar c1 = m[0][0]*m[1][1] - m[1][0]*m[1][0] + m[0][0]*m[2][2] - m[2][0]*m[2][0] + m[1][1]*m[2][2] - m[2][1]*m[2][1];
    Scalar c2 = m[0][0] + m[1][1] + m[2][2];
 
    // Construct the parameters used in classifying the roots of the equation
    // and in solving the equation for the roots in closed form.
    Scalar c2_over_3 = c2*s_inv3;
    Scalar a_over_3 = (c2*c2_over_3 - c1)*s_inv3;
    a_over_3 = std::max(a_over_3, 0.0);
 
    Scalar half_b = 0.5*(c0 + c2_over_3*(2.0*c2_over_3*c2_over_3 - c1));
 
    Scalar q = a_over_3*a_over_3*a_over_3 - half_b*half_b;
    q = std::max(q, 0.0);
 
    // Compute the eigenvalues by solving for the roots of the polynomial.
    Scalar rho = std::sqrt(a_over_3);
    Scalar theta = std::atan2(std::sqrt(q),half_b)*s_inv3;  // since sqrt(q) > 0, atan2 is in [0, pi] and theta is in [0, pi/3]
    Scalar cos_theta = std::cos(theta);
    Scalar sin_theta = std::sin(theta);
    // roots are already sorted, since cos is monotonically decreasing on [0, pi]
    roots[0] = c2_over_3 - rho*(cos_theta + s_sqrt3*sin_theta); // == 2*rho*cos(theta+2pi/3)
    roots[1] = c2_over_3 - rho*(cos_theta - s_sqrt3*sin_theta); // == 2*rho*cos(theta+ pi/3)
    roots[2] = c2_over_3 + 2.0*rho*cos_theta;
  }
 


  template <class Scalar>
  static inline void extract_kernel3(const Scalar mat[3][3],
                                     Scalar res[3],
                                     Scalar representative[3])
  {
    using std::abs;

    // Find non-zero column i0 (by construction, there must exist a
    // non zero coefficient on the diagonal):
    int i0 = 0;
    Scalar max_diag = abs(mat[0][0]);
    for (int d=1; d<3; d++) {
      if (abs(mat[d][d]) > max_diag) {
        i0 = d;
        max_diag = abs(mat[d][d]);
      }
    }

    Scalar matT[3][3];
    Transpose3(mat, matT);

    // matT[i0] = i0'th column of mat  (ie, mat[0][i0], mat[1][i0], mat[2][i0])
    // This is a good candidate for an orthogonal vector to 
    // the current eigenvector, so let's save it:
    for (int d=0; d<3; d++) {
      representative[d] = matT[i0][d]; // = mat[d][i0];
    }

    Scalar c0[3], c1[3];

    CrossProduct(representative, matT[ (i0+1)%3 ], c0);
    CrossProduct(representative, matT[ (i0+2)%3 ], c1);
    Scalar n0 = SquaredNorm3(c0);
    Scalar n1 = SquaredNorm3(c1);

    if(n0 > n1) {
      Scalar norm_c = 1.0 / std::sqrt(n0);
      for (int d=0; d<3; d++)
        res[d] = c0[d] * norm_c;
    }
    else {
      Scalar norm_c = 1.0 / std::sqrt(n1);
      for (int d=0; d<3; d++)
        res[d] = c1[d] * norm_c;
    }
 
  } //extract_kernel3()
 


  template <class Scalar>
  inline void Diagonalize3(const Scalar mat[3][3],
                           Scalar eivals[3], // store eigenvalues here
                           // If the user supplies a 2D eivects array,
                           // then the 3 eigenvectors will be stored in
                           // the 3 rows of the "eivects" 2D array
                           // (ie, eivects[0], eivects[1], eivects[2])
                           // Note: eivects
                           Scalar eivects[][3] = NULL)
  {
    const Scalar EPSILON = std::numeric_limits<Scalar>::epsilon();
    // Shift the matrix to the mean eigenvalue and map the matrix
    // coefficients to [-1:1] to avoid over- and underflow.
    Scalar shift = trace3(mat) / 3.0;
    Scalar scaledMat[3][3];
    Copy3(mat, scaledMat);
    // Apply the shift
    for (int d=0; d<3; d++)
      scaledMat[d][d] -= shift;
    // Scale the matrix so that all entries lie between -1 and 1
    Scalar scale = -1.0;
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        if (abs(scaledMat[i][j]) > scale)
          scale = abs(scaledMat[i][j]);
    assert(scale >= 0.0);
    if (scale > 0) {
      Scalar scale_inv = 1.0/scale;
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
          scaledMat[i][j] *= scale_inv;
    }
 
    // compute the eigenvalues
    computeRoots3(scaledMat, eivals);
 
    // compute the eigenvectors ?
    if (eivects)
    {
      // If non-null, eivects is assumed to be a 3x3 matrix
      // The i'th row of this matrix will contain the the eigenvector
      // corresponding to the i'th eigenvalue.
      if((eivals[2]-eivals[0])<=EPSILON)
      {
        for (int i=0; i<3; i++)
          for (int j=0; j<3; j++)
            // If all three eigenvalues are numerically the same, then
            // by default, set the eigenvectors to the x, y, z axes.
            eivects[i][j] = ((i == j) ? 1.0 : 0.0);
      }
      else
      {
        // Compute the eigenvector of the most distinct eigenvalue
        Scalar d0 = eivals[2] - eivals[1];
        Scalar d1 = eivals[1] - eivals[0];
        int k=0;
        int l=2;
        if(d0 > d1)
        {
          d0 = d1;
          // swap k and l
          int _k = l;
          l = k;
          k = _k;
        }
 
        Scalar tmp[3][3];
        Copy3(scaledMat, tmp);
        // Compute the eigenvector of index k
        {
          // subtract eivals[k] * Identity from "tmp"
          for (int d=0; d<3; d++)
            tmp[d][d] -= eivals[k];
          // By construction, 'tmp' is of rank 2, and its kernel
          // corresponds to the respective eigenvector.
          extract_kernel3(tmp, eivects[k], eivects[l]);
        }
 
        // Compute eigenvector of index l
        if(d0<=2*EPSILON*d1)
        {
          // If d0 is too small, then the two other eigenvalues
          // are numerically the same,
          // and thus we only have to ortho-normalize the near
          // orthogonal vector we saved above.
          Scalar k_dot_l = DotProduct3(eivects[k], eivects[l]);
          for (int d=0; d<3; d++)
            eivects[l][d] -= k_dot_l * eivects[l][d];
          Normalize3(eivects[l]);
        }
        else
        {
          Copy3(scaledMat, tmp);
          for (int d=0; d<3; d++)
            tmp[d][d] -= eivals[l];

          Scalar dummy[3];
          extract_kernel3(tmp, eivects[l], dummy);
        }

        // Compute last eigenvector (eivects[1]) from the other two
        CrossProduct(eivects[2], eivects[0], eivects[1]);
        Normalize3(eivects[1]);
      }
    }
 
    // Rescale back to the original size.
    for (int d=0; d<3; d++) {
      eivals[d] *= scale;
      eivals[d] += shift;
    }

  } //Diagonalize3()


} // namespace selfadjoint_eigenvalues3




#endif //#ifndef _EIGEN3_SIMPLE_H
