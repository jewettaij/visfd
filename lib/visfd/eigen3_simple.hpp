#ifndef _EIGEN3_SIMPLE_HPP
#define _EIGEN3_SIMPLE_HPP

// This file contains code from "Eigen", a lightweight C++ template library
// for linear algebra.
// (Eigen is quite large.  The code in this directory is a simplified version
//  which removes the many dependencies on other files from the Eigen library.)
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
#include <type_traits>
using namespace std;


#include "lin3_utils.hpp"


namespace visfd {


namespace selfadjoint_eigen3 {


  typedef enum eEigenOrderType {
    INCREASING_EIVALS,
    DECREASING_EIVALS,
    INCREASING_ABS_EIVALS,
    DECREASING_ABS_EIVALS,
    INCREASINGLY_DISTINCT_EIVALS,
    DECREASINGLY_DISTINCT_EIVALS,
  } EigenOrderType;



  template <class Scalar>
  static void
  computeRoots3(const Scalar m[3][3],
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
    a_over_3 = std::max(a_over_3, static_cast<Scalar>(0.0));
 
    Scalar half_b = 0.5*(c0 + c2_over_3*(2.0*c2_over_3*c2_over_3 - c1));
 
    Scalar q = a_over_3*a_over_3*a_over_3 - half_b*half_b;
    q = std::max(q, static_cast<Scalar>(0.0));
 
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
  static void
  extract_kernel3(const Scalar mat[3][3],
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
  void
  DiagonalizeSym3(const Scalar mat[3][3],
                  Scalar eivals[3], // store eigenvalues here
                  // If the user supplies a 2D eivects array,
                  // then the 3 eigenvectors will be stored in
                  // the 3 rows of the "eivects" 2D array
                  // (ie, eivects[0], eivects[1], eivects[2])
                  Scalar eivects[][3] = NULL,
                  EigenOrderType eival_order = INCREASING_EIVALS)
  {
    const Scalar EPSILON = std::numeric_limits<Scalar>::epsilon();
    // Shift the matrix to the mean eigenvalue and map the matrix
    // coefficients to [-1:1] to avoid over- and underflow.
    Scalar shift = Trace3(mat) / 3.0;
    Scalar scaledMat[3][3];
    Copy3x3(mat, scaledMat);
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
          swap(k, l);
        }
 
        Scalar tmp[3][3];
        Copy3x3(scaledMat, tmp);
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
          Copy3x3(scaledMat, tmp);
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


    // At this point, the eigenvalues should be montonically increasing.

    assert((eivals[0] <= eivals[1]) && (eivals[1] <= eivals[2]));

    if (((eival_order == INCREASING_EIVALS) && (eivals[0] > eivals[2])) ||
        ((eival_order == DECREASING_EIVALS) && (eivals[0] < eivals[2])) ||
        ((eival_order == INCREASING_ABS_EIVALS) && (abs(eivals[0]) > abs(eivals[2]))) ||
        ((eival_order == DECREASING_ABS_EIVALS) && (abs(eivals[0]) < abs(eivals[2]))) ||
        ((eival_order == INCREASINGLY_DISTINCT_EIVALS) && (eivals[1]-eivals[0] > eivals[2]-eivals[1])) ||
        ((eival_order == DECREASINGLY_DISTINCT_EIVALS) && (eivals[1]-eivals[0] < eivals[2]-eivals[1])))
    {
      // Swap the first and last eigenvalues:
      swap(eivals[0], eivals[2]);
      // Swap the first and last eigenvectors:
      for (int d=0; d<3; d++)
        swap(eivects[0][d], eivects[2][d]);
    }

  } //DiagonalizeSym3()




  template <class FlatSymMatrix, class ConstFlatSymMatrix>
  void
  DiagonalizeFlatSym3(ConstFlatSymMatrix source,
                      FlatSymMatrix dest,
                      EigenOrderType eival_order = INCREASING_EIVALS)
  {
    // We need to define some temporary variables that store calculations.
    // These variables will be of type "Scalar" which should be the same
    // numeric type (ie float, double) as the entries in the "dest" argument
    // (which is array-like and supports subscripting).
    // Unfortunately, the next line of hellish C++ failed to infer the type
    // of the entries in the "dest" argument correctly:
    // typedef std::remove_reference<decltype(dest[0])>::type Scalar;
    // ...so I'm just goint to use doubles instead for these temporary variables
    typedef double Scalar;

    // Convert the "source" from type "FlatSymMatrix"
    // ...to an ordinary 3x3 array ("matrix"):
    Scalar matrix[3][3];
    for (int di=0; di<3; di++)
      for (int dj=0; dj<3; dj++)
        matrix[di][dj] = source[ MapIndices_3x3_to_linear[di][dj] ];

    Scalar eivals[3];
    Scalar eivects[3][3];

    // Now diagonalize this 3x3 array
    DiagonalizeSym3(matrix, eivals, eivects, eival_order);

    // Convert to quaternions:
    // There are 3 eigenvectors, each containing 3 numbers (9 total).
    // We will store the eigenvectors and eigenvalues for every voxel
    // in the image in a large (4-dimensional) array.
    // Allocating space for these huge arrays is a serious problem.
    // Rotation matrices can be represented using quaternions.
    // So we convert the eigenvects into quaternion format in order
    // to reduce the space needed from 9=3x3 numbers down to 4 (per voxel)
    // (I suppose I could use Euler angles or Shoemake coordinates,
    //  to reduce this down to 3, but I'm too lazy to be bothered.)

    // Problem:
    // The eigenvectors are orthonormal, but not necessarily a rotation.
    // If the determinant is negative, then flip one of the eigenvectors
    // to insure the determinant is positive.  Handle this below:

    if (Determinant3(eivects) < 0.0) {
      for (int d=0; d<3; d++)
        eivects[0][d] *= -1.0;
    }

    // Note: Each eigenvector is a currently row-vector in eivects[3][3];
    // It's optional, but I prefer to transpose this, because I think of
    // each eigenvector as a column vector.  Either way should work.

    Transpose3(eivects);

    // convert the 3x3 matrix of eigenvector components down to 3 numbers
    // ("Shoemake" coordinates representing the rotation corresponding
    //  to the 3 eigenvectors of the 3x3 matrix.)

    Scalar shoemake[3];
    Matrix2Shoemake(eivects, shoemake);

    // Store the eigenvalues and eigevectors in the following format:
    dest[0] = eivals[0];
    dest[1] = eivals[1];
    dest[2] = eivals[2];
    dest[3] = shoemake[0];
    dest[4] = shoemake[1];
    dest[5] = shoemake[2];

  } // DiagonalizeFlatSym3()



  template <class FlatSymMatrix, class ConstFlatSymMatrix>
  void
  UndiagonalizeFlatSym3(ConstFlatSymMatrix source,
                        FlatSymMatrix dest)
  {
    // We need to define some temporary variables that store calculations.
    // These variables will be of type "Scalar" which should be the same
    // numeric type (ie float, double) as the entries in the "dest" argument
    // (which is array-like and supports subscripting).
    // Unfortunately, the next line of hellish C++ failed to infer the type
    // of the entries in the "dest" argument correctly:
    // typedef std::remove_reference<decltype(dest[0])>::type Scalar;
    // ...so I'm just goint to use doubles instead for these temporary variables
    typedef double Scalar;

    Scalar eivals[3];
    eivals[0] = source[0];
    eivals[1] = source[1];
    eivals[2] = source[2];
    Scalar shoemake[3];
    shoemake[0] = source[3];
    shoemake[1] = source[4];
    shoemake[2] = source[5];
    Scalar eivects[3][3];
    Shoemake2Matrix(shoemake, eivects);

    // Here we reverse the result of an earlier diagonalization
    // Mij = sum_d  lambda[d] * V_d[i] * V_d[j]
    Scalar matrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
    for (int d=0; d<3; d++) {
      for (int di=0; di<3; di++) {
        for (int dj=0; dj<3; dj++) {
          matrix[di][dj] += eivals[d] * eivects[d][di] * eivects[d][dj];
        }
      }
    }
    for (int di=0; di<3; di++) {
      for (int dj=0; dj<3; dj++) {
        dest[ MapIndices_3x3_to_linear[di][dj] ] = matrix[di][dj];
      }
    }

  } // UndiagonalizeFlatSym3()


  template <class Scalar>
  void ConvertFlatSym2Evects3(const Scalar m[6],    //!< a symmetrix 3x3 matrix which has been diagonalized and flattened
                              Scalar eivals[3],     //!< store eigenvalues here
                              Scalar eivects[3][3], //!< store eigenvectors (as rows) here
                              EigenOrderType eival_order = INCREASING_EIVALS //!< choose the order of the eigenvalues and eigenvectors
                              )
  {
    Scalar m_diag[6];
    DiagonalizeFlatSym3(m,
                        m_diag,
                        eival_order);
    ConvertDiagFlatSym2Evects3(m_diag,  //a matrix which has been diagonalized and flattened
                               eivals,
                               eivects);
  }


} //namespace selfadjoint_eigenvalues3


} //namespace visfd


#endif //#ifndef _EIGEN3_SIMPLE_HPP
