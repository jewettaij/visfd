#ifndef _LIN3_UTILS_H
#define _LIN3_UTILS_H

template <class Scalar>
static inline Scalar Trace3(const Scalar m[3][3]) {
  return m[0][0]+m[1][1]+m[2][2];
}

template <class Scalar>
static inline void CrossProduct(const Scalar a[3],
                                const Scalar b[3],
                                Scalar dest[3])
{
  dest[2] = a[0]*b[1] - a[1]*b[0];
  // ...and cyclic permutations...
  dest[0] = a[1]*b[2] - a[2]*b[1];
  dest[1] = a[2]*b[0] - a[0]*b[2];
}

template <class Scalar>
static inline Scalar DotProduct3(const Scalar a[3], const Scalar b[3]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


template <class Scalar>
static inline Scalar SquaredNorm3(const Scalar a[3]) {
  return DotProduct3(a, a);
}

template <class Scalar>
static inline void Normalize3(Scalar a[3]) {
  Scalar L_inv = sqrt(SquaredNorm3(a));
  if (L_inv > 0.0) {
    L_inv = 1.0 / L_inv;
    for (int d=0; d<3; d++)
      a[d] *= L_inv;
  }
  else {
    a[0] = 1.0;
    a[1] = 0.0;
    a[2] = 0.0;
  }
}

template <class Scalar>
static inline void Copy3(const Scalar source[3][3], Scalar dest[3][3]) {
  //for (int i=0; i<3; i++)
  //  for (int j=0; j<3; j++)
  //    dest[i][j] = source[i][j];
  // Use memcpy() instead:
  memcpy(dest, source,  3*3*sizeof(Scalar));
}

template <class Scalar>
static inline void Transpose3(const Scalar source[3][3], Scalar dest[3][3]) {
  assert(source != dest);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      dest[i][j] = source[j][i];
}

template <class Scalar>
static inline void Transpose3(Scalar m[3][3]) {
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      Scalar tmp = m[i][j];
      m[i][j] = m[j][i];
      m[j][i] = tmp;
    }
  }
}


template <class Scalar>
static inline Scalar Determinant3(const Scalar m[3][3]) {
  Scalar v12[3];
  CrossProduct(m[0], m[1], v12);
  return DotProduct3(m[2], v12);
}


/// @brief  Convert a 3x3 rotation matrix (M) to a quaternion (q)
/// http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/

template<class RealNum>
static inline void Matrix2Quaternion(const RealNum M[3][3], RealNum q[4])
{
  RealNum S;
  RealNum qw, qx, qy, qz;
  RealNum tr = Trace3(M);  // = M[0][0] + M[1][1] + M[2][2];
  if (tr > 0) {
    S = sqrt(tr+1.0) * 2;                        // S=4*qw 
    qw = 0.25 * S;
    qx = (M[2][1] - M[1][2]) / S;
    qy = (M[0][2] - M[2][0]) / S;
    qz = (M[1][0] - M[0][1]) / S;
  }
  else if ((M[0][0] > M[1][1]) and (M[0][0] > M[2][2])) {
    S = sqrt(1.0 + M[0][0] - M[1][1] - M[2][2]) * 2;   // S=4*qx 
    qw = (M[2][1] - M[1][2]) / S;
    qx = 0.25 * S;
    qy = (M[0][1] + M[1][0]) / S;
    qz = (M[0][2] + M[2][0]) / S;
  }
  else if (M[1][1] > M[2][2]) {
    S = sqrt(1.0 + M[1][1] - M[0][0] - M[2][2]) * 2;   // S=4*qy
    qw = (M[0][2] - M[2][0]) / S;
    qx = (M[0][1] + M[1][0]) / S;
    qy = 0.25 * S;
    qz = (M[1][2] + M[2][1]) / S;
  }
  else {
    S = sqrt(1.0 + M[2][2] - M[0][0] - M[1][1]) * 2;   // S=4*qz
    qw = (M[1][0] - M[0][1]) / S;
    qx = (M[0][2] + M[2][0]) / S;
    qy = (M[1][2] + M[2][1]) / S;
    qz = 0.25 * S;
  }
  q[0] = qw;
  q[1] = qx;
  q[2] = qy;
  q[3] = qz;

} //Matrix2Quaternion(const RealNum M[3][3], RealNum q[4])



/// @brief  Convert a quaternion (q) to a 3x3 rotation matrix (M)
template<class RealNum>
static inline void Quaternion2Matrix(const RealNum q[4], RealNum M[3][3])
{
  M[0][0] =  (q[0]*q[0])-(q[1]*q[1])-(q[2]*q[2])+(q[3]*q[3]);
  M[1][1] = -(q[0]*q[0])+(q[1]*q[1])-(q[2]*q[2])+(q[3]*q[3]);
  M[2][2] = -(q[0]*q[0])-(q[1]*q[1])+(q[2]*q[2])+(q[3]*q[3]);
  M[0][1] = 2*(q[0]*q[1] - q[2]*q[3]);
  M[1][0] = 2*(q[0]*q[1] + q[2]*q[3]);
  M[1][2] = 2*(q[1]*q[2] - q[0]*q[3]);
  M[2][1] = 2*(q[1]*q[2] + q[0]*q[3]);
  M[0][2] = 2*(q[0]*q[2] + q[1]*q[3]);
  M[2][0] = 2*(q[0]*q[2] - q[1]*q[3]);
}




#endif //#ifndef _LIN3_UTILS_H
