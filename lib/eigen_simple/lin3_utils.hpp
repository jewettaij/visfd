#ifndef _LIN3_UTILS_HPP
#define _LIN3_UTILS_HPP

#include <cmath>
#include <cassert>
using namespace std;


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
  Scalar L_inv = std::sqrt(SquaredNorm3(a));
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
  if (source != dest) {
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        dest[i][j] = source[j][i];
  }
  else {
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        swap(dest[i][j], dest[j][i]);
  }
}


template <class Scalar>
static inline void Transpose3(Scalar m[3][3]) {
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      swap(m[i][j], m[j][i]);
}


template <class Scalar>
static inline Scalar Determinant3(const Scalar m[3][3]) {
  Scalar v12[3];
  CrossProduct(m[0], m[1], v12);
  return DotProduct3(m[2], v12);
}


/// @brief  Convert a 3x3 rotation matrix (M) to a quaternion (q)
/// http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/

template<class Scalar>
static inline void Matrix2Quaternion(const Scalar M[3][3], Scalar q[4])
{
  Scalar S;
  Scalar qw, qx, qy, qz;
  Scalar tr = Trace3(M);  // = M[0][0] + M[1][1] + M[2][2];
  if (tr > 0) {
    S = std::sqrt(tr+1.0) * 2;                        // S=4*qw 
    qw = 0.25 * S;
    qx = (M[2][1] - M[1][2]) / S;
    qy = (M[0][2] - M[2][0]) / S;
    qz = (M[1][0] - M[0][1]) / S;
  }
  else if ((M[0][0] > M[1][1]) and (M[0][0] > M[2][2])) {
    S = std::sqrt(1.0 + M[0][0] - M[1][1] - M[2][2]) * 2;   // S=4*qx 
    qw = (M[2][1] - M[1][2]) / S;
    qx = 0.25 * S;
    qy = (M[0][1] + M[1][0]) / S;
    qz = (M[0][2] + M[2][0]) / S;
  }
  else if (M[1][1] > M[2][2]) {
    S = std::sqrt(1.0 + M[1][1] - M[0][0] - M[2][2]) * 2;   // S=4*qy
    qw = (M[0][2] - M[2][0]) / S;
    qx = (M[0][1] + M[1][0]) / S;
    qy = 0.25 * S;
    qz = (M[1][2] + M[2][1]) / S;
  }
  else {
    S = std::sqrt(1.0 + M[2][2] - M[0][0] - M[1][1]) * 2;   // S=4*qz
    qw = (M[1][0] - M[0][1]) / S;
    qx = (M[0][2] + M[2][0]) / S;
    qy = (M[1][2] + M[2][1]) / S;
    qz = 0.25 * S;
  }
  q[0] = qw;
  q[1] = qx;
  q[2] = qy;
  q[3] = qz;

} //Matrix2Quaternion(const Scalar M[3][3], Scalar q[4])






/// @brief  Convert a quaternion (q) to a 3x3 rotation matrix (M)
/// http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm

template<class RealNum>
static inline void Quaternion2Matrix(const RealNum q[4], RealNum M[3][3])
{
  //M[0][0] =  (q[0]*q[0])-(q[1]*q[1])-(q[2]*q[2])+(q[3]*q[3]);
  //M[1][1] = -(q[0]*q[0])+(q[1]*q[1])-(q[2]*q[2])+(q[3]*q[3]);
  //M[2][2] = -(q[0]*q[0])-(q[1]*q[1])+(q[2]*q[2])+(q[3]*q[3]);
  //M[0][1] = 2*(q[0]*q[1] - q[2]*q[3]);
  //M[1][0] = 2*(q[0]*q[1] + q[2]*q[3]);
  //M[1][2] = 2*(q[1]*q[2] - q[0]*q[3]);
  //M[2][1] = 2*(q[1]*q[2] + q[0]*q[3]);
  //M[0][2] = 2*(q[0]*q[2] + q[1]*q[3]);
  //M[2][0] = 2*(q[0]*q[2] - q[1]*q[3]);

  M[0][0] =  1.0 - 2*(q[2]*q[2]) - 2*(q[3]*q[3]);
  M[1][1] =  1.0 - 2*(q[1]*q[1]) - 2*(q[3]*q[3]);
  M[2][2] =  1.0 - 2*(q[1]*q[1]) - 2*(q[2]*q[2]);
  M[0][1] = 2*(q[1]*q[2] - q[3]*q[0]);
  M[1][0] = 2*(q[1]*q[2] + q[3]*q[0]);
  M[1][2] = 2*(q[2]*q[3] - q[1]*q[0]);
  M[2][1] = 2*(q[2]*q[3] + q[1]*q[0]);
  M[0][2] = 2*(q[1]*q[3] + q[2]*q[0]);
  M[2][0] = 2*(q[1]*q[3] - q[2]*q[0]);
}


/// @brief  Convert a 3-component Shoemake coordinate list to a quaternion (q)
///         Shoemake, Graphics Gems III (1992) pp. 124-132
template <class Scalar>
static inline void Shoemake2Quaternion(const Scalar sm[3], Scalar q[4])
{
  const Scalar M_2PI = 6.283185307179586;
  Scalar X0 = sm[0];
  Scalar X1 = sm[1];
  Scalar X2 = sm[2];
  Scalar theta1 = M_2PI * X1;
  Scalar theta2 = M_2PI * X2;
  Scalar r1 = std::sqrt(1.0-X0);
  Scalar r2 = std::sqrt(X0);
  Scalar s1 = std::sin(theta1);
  Scalar c1 = std::cos(theta1);
  Scalar s2 = std::sin(theta2);
  Scalar c2 = std::cos(theta2);

  // Alternative quaternion convention, where q[3] is real instead of q[0]
  q[0] = s1*r1;
  q[1] = c1*r1;
  q[2] = s2*r2;
  q[3] = c2*r2;

  // Alternative quaternion convention, where q[0] is real instead of q[3]
  //q[3] = s1*r1;
  //q[0] = c1*r1;
  //q[1] = s2*r2;
  //q[2] = c2*r2;
}



/// @brief  Convert a quaternion (q) to a 3-component Shoemake coordinate list
///         Shoemake, Graphics Gems III (1992) pp. 124-132
template <class Scalar>
static inline void Quaternion2Shoemake(const Scalar q[4], Scalar sm[3])
{
  const Scalar M_2PI = 6.283185307179586;

  // Alternative quaternion convention, where q[3] is real instead of q[0]
  Scalar r1 = std::sqrt(q[0]*q[0] + q[1]*q[1]);
  Scalar r2 = std::sqrt(q[2]*q[2] + q[3]*q[3]);
  Scalar X0 = r2*r2;
  Scalar theta1 = 0.0;
  if (r1 > 0)
    theta1 = std::atan2(q[0], q[1]);
  Scalar theta2 = 0.0;
  if (r2 > 0)
    theta2 = std::atan2(q[2], q[3]);

  // Alternative quaternion convention, where q[0] is real instead of q[3]
  //Scalar r1 = std::sqrt(q[3]*q[3] + q[0]*q[0]);
  //Scalar r2 = std::sqrt(q[1]*q[1] + q[2]*q[2]);
  //Scalar X0 = r2;
  //Scalar theta1 = 0.0;
  //if (r1 > 0)
  //  theta1 = std::atan2(q[3], q[0]);
  //Scalar theta2 = 0.0;
  //if (r2 > 0)
  //  theta2 = std::atan2(q[1], q[2]);

  Scalar X1 = theta1 / M_2PI;
  Scalar X2 = theta2 / M_2PI;
  sm[0] = X0;
  sm[1] = X1;
  sm[2] = X2;
}



/// @brief  Convert 3 Shoemake coordinates to a 3x3 rotation matrix (M)
template <class Scalar>
static inline void Shoemake2Matrix(const Scalar sm[3], Scalar M[3][3]) {
  Scalar q[4];
  Shoemake2Quaternion(sm, q);
  Quaternion2Matrix(q, M);
}


/// @brief  Convert 3 Shoemake coordinates to a 3x3 rotation matrix (M)
template <class Scalar>
static inline void Matrix2Shoemake(const Scalar M[3][3], Scalar sm[3]) {
  Scalar q[4];
  Matrix2Quaternion(M, q);
  Quaternion2Shoemake(q, sm);
}



/// Lookup tables for representing symmetric 3x3 matrices as compact 1D tables

static const int MapIndices_3x3_to_linear[][3] = 
  { {0, 3, 5},
    {3, 1, 4},
    {5, 4, 2} };

static const int MapIndices_linear_to_3x3[][2] =
  { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };


/// @brief  A simple class for storing 3x3 symmetric arrays
///         A symmetric 3x3 matrix contains only 6 unique numbers.
///         Internally, this is stored as a 1-D array with 6 entries.

template<class Scalar>

class SymmetricMatrix3x3 {

private:
  Scalar data[6];

public:

  /// @brief  Access the contents of this object as if it were a 2-dimensional
  ///         matrix.
  /// @returns  data[i], where i is 1-to-1 with (di,dj), assuming di<=dj

  Scalar& operator() (int di, int dj) {
    assert((0 <= di < 3) && (0 <= dj < 3));
    return data [ MapIndices_3x3_to_linear[di][dj] ];
  }

  /// @brief  Sometimes it's useful to be able to pass a "SymmetricMatrix3x3"
  ///         object to functions which expect an object which can be
  ///         subscripted as a single array (in this case, with 6 entries).
  ///         The lookup functions defined below can be used to make sense of
  ///         the order of the entries in this 1D table.
  /// @returns  data[i]

  Scalar& operator[](int i) {
    assert(0 <= i < 6);
    return data[i];
  }



  /// @brief  A symmetric 3x3 matrix contains only 6 unique numbers.
  ///         To save memory, save these 6 numbers in a 1-dimensional array
  ///         and define a 1-to-1 mapping between the 1-D and 3-D arrays.
  ///         These maps are made public so that the caller can convert
  ///         between one representation and the other.
  ///         This first function, converts di, dj -> i,
  ///         where di, and dj are indices into a symmetric 3x3 array,
  ///         and "i" is the index into the corresponding 1-D array.
  /// @code
  /// inputs:  di ∈ [0,2]
  ///          dj ∈ [0,2]
  /// @endcode
  /// @return:  i ∈ [0,5]
  
  const int MapIndices(int di, int dj)
  { return MapIndices_3x3_to_linear[di][dj]; }



  /// @brief  The second function is the inverse lookup.
  ///        convert i -> di, dj, where "i" is the index into the 1-D array, and
  ///        di, and dj are indices into the corresponding symmetric 3x3 array.
  /// 
  /// input:    i ∈ [0,5]
  /// outputs: di ∈ [0,2]
  ///          dj ∈ [0,2]
  /// @returns   a pointer to a (size 2) array containing {di, dj}, where di<=dj
  ///            (this array is statically allocated so you don't have
  ///            to invoke delete[] on it later)
  static const int *MapIndices(int i)
  { return MapIndices_linear_to_3x3[i]; }


}; //class SymmetricMatrix3x3





#endif //#ifndef _LIN3_UTILS_HPP
