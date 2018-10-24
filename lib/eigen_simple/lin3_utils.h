#ifndef _LIN3_UTILS_H
#define _LIN3_UTILS_H

template <class Scalar>
static inline Scalar trace3(const Scalar m[3][3]) {
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


#endif //#ifndef _LIN3_UTILS_H
