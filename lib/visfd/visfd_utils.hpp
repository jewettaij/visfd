///   @file visfd_utils.hpp
///   @brief a collection of functions used repeatedly in different places
///          in the VISFD library.
///   @author Andrew Jewett
///   @date 2019-4-15


#ifndef _VISFD_UTILS_HPP
#define _VISFD_UTILS_HPP

#include <cassert>
#include <limits>
using namespace std;
#include <alloc3d.hpp>


#include <eigen3_simple.hpp>  //defines namespace selfadjoint_eigen3
#include <lin3_utils.hpp>  //defines namespace selfadjoint_eigen3
using namespace visfd::selfadjoint_eigen3;




namespace visfd {




template<typename Scalar>
static inline Scalar SQR(Scalar x) { return x*x; }

template <typename Scalar>
static int SGN(Scalar val) {
  return (static_cast<Scalar>(0) < val) - (val < static_cast<Scalar>(0));
}



/// @brief invert a permutation
template<typename T, typename Integer>
void
invert_permutation(const vector<Integer>& p,
                   vector<T>& p_inv)
{
  p_inv.resize(p.size());
  for (Integer i = 0; i < p.size(); i++)
    p_inv[p[i]] = i;
}



/// @brief apply a permutation to a std::vector in-place
template<typename T, typename Integer>
void
apply_permutation(const vector<Integer>& p,
                  vector<T>& v)
{
  assert(p.size() == v.size());
  int n = p.size();
  vector<T> v_copy(v);
  for (Integer i = 0; i < n; i++) {
    Integer j = p[i];
    assert((0 <= j) && (j < n));
    v[i] = v_copy[j];
  }
}




/// @brief  This function was intended to be used as a way to allow end-users
///         to select certain voxels for consideration.
///         However, because it is difficult to click on exactly the
///         feature you are looking for (which is normally how the
///         "nearby_location" coordinates are obtained),
///         we often need to search for nearby
///         voxels that have the feature we want.
///
///         This function finds the voxel in the image whose indices
///         are closest to the "nearby_location" argument.
///         Only voxels whose corresponding entry in the aaaiVoxels[][][] array
///         belong to "select_these_voxel_types"
///         (AND whose corresponding entry in aaafMask[][][] is not zero)
///         will be considered (unless "invert_selection" is true).
/// @return The function has no return value, however the coordinates of the
///         nearest voxel will be stored in the "nearest_location" argument.
/// @note   Zero-indexing is used.  In other words image indices (ie, entries
///         in the location arguments) are assumed to begin at 0 not 1.
/// @note   In this variant of the function, the "nearby_location" and
///         "location" arguments are both of type C++-style std::array
///         instead of C-style pointers.

template<typename Scalar, typename Label, typename Coordinate>

void
FindNearestVoxel(int const image_size[3],                   //!< #voxels in xyz
                 const array<Coordinate, 3> &nearby_location, //!< find the voxel in aaaiVoxels closest to this
                 array<Coordinate, 3> &nearest_location, //!< and store the location of that voxel here.
                 Label const *const *const *aaaiVoxels, //!< some property associated with each voxel
                 Scalar const *const *const *aaafMask,    //!< optional: Ignore voxels whose mask value is 0
                 set<Label> select_these_voxel_types, //!< voxels must have one of these properties
                 bool invert_selection=false //!< (...or NOT one of these properties)
                 )
{
  nearest_location[0] = -1; // an impossible initial value
  nearest_location[1] = -1; // an impossible initial value
  nearest_location[2] = -1; // an impossible initial value
  Coordinate r_min_sq = -1; //special (uninitialized) impossible value
  // Note: I'm too lazy to write this a faster, smarter way.  It's fast enough.
  for (int iz=0; iz<image_size[2]; iz++) {
    for (int iy=0; iy<image_size[1]; iy++) {
      for (int ix=0; ix<image_size[0]; ix++) {
        if (aaafMask && aaafMask[iz][iy][ix] == 0.0)
          continue;
        array<double, 3> rv;
        rv[0] = nearby_location[0] - ix;
        rv[1] = nearby_location[1] - iy;
        rv[2] = nearby_location[2] - iz;
        double r_sq = SquaredNorm3(rv); // length of rv squared
        bool selected = (select_these_voxel_types.find(aaaiVoxels[iz][iy][ix])
                         != select_these_voxel_types.end());
        if (invert_selection)
          selected = !selected;
        if (! selected)
          continue;
        if ((r_min_sq == -1) || (r_sq < r_min_sq)) {
          r_min_sq = r_sq;
          nearest_location[0] = ix;
          nearest_location[1] = iy;
          nearest_location[2] = iz;
        }
      }
    }
  }
} //FindNearestVoxel()




/// @brief  This function was intended to be used as a way to allow end-users
///         to select certain voxels for consideration.
///         However, because it is difficult to click on exactly the
///         feature you are looking for (which is normally how the
///         "nearby_location" coordinates are obtained),
///         we often need to search for nearby
///         voxels that have the feature we want.
///
///         This function finds the voxel in the image whose indices
///         are closest to the "nearby_location" argument.
///         Only voxels whose corresponding entry in the aaaiVoxels[][][] array
///         belong to "select_these_voxel_types"
///         (AND whose corresponding entry in aaafMask[][][] is not zero)
///         will be considered (unless "invert_selection" is true).
/// @return The function has no return value, however the coordinates of the
///         nearest voxel will be stored in the "nearest_location" argument.
/// @note   Zero-indexing is used.  In other words image indices (ie, entries
///         in the location arguments) are assumed to begin at 0 not 1.
/// @note   In this variant of the function, the "nearby_location" and
///         "location" arguments are both C-style pointers,
///         instead of C++-style std::array.

template<typename Scalar, typename Label, typename Coordinate>

void
FindNearestVoxel(int const image_size[3],             //!< #voxels in xyz
                 const Coordinate nearby_location[3], //!< find voxel in aaaiVoxels closest to this
                 Coordinate nearest_location[3], //!< find voxel in aaaiVoxels closest to this
                 Label const *const *const *aaaiVoxels, //!< some property associated with each voxel
                 Scalar const *const *const *aaafMask,    //!< optional: Ignore voxels whose mask value is 0
                 set<Label> select_these_voxel_types, //!< voxels must have this property
                 bool invert_selection=true //invert selection (skip over them)
                 )
{
  assert(nearby_location);
  assert(nearest_location);

  array<Coordinate, 3> aNearbyLocation;
  array<Coordinate, 3> aNearestLocation;

  for (int d = 0; d < 3; d++)
    aNearbyLocation[d] = nearby_location[d];

  FindNearestVoxel(image_size,
                   aNearbyLocation,
                   aNearestLocation,
                   aaaiVoxels,
                   aaafMask,
                   select_these_voxel_types,
                   invert_selection);

  for (int d = 0; d < 3; d++)
    nearest_location[d] = aNearestLocation[d];
} //FindNearestVoxel()





// @brief  This function chooses the "score" threshold which maximizes the
//         number of correctly labelled training data.
//         (One can think of this threshold as a crude 1D linear SVM.)
//         If theshold_is_lower_bound == true, then data whose scores are above
//         this threshold will be classified as accepted.  (below->rejected.)
//         The reverse is true if theshold_is_lower_bound==false.
//         (This function can be run twice with threshold_is_lower_bound
//          set to both true and false to find both a lower and upper bound.)

template<typename Scalar>
Scalar
ChooseThreshold1D(const vector<Scalar>& training_set_scores, //!< a list of scores in sorted order
                  const vector<bool>& training_set_accepted,//!< was this datum accepted or rejected (positive or negative)
                  bool threshold_is_lower_bound = true //!< should data with scores ABOVE the threshold be accepted? (or BELOW?)
                  )
{
  ptrdiff_t N = training_set_scores.size();
  assert(N == training_set_accepted.size());

  // Create a local copy of the arrays that we can modify:
  vector<Scalar> _training_set_scores = training_set_scores;
  vector<bool> _training_set_accepted  = training_set_accepted;
  
  ptrdiff_t Nn = 0;
  ptrdiff_t Np = 0;
  for (ptrdiff_t i = 0; i < N; i++) {
    if (_training_set_accepted[i])
      Np++;
    else
      Nn++;
  }
  assert(Nn + Np == N);

  Scalar SGN = 1.0;
  if (! threshold_is_lower_bound)
    SGN = -1.0;

  // sort the training set according to their scores:
  vector<tuple<Scalar, ptrdiff_t> > score_index(N);
  for (ptrdiff_t i = 0; i < N; i++)
    score_index[i] = make_tuple(_training_set_scores[i], i);
  if (N > 0) {
    if (threshold_is_lower_bound)
      // then sort the list (of data according to scores) in increasing order
      sort(score_index.begin(),
           score_index.end());
    else
      // then sort the list (of data according to scores) in decreasing order
      sort(score_index.rbegin(),
           score_index.rend());
    vector<ptrdiff_t> permutation(N);
    for (ptrdiff_t i = 0; i < score_index.size(); i++)
      permutation[i] = get<1>(score_index[i]);
    score_index.clear();
    apply_permutation(permutation, _training_set_scores);
    apply_permutation(permutation, _training_set_accepted);
  }

  // The negative training data are assumed to have lower scores
  // than the positive training data (after multiplication by SGN).
  // If the threshold for deciding if a blob is good or bad
  // is set below the scores of all of the negative training data
  // then the number of mistakes is Nn, the number of negative traing data,
  // because if we set the threshold there, every one of them will be
  // predicted to be positive, when in fact they are negative.
  // (Similarly, if the threshold is set above the highest scores,
  // then the number of mistakes is Np, the number of positive traing data.)

  // In the loop below, we start with a threshold below the minimum score, and
  // increase the threshold until we exceed the maximum score.
  // We keep track of which thresholds corresponded to the mininum number of
  // mistakes.

  ptrdiff_t min_num_mistakes;
  min_num_mistakes = Nn;

  {
    ptrdiff_t num_mistakes;
    num_mistakes = Nn;
    int i = -1;
    while (i < N) {
      if (i >= 0) {
        if (_training_set_accepted[i])
          num_mistakes++;
        else
          num_mistakes--;
      }
      if (num_mistakes < min_num_mistakes)
        min_num_mistakes = num_mistakes;
      i++;
    }
  }

  assert(min_num_mistakes <= std::max(Nn, Np));

  // at what threshold value(s) was that same number of mistakes made?
  vector<ptrdiff_t> indices_min_mistakes;

  {
    ptrdiff_t num_mistakes;
    if (threshold_is_lower_bound)
      num_mistakes = Nn;  // (see comment above)
    else
      num_mistakes = Np;  // (see comment above)
    int i = -1;
    while (i < N) {
      if (i >= 0) {
        if (_training_set_accepted[i])
          num_mistakes++;
        else
          num_mistakes--;
      }
      if (num_mistakes == min_num_mistakes)
        indices_min_mistakes.push_back(i);
      i++;
    }
  }
  assert(indices_min_mistakes.size() > 0);

  Scalar threshold; // the threshold returned to the caller

  //then pick the corresponding nearby score 

  //choose the median index into this array of threshold values
  ptrdiff_t i_threshold = indices_min_mistakes[indices_min_mistakes.size()/2];

  if (i_threshold == -1) {
    // ...then you get the minimum number of mistakes
    // if you set the threshold so that all of the training data is accepted.
    // This means the caller only provided examples of "positive" training data
    // (ie. blobs which were accepted.)  
    // Setting the threshold to -infinity will reproduce this result.
    threshold = -SGN * std::numeric_limits<Scalar>::infinity();
  }
  else if (i_threshold == N-1) {
    // ...then you get the minimum number of mistakes
    // if you set the threshold so that all of the training data is rejected.
    // This means the caller only provided examples of "negative" training 
    // data (ie. blobs which were rejected.)
    // Setting the threshold to infinity will reproduce this result.
    threshold = SGN * std::numeric_limits<Scalar>::infinity();
  }
  else {
    // Otherwise choose the threshold which minimizes the number of mistakes:
    threshold = _training_set_scores[i_threshold];
    if (i_threshold < N-1) {
      // for a more robust estimate, 
      // choose the average of this score and the one above it
      threshold = 0.5 * (_training_set_scores[i_threshold] +
                         _training_set_scores[i_threshold+1]);
    }
  }

  return threshold;

} //ChooseThreshold1D()




/// @brief  Calculate the hessian of a 3D image at a particular position 
///         ix, iy, iz, from the differences between neighboring voxels.
///         It is the responsibility of the caller to smooth the source image
///         (if necessary) before this function is called.
/// @note   You must insure that ix,iy,iz (and their neighbors!) lie within
///         the boundaries of the image.  THERE IS NO BOUNDS CHECKING.
template<typename Scalar>
void
CalcHessianFiniteDifferences(Scalar const *const *const *aaafSource, //!< source image
                             int ix, int iy, int iz, //!< location in the image where you wish to calculate the Hessian
                             Scalar (*hessian)[3]  //!<store resulting 3x3 matrixhere
                             )
{
  assert(aaafSource);
  assert(hessian);

  hessian[0][0] = (aaafSource[iz][iy][ix+1] + 
                   aaafSource[iz][iy][ix-1] - 
                   2*aaafSource[iz][iy][ix]);
  hessian[1][1] = (aaafSource[iz][iy+1][ix] + 
                   aaafSource[iz][iy-1][ix] - 
                   2*aaafSource[iz][iy][ix]);
  hessian[2][2] = (aaafSource[iz+1][iy][ix] + 
                   aaafSource[iz-1][iy][ix] - 
                   2*aaafSource[iz][iy][ix]);

  hessian[0][1] = 0.25 * (aaafSource[iz][iy+1][ix+1] + 
                          aaafSource[iz][iy-1][ix-1] - 
                          aaafSource[iz][iy-1][ix+1] - 
                          aaafSource[iz][iy+1][ix-1]);
  hessian[1][0] = hessian[0][1];

  hessian[1][2] = 0.25 * (aaafSource[iz+1][iy+1][ix] + 
                          aaafSource[iz-1][iy-1][ix] - 
                          aaafSource[iz-1][iy+1][ix] - 
                          aaafSource[iz+1][iy-1][ix]);
  hessian[2][1] = hessian[1][2];

  hessian[2][0] = 0.25 * (aaafSource[iz+1][iy][ix+1] + 
                          aaafSource[iz-1][iy][ix-1] - 
                          aaafSource[iz+1][iy][ix-1] - 
                          aaafSource[iz-1][iy][ix+1]);
  hessian[0][2] = hessian[2][0];
} //CalcHessianFiniteDifferences()



/// @brief  Calculate the hessian of a 3D image at a particular position 
///         ix, iy, iz, from the differences between neighboring voxels.
///         It is the responsibility of the caller to smooth the source image
///         (if necessary) before this function is called.
/// @note   This version of the function provides bounds checking for ix,iy,iz.
///         (At the boundaries of the image, it will substitute the voxel
///          brightnesses from the nearest neighbor.)
/// @note   This function was not intended for public use.

template<typename Scalar>
void
CalcHessianFiniteDifferences(Scalar const *const *const *aaafSource, //!< source image
                             int ix, int iy, int iz, //!< location in the image where you wish to calculate the Hessian
                             Scalar (*hessian)[3],  //!< store resulting 3x3 matrixhere
                             const int image_size[3] //!< number of voxels in xyz directions (for bounds checking)
                             )
{
  assert(aaafSource);
  assert(hessian);

  int _ix = ix;
  int _iy = iy;
  int _iz = iz;

  if (image_size) {
    assert(image_size[0] >= 3);
    assert(image_size[1] >= 3);
    assert(image_size[2] >= 3);

    if (_ix == 0)
      _ix++;
    else if (_ix == image_size[0]-1)
      _ix--;

    if (_iy == 0)
      _iy++;
    else if (_iy == image_size[1]-1)
      _iy--;

    if (_iz == 0)
      _iz++;
    else if (_iz == image_size[2]-1)
      _iz--;
  }

  CalcHessianFiniteDifferences(aaafSource,
                               _ix, _iy, _iz,
                               hessian);
}



/// @brief  Calculate the hessian of a 3D image at a particular position 
///         ix, iy, iz, from the differences between neighboring voxels.
///         It is the responsibility of the caller to smooth the source image
///         (if necessary) before this function is called.
/// @note   This version of the function provides bounds checking for ix,iy,iz.
///         (At the boundaries of the image, it will substitute the voxel
///          brightnesses from the nearest neighbor.)
/// @note   This function was not intended for public use.

template<typename Scalar>
void
CalcGradientFiniteDifferences(Scalar const *const *const *aaafSource, //!< source image
                              int ix, int iy, int iz, //!< location in the image where you wish to calculate the Hessian
                              Scalar *gradient,  //!< store resulting 3x3 matrixhere
                              const int image_size[3] //!< number of voxels in xyz directions (for bounds checking)
                              )
{
  assert(aaafSource);
  assert(gradient);

  int _ix = ix;
  int _iy = iy;
  int _iz = iz;

  if (image_size) {
    assert(image_size[0] >= 3);
    assert(image_size[1] >= 3);
    assert(image_size[2] >= 3);

    if (_ix == 0)
      _ix++;
    else if (_ix == image_size[0]-1)
      _ix--;
    if (_iy == 0)
      _iy++;
    else if (_iy == image_size[1]-1)
      _iy--;
    if (_iz == 0)
      _iz++;
    else if (_iz == image_size[2]-1)
      _iz--;
  }

  gradient[0]=0.5*(aaafSource[_iz][_iy][_ix+1] - 
                   aaafSource[_iz][_iy][_ix-1]);
  gradient[1]=0.5*(aaafSource[_iz][_iy+1][_ix] - 
                   aaafSource[_iz][_iy-1][_ix]);
  gradient[2]=0.5*(aaafSource[_iz+1][_iy][_ix] - 
                   aaafSource[_iz-1][_iy][_ix]);
}




/// @brief 
/// Compute a weighted average of the entries in a 3-dimensional array.
///         This function was not intended for public use.
/// @param  array_size  an array of 3 integers storing the size of the array
/// @param  aaafH       the array whose entries will be averaged
/// @param  aaafW       an optional array of weights (for weighted averages)
/// @return the weighted average


template<typename Scalar, typename Integer>

Scalar AverageArr(Integer const array_size[3],
                  Scalar const *const *const *aaafH,
                  Scalar const *const *const *aaafW = nullptr) 
{
  Scalar total = 0.0;
  Scalar denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        Scalar h = aaafH[iz][iy][ix];
        if (aaafW) {
          Scalar w = aaafW[iz][iy][ix];
          h *= w;
          denom += w;
        }
        else
          denom += 1.0;
        total += h;
      }
    }
  }
  return static_cast<Scalar>(total / denom);
} //void AverageArr()



/// @brief 
/// Compute a weighted average of the square of entries in a 3D array.
///         This function was not intended for public use.
/// @param  array_size  an array of 3 integers storing the size of the array
/// @param  aaafH       the array whose entries will be squared and averaged
/// @param  aaafW       an optional array of weights (for weighted averages)
/// @return the weighted average of the squared entries in the array

template<typename Scalar, typename Integer>

static
Scalar _AveSqrArr(Integer const array_size[3],
                  Scalar const *const *const *aaafH,
                  Scalar const *const *const *aaafW = nullptr) 
{
  Scalar total = 0.0;
  Scalar denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        Scalar h = aaafH[iz][iy][ix];
        h *= h;

        if (aaafW) {
          Scalar w = aaafW[iz][iy][ix];
          h *= w;
          denom += w;
        }
        else
          denom += 1.0;
        total += h;
      }
    }
  }
  return static_cast<Scalar>(sqrt(total / denom));
} //void _AveSqrArr()




/// @brief 
/// Compute the weighted standard deviation of entries in a 3-dimensional array.
///         This function was not intended for public use.
/// @param  array_size  an array of 3 integers storing the size of the array
/// @param  aaafH       we will compute the standard deviation of entries
///                     in this array
/// @param  aaafW       an optional array of weights (for weighted averages)
/// @return the weighted standard deviation

template<typename Scalar, typename Integer>

Scalar StdDevArr(Integer const array_size[3],
                 Scalar const *const *const *aaafH,
                 Scalar const *const *const *aaafW = nullptr) 
{
  Scalar ave = AverageArr(array_size, aaafH, aaafW);
  Scalar total = 0.0;
  Scalar denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        Scalar h = (aaafH[iz][iy][ix] - ave);
        h *= h;

        if (aaafW) {
          Scalar w = aaafW[iz][iy][ix];
          h *= w;
          denom += w;
        }
        else
          denom += 1.0;
        total += h;
      }
    }
  }
  return static_cast<Scalar>(sqrt(total / denom));
} //void StdDevArr()



/// @brief Compute a (weighted) sum of the entries in a 3D array.
///         This function was not intended for public use.
/// @param array_size  an array of 3 integers storing the size of the array
/// @param  aaafH       we will compute the sum of the squares of the 
///                     entries in this array
/// @param  aaafW       an optional array of weights (which multiply the
///                     the corresponding entry value)
/// @return the (weighted) sum

template<typename Scalar, typename Integer>

static
Scalar _SumArr(Integer const array_size[3],
               Scalar const *const *const *aaafH,
               Scalar const *const *const *aaafW = nullptr) 
{
  Scalar total = 0.0;
  Scalar denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        Scalar h = aaafH[iz][iy][ix];
        total += h;
        if (aaafW)
          h *= aaafW[iz][iy][ix];
      }
    }
  }
  return static_cast<Scalar>(total);
} //void _SumArr()



/// @brief
/// Compute a (weighted) sum of the squares of the entries in a 3D array
///         This function was not intended for public use.
/// @param array_size  an array of 3 integers storing the size of the array
/// @param  aaafH       we will compute the sum of the squares of the 
///                     entries in this array
/// @param  aaafW       an optional array of weights (which multiply the
///                     the corresponding squared entry value)
/// @return the (weighted) sum of squares

template<typename Scalar, typename Integer>

static
Scalar _SumSqrArr(Integer const array_size[3],
                  Scalar const *const *const *aaafH,
                  Scalar const *const *const *aaafW = nullptr) 
{
  Scalar total = 0.0;
  Scalar denom = 0.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        Scalar h = aaafH[iz][iy][ix];
        h *= h;
        if (aaafW)
          h *= aaafW[iz][iy][ix];
        total += h;
      }
    }
  }
  return static_cast<Scalar>(total);
} //void _SumSqrArr()



/// @brief  Add a scalar offset to all of the entries in a 3D array.
///         This function was not intended for public use.
/// @param  offset       the number that will be added to all array entries
/// @param  array_size  an array of 3 integers storing the size of the array
/// @param  aaafH       the 3D array containing the entries to be modified

template<typename Scalar, typename Integer>

void AddScalarArr(Scalar offset,
                  Integer const array_size[3],
                  Scalar ***aaafH)
{
  for (Integer iz = 0; iz < array_size[2]; iz++)
    for (Integer iy = 0; iy < array_size[1]; iy++)
      for (Integer ix = 0; ix < array_size[0]; ix++)
        aaafH[iz][iy][ix] += offset;
}



/// @brief  Multiply all the entries in a 3D array by a scalar.
///         This function was not intended for public use.
/// @param  scale       the number that will be multiplied by the array entries
/// @param  array_size  an array of 3 integers storing the size of the array
/// @param  aaafH       the 3D array containing the entries to be modified

template<typename Scalar, typename Integer>

void MultiplyScalarArr(Scalar scale,
                       Integer const array_size[3],
                       Scalar ***aaafH)
{
  for (Integer iz = 0; iz < array_size[2]; iz++)
    for (Integer iy = 0; iy < array_size[1]; iy++)
      for (Integer ix = 0; ix < array_size[0]; ix++)
        aaafH[iz][iy][ix] *= scale;
}





/// @brief  Find the minimum entries in a 3D array.
///         This function was not intended for public use.
/// @param  array_size  an array of 3 integers storing the size of the array
/// @param  aaafI       the 3D array containing the entries
/// @param  aaafMask    (optional) If aaafMask[i][j][k]==0 ignore this entry
/// @param  afLocation  (optional) report the location of the minima (ix,iy,iz)
/// @return the minimum entry in the aaafI array (not ignored by the mask)
///         (or std::numeric_limits::infinty(), if the aaafMask array has no non-zero entries)

template<typename Scalar, typename Integer>

static
Scalar _MinArr(Integer const array_size[3],
               Scalar const *const *const *aaafI,
               Scalar const *const *const *aaafMask = nullptr,
               Integer *afLocation = nullptr)
{
  Scalar min_I = std::numeric_limits<Scalar>::infinity(); //return this suspicious value if mask is empty
  bool first_iter = true;
  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        if ((aaafMask) and (aaafMask[iz][iy][ix] == 0.0))
          continue;
        if ((aaafI[iz][iy][ix] < min_I) || (first_iter)) {
          min_I = aaafI[iz][iy][ix];
          if (afLocation) {
            afLocation[0] = ix;
            afLocation[1] = iy;
            afLocation[2] = iz;
          }
        }
        first_iter = false;
      }
    }
  }
  return min_I;
} //void _MinArr()


/// @brief  Find the maximum and maximum entries in a 3D array.
///         This function was not intended for public use.
/// @param  array_size  an array of 3 integers storing the size of the array
/// @param  aaafI       the 3D array containing the entries
/// @param  aaafMask    (optional) If aaafMask[i][j][k]==0 ignore this entry
/// @param  afLocation  (optional) report the location of the maxima (ix,iy,iz)
/// @return the minimum entry in the aaafI array (not ignored by the mask)
///         (or std::numeric_limits::infinty(), if the aaafMask array has no non-zero entries)

template<typename Scalar, typename Integer>

static
Scalar _MaxArr(Integer const array_size[3],
               Scalar const *const *const *aaafI,
               Scalar const *const *const *aaafMask = nullptr,
               Integer *afLocation = nullptr)
{
  Scalar max_I = -std::numeric_limits<Scalar>::infinity(); //return this suspicious value if mask is empty
  bool first_iter = true;
  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        if ((aaafMask) and (aaafMask[iz][iy][ix] == 0.0))
          continue;
        if ((aaafI[iz][iy][ix] < max_I) || (first_iter)) {
          max_I = aaafI[iz][iy][ix];
          if (afLocation) {
            afLocation[0] = ix;
            afLocation[1] = iy;
            afLocation[2] = iz;
          }
        }
        first_iter = false;
      }
    }
  }
  return max_I;
} //void _MaxArr()



/// @brief  Generate a histogram of the values contained in a 3D numeric array
/// @param  paHistX     store entry values corresponding to the center of each bin here
/// @param  paHistY     store the number of entries that fall into that bin here
/// @param  nbins       specify number of bins (if positive)
/// @param  bin_width   alternatively, (if nbins<=0), specify the width of each bin
/// @param  aaafI       the 3D array containing the entries
/// @param  aaafMask    (optional) If aaafMask[i][j][k]==0 ignore this entry

template<typename Scalar, typename Integer>

void
HistogramArr(Scalar **paHistX,  
             size_t  **paHistY,
             Integer &nbins,
             Scalar &bin_width,
             Integer const array_size[3],
             Scalar ***aaafI,
             Scalar const *const *const *aaafMask = nullptr) 
{
  Scalar hmin = 0.0;
  Scalar hmax = -1.0;

  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        if (hmin > hmax) {
          hmin = aaafI[iz][iy][ix];
          hmax = aaafI[iz][iy][ix];
        }
        if (aaafI[iz][iy][ix] < hmin)
          hmin = aaafI[iz][iy][ix];
        if (aaafI[iz][iy][ix] > hmax)
          hmax = aaafI[iz][iy][ix];
      }
    }
  }

  if (nbins > 1) {
    bin_width = (hmax - hmin) / (nbins - 1);
  }
  else if (bin_width > 0) {
    nbins = 1 + ceil((hmax - hmin)/bin_width);
  }

  assert(paHistX);
  *paHistX = new Scalar [nbins];
  for (int i=0; i<nbins; i++)
    (*paHistX)[i] = hmin + i*bin_width;

  assert(paHistY);
  *paHistY = new size_t [nbins];
  for (int i=0; i < nbins; i++)
    (*paHistY)[i] = 0;
  for (Integer iz = 0; iz < array_size[2]; iz++) {
    for (Integer iy = 0; iy < array_size[1]; iy++) {
      for (Integer ix = 0; ix < array_size[0]; ix++) {
        Scalar h = aaafI[iz][iy][ix];
        if ((! aaafMask) || (aaafMask[iz][iy][ix] != 0.0)) {
          Scalar h = aaafI[iz][iy][ix];
          Integer i = floor( (h - hmin)/bin_width );
          assert((0 <= i) && (i < nbins));
          (*paHistY)[i]++;
        }
      }
    }
  }
} //void HistogramArr()



} //namespace visfd



#endif //#ifndef _VISFD_UTILS_HPP
