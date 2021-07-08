///   @file feature.hpp
///   @brief  IMPLEMENTATION DETAILS.  THIS CODE IS NOT INTENDED FOR PUBLIC USE
///           (The functions defined here are invoked by "feature.hpp")
///   @author Andrew Jewett
///   @date 2019-4-17

#ifndef _FEATURE_IMPLEMENTATION_HPP
#define _FEATURE_IMPLEMENTATION_HPP

#include <cstring>
#include <cassert>
#include <cmath>
#include <limits>
#include <ostream>
#include <vector>
#include <tuple>
#include <set>
#include <queue>
using namespace std;
#include <err_visfd.hpp> // defines the "VisfdErr" exception type
#include <eigen3_simple.hpp>  // defines matrix diagonalizer (DiagonalizeSym3())
#include <visfd_utils.hpp>    // defines invert_permutation(), AveArray(), ...
#include <alloc2d.hpp>    // defines Alloc2D() and Dealloc2D()
#include <alloc3d.hpp>    // defines Alloc3D() and Dealloc3D()
#include <filter1d.hpp>   // defines "Filter1D" (used in ApplySeparable())
#include <filter3d.hpp>   // defines common 3D image filters



namespace visfd {



/// @brief Figure out the score for blobs located at each position in crds[].
///        Terminology:  Every blob has a position, diameter, and a score.
///        The j'th "blob" is a sphere, centerd at blob_crds[j], with diameter
///        blob_diameters[j].  Associated with every blob is a "score"
///        (a number) stored in blob_scores[j].
///        If crds[i] lies inside one of the blobs, store the corresponding
///        score for that blob in scores[i].  (If crds[i] lies inside more
///        than one sphere, give priority spheres occuring later in the list.)
///        If crds[i] lies outside any of the spheres, store -infinity there.
///
/// @returns void.  Results are stored in "scores".
///
/// @note  THIS FUNCTION WAS NOT INTENDED FOR PUBLIC USE

template<typename Scalar, typename Integer>
static void
_FindBlobScores(const vector<array<Scalar,3> >& crds, //!< locations of blob-like things we are looking for
                vector<Scalar>& scores, //!< stores the score of the blob (sphere) contains that position
                vector<Integer>& sphere_ids, //!< which blob (sphere) contains this position?
                const vector<array<Scalar,3> >& blob_crds, //!< location of the center of each spherical blob (sorted in order of increasing priority)
                const vector<Scalar>& blob_diameters,  //!< diameter of each blob
                const vector<Scalar>& blob_scores, //!< "score" of each blob (a number)
                SortCriteria sort_blob_criteria = SORT_DECREASING_MAGNITUDE, //!< give priority to high or low scoring blobs?
                ostream *pReportProgress = nullptr //!< report progress back to the user?
                )
{
  // Sort the blobs by score in order of increasing priority
  vector<array<Scalar,3> > blob_crds_sorted = blob_crds;
  vector<Scalar> blob_diameters_sorted = blob_diameters;
  vector<Scalar> blob_scores_sorted = blob_scores;

  SortBlobs(blob_crds_sorted,
            blob_diameters_sorted, 
            blob_scores_sorted,
            sort_blob_criteria,
            true,
            nullptr,
            pReportProgress);

  // Figure out which sphere is located at each position
  FindSpheres(crds,
              sphere_ids,
              blob_crds_sorted,
              blob_diameters_sorted,
              pReportProgress);

  assert(sphere_ids.size() == crds.size());
  scores.resize(sphere_ids.size(),
                -std::numeric_limits<Scalar>::infinity()); //<-impossible score
  const size_t UNOCCUPIED = 0;
  for (size_t i=0; i < sphere_ids.size(); i++) {
    size_t which_blob = sphere_ids[i];
    if (which_blob != UNOCCUPIED)
      scores[i] = blob_scores_sorted[which_blob-1];
  }
} //_FindBlobScores()






/// @brief  Throw an exception if either Nn or Np = 0.
///         Include a (verbose) error message.
template<typename Integer>
static void
_ComplainIfTrainingDataEmpty(Integer Nn, Integer Np) {
  if (Nn == 0)
    throw VisfdErr("Error: Empty list of negative training examples.\n"
                   "       Either you have provided no examples of data that you want to discard\n"
                   "          (IE. Your file of negative training examples is empty),\n"
                   "       OR none of the examples that you provided are sufficiently close to ANY\n"
                   "       of the blobs in the list of blobs that you are considering discarding.\n"
                   "       Either provide more negative training examples, or specify your\n"
                   "       selection threshold(s) manually.\n");
  if (Np == 0)
    throw VisfdErr("Error: Empty list of positive training examples.\n"
                   "       Either you have provided no examples of data that you want to keep\n"
                   "          (IE. Your file of positive training examples is empty),\n"
                   "       OR none of the examples that you provided are sufficiently close to ANY\n"
                   "       of the blobs in the list of blobs that you are considering discarding.\n"
                   "       Either provide more positive training examples, or specify your\n"
                   "       selection threshold(s) manually.\n");
} //_ComplainIfTrainingDataEmpty(Integer Nn, Integer Np)





/// @brief find either an upper or a lower bound for a list
///        of 1-D features (scores).
///        (Although both upper and lower bounds are calculated, 
///         only one of them should have a non-infinite value.
///         Otherwise your data is not one-sided.)
///
/// @return This function returns void.
///         Assuming they are not nullptr, the results are stored in:
///           *pthreshold_lower_bound
///           *pthreshold_uppwer_bound
///
/// @note THIS FUNCTION WAS NOT INTENDED FOR PUBLIC USE

template<typename Scalar>

static void
_ChooseThresholdInterval(const vector<Scalar>& training_scores, //!< a 1-D feature describing each training data
                         const vector<bool>& training_accepted, //!< classify each training data as "accepted" (true) or "rejected" (false)
                         Scalar *pthreshold_lower_bound = nullptr,  //!< return threshold to the caller (if not nullptr)
                         Scalar *pthreshold_upper_bound = nullptr,  //!< return threshold to the caller (if not nullptr)
                         ostream *pReportProgress = nullptr //!< report progress back to the user?
                         )
{
  size_t N = training_scores.size();
  assert(N == training_accepted.size());

  Scalar threshold_lower_bound = -1.0;
  Scalar threshold_upper_bound = -1.0;

  { // calculate threshold_lower_bound, threshold_upper_bound

    // Choose the upper and lower bounds
    // I assumed that the accepted results lie between
    //    thresh_lower_bound  and  thresh_upper_bound
    // (This might be a bad assumption.  It could be the inverse of this.)
    // Another issue:
    //Sometimes if we choose the lower bound first, we get the wrong upper bound
    //Sometimes if we choose the upper bound first, we get the wrong lower bound
    //So we try both ways, and pick which way minimizes the number of mistakes.
    //(Note: This still might not give us both optimal upper and lower bounds.
    //       However if either lower_bound=-infinity OR upper_bound=infinity,
    //       this will give us the optimal threshold for the other boundary.)
    size_t num_mistakes_lower_bound_first;
    Scalar choose_threshold_lower_bound_first;
    Scalar choose_threshold_upper_bound_second;
    {
      // choose the lower bound first:
      choose_threshold_lower_bound_first =
        ChooseThreshold1D(training_scores,
                          training_accepted,
                          true);
      vector<bool> training_accepted_remaining;
      vector<Scalar> training_scores_remaining;
      for (size_t i = 0; i < N; i++) {
        if (training_scores[i] >= choose_threshold_lower_bound_first) {
          training_scores_remaining.push_back(training_scores[i]);
          training_accepted_remaining.push_back(training_accepted[i]);
        }
      }
      choose_threshold_upper_bound_second =
        ChooseThreshold1D(training_scores_remaining,
                          training_accepted_remaining,
                          false);
      num_mistakes_lower_bound_first = 0;
      for (size_t i = 0; i < N; i++) {
        if (training_accepted[i] !=
            ((training_scores[i] >= choose_threshold_lower_bound_first) &&
             (training_scores[i] <= choose_threshold_upper_bound_second)))
          num_mistakes_lower_bound_first++;
      }
    } // calculate choose_threshold_lower_bound_first, choose_threshold_upper_bound_second

    size_t num_mistakes_upper_bound_first;
    Scalar choose_threshold_upper_bound_first;
    Scalar choose_threshold_lower_bound_second;
    {
      // choose the upper bound first:
      choose_threshold_upper_bound_first =
        ChooseThreshold1D(training_scores,
                          training_accepted,
                          false);
      vector<bool> training_accepted_remaining;
      vector<Scalar> training_scores_remaining;
      for (size_t i = 0; i < N; i++) {
        if (training_scores[i] <= choose_threshold_upper_bound_first) {
          training_scores_remaining.push_back(training_scores[i]);
          training_accepted_remaining.push_back(training_accepted[i]);
        }
      }
      choose_threshold_lower_bound_second =
        ChooseThreshold1D(training_scores_remaining,
                          training_accepted_remaining,
                          true);
      num_mistakes_upper_bound_first = 0;
      for (size_t i = 0; i < N; i++) {
        if (training_accepted[i] !=
            ((training_scores[i] >= choose_threshold_lower_bound_second) &&
             (training_scores[i] <= choose_threshold_upper_bound_first)))
          num_mistakes_upper_bound_first++;
      }
    } // calculate choose_threshold_upper_bound_first, choose_threshold_lower_bound_second

    if (num_mistakes_lower_bound_first <= num_mistakes_upper_bound_first) {
      threshold_lower_bound = choose_threshold_lower_bound_first;
      threshold_upper_bound = choose_threshold_upper_bound_second;
    }
    else {
      threshold_lower_bound = choose_threshold_lower_bound_second;
      threshold_upper_bound = choose_threshold_upper_bound_first;
    }
  } // calculate threshold_lower_bound, threshold_upper_bound

  if (pthreshold_lower_bound)
    *pthreshold_lower_bound = threshold_lower_bound;
  if (pthreshold_upper_bound)
    *pthreshold_upper_bound = threshold_upper_bound;

  if (pReportProgress) {
    *pReportProgress
      << "  threshold lower bound: " << threshold_lower_bound << "\n"
      << "  threshold upper bound: " << threshold_upper_bound << "\n";
    // Optional: Also report the number of mistakes to the user
    size_t num_false_negatives = 0;
    size_t num_false_positives = 0;
    for (size_t i=0; i < N; i++) {
      if ((training_scores[i] >= threshold_lower_bound) &&
          (training_scores[i] <= threshold_upper_bound)) {
        if (! training_accepted[i])
          num_false_positives++;
      }
      else {
        if (training_accepted[i])
          num_false_negatives++;
      }
    }

    size_t Nn = 0;
    size_t Np = 0;
    for (size_t i=0; i < N; i++) {
      if (training_accepted[i])
        Np++;
      else
        Nn++;
    }
    *pReportProgress
      << "  number of false positives: " << num_false_positives
      << " (out of " << Nn << " negatives)\n"
      << "  number of false negatives: " << num_false_negatives
      << " (out of " << Np << " positives)\n"
      << endl;
  } //if (pReportProgress)

} //_ChooseThresholdInterval()





/// @brief  This function calculates a threshold score above which (or below
///         which) the training data is usually accepted (or rejected).
///         The threshold is chosen to maximize the accuracy of training
///         data provided by the caller.
///         (It minimizes the number of times that either the positive training
///          set has a score outside this range, and the negative training set
///          has scores inside this range.  Equal weight is given to to
///          either false positives or false negatives.)
///
/// @note   THIS FUNCTION WAS NOT INTENDED FOR PUBLIC USE.
///
/// @note:  This function was intended to be used when it is possible to use
///         a single score threshold to distinguish good blobs from bad ones.
///         It was not not intended to be used if there is a narrow interval
///         of good scores.  (IE having both an upper and a lower bound.)
///         Although this function calculates both upper and lower bounds
///         for the score and returns them to the caller, ...
///            (*pthreshold_lower_bound and *pthreshold_lower_bound)
///         ...usually, only one of them is set to a meaningful value.
///         The other threshold should be set to either -infinity or +infinity.
///         If this is not the case, then your training data does not
///         fit the assumptions used by this function, and you should
///         discard the results.
///
/// @return This function does not return anything.
///         The two thresholds are returned to the caller using the 
///         pthreshold_lower_bound and pthreshold_upper_bound arguments.

template<typename Scalar>
static void
_ChooseBlobScoreThresholds(const vector<array<Scalar,3> >& blob_crds, //!< location of each blob (in voxels, sorted by score in increasing priority)
                           const vector<Scalar>& blob_diameters,  //!< diameger of each blob (sorted by score in increasing priority)
                           const vector<Scalar>& blob_scores, //!< priority of each blob (sorted by score in increasing priority)
                           const vector<array<Scalar,3> >& training_crds, //!< locations of blob-like things
                           const vector<bool>& training_accepted, //!< classify each blob as "accepted" (true) or "rejected" (false)
                           Scalar *pthreshold_lower_bound = nullptr, //!< return threshold to the caller
                           Scalar *pthreshold_upper_bound = nullptr, //!< return threshold to the caller
                           SortCriteria sort_blob_criteria = SORT_DECREASING_MAGNITUDE, //!< give priority to high or low scoring blobs?
                           ostream *pReportProgress = nullptr //!< report progress back to the user?
                           )
{
  vector<array<Scalar,3> > final_training_crds;
  vector<bool> final_training_accepted;
  vector<Scalar> final_training_scores;

  FindBlobScores(training_crds,
                 training_accepted,
                 final_training_crds,
                 final_training_accepted,
                 final_training_scores,
                 blob_crds,
                 blob_diameters,
                 blob_scores,
                 sort_blob_criteria,
                 pReportProgress);
  size_t N = final_training_crds.size();
  assert(N == final_training_scores.size());
  assert(N == final_training_accepted.size());

  // make sure that both positive and negative training data sets are non-empty
  size_t Nn = 0;
  size_t Np = 0;
  for (size_t i=0; i < N; i++) {
    if (final_training_accepted[i])
      Np++;
    else
      Nn++;
  }
  _ComplainIfTrainingDataEmpty(Nn, Np);
  
  if (pReportProgress)
    *pReportProgress
      << "  examining training data to determine optimal thresholds\n";

  _ChooseThresholdInterval(final_training_scores,
                           final_training_accepted,
                           pthreshold_lower_bound,
                           pthreshold_upper_bound,
                           pReportProgress);
                             
} //_ChooseBlobScoreThresholds()




/// @brief  The following version of the function works with multiple
///         independent training sets.  (That is, training sets corresponding
///         to independent sets of overlapping blob_crds, blob_diameters, ...
///         In practice, these different training sets and blob lists are taken
///         from different images.  The size of the following arguments should
///         equal this number of images:   blob_crds, blob_diameters,
///         blob_scores, training_crds, training_accepted
///
/// @note   THIS FUNCTION WAS NOT INTENDED FOR PUBLIC USE.

template<typename Scalar>
static void
_ChooseBlobScoreThresholdsMulti(const vector<vector<array<Scalar,3> > >& blob_crds, //!< location of each blob (in voxels, sorted by score in increasing priority)
                                const vector<vector<Scalar> >& blob_diameters,  //!< diameger of each blob (sorted by score in increasing priority)
                                const vector<vector<Scalar> >& blob_scores, //!< priority of each blob (sorted by score in increasing priority)
                                const vector<vector<array<Scalar,3> > >& training_crds, //!< locations of blob-like things
                                const vector<vector<bool> >& training_accepted, //!< classify each blob as "accepted" (true) or "rejected" (false)
                                Scalar *pthreshold_lower_bound = nullptr, //!< return threshold to the caller
                                Scalar *pthreshold_upper_bound = nullptr, //!< return threshold to the caller
                                SortCriteria sort_blob_criteria = SORT_DECREASING_MAGNITUDE, //!< give priority to high or low scoring blobs?
                                ostream *pReportProgress = nullptr //!< report progress back to the user?
                                )
{
  int Nsets = training_crds.size();
  assert(Nsets == training_accepted.size());
  assert(Nsets == blob_crds.size());
  assert(Nsets == blob_diameters.size());
  assert(Nsets == blob_scores.size());

  vector<array<Scalar,3> > final_training_crds;
  vector<bool> final_training_accepted;
  vector<Scalar> final_training_scores;

  // Loop over all of the different training sets:
  for (int I = 0; I < Nsets; I++) {

    assert(blob_crds[I].size() == blob_diameters[I].size());
    assert(blob_crds[I].size() == blob_scores[I].size());

    vector<array<Scalar,3> > training_crds_I;
    vector<bool> training_accepted_I;
    vector<Scalar> training_scores_I;
    FindBlobScores(training_crds[I],
                   training_accepted[I],
                   training_crds_I,
                   training_accepted_I,
                   training_scores_I,
                   blob_crds[I],
                   blob_diameters[I],
                   blob_scores[I],
                   sort_blob_criteria,
                   pReportProgress);

    // Concatinate all of the training data together.
    final_training_crds.insert(final_training_crds.end(),
                               training_crds_I.begin(),
                               training_crds_I.end());
    final_training_accepted.insert(final_training_accepted.end(),
                                   training_accepted_I.begin(),
                                   training_accepted_I.end());
    final_training_scores.insert(final_training_scores.end(),
                                 training_scores_I.begin(),
                                 training_scores_I.end());

  } //for (int I = 0; I < Nsets; I++)

  size_t N = final_training_crds.size();
  assert(N == final_training_scores.size());
  assert(N == final_training_accepted.size());

  // make sure that both positive and negative training data sets are non-empty
  size_t Nn = 0;
  size_t Np = 0;
  for (size_t i=0; i < N; i++) {
    if (final_training_accepted[i])
      Np++;
    else
      Nn++;
  }
  _ComplainIfTrainingDataEmpty(Nn, Np);
  
  if (pReportProgress)
    *pReportProgress
      << "  examining training data to determine optimal thresholds\n";

  _ChooseThresholdInterval(final_training_scores,
                           final_training_accepted,
                           pthreshold_lower_bound,
                           pthreshold_upper_bound,
                           pReportProgress);
                             
} //_ChooseBlobScoreThresholdsMulti()





} //namespace visfd



#endif //#ifndef _FEATURE_IMPLEMENTATION_HPP
