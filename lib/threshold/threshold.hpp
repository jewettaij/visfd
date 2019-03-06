#ifndef _THRESHOLD_HPP
#define _THRESHOLD_HPP

#include <cmath>
#include <cassert>
using namespace std;


template<class X>
bool IsBetween(X x, X a, X b) {
  return ((a <= x) && (x < b)) || ((b < x) && (x <= a));
}


/// @brief  Read an intensity value and return a new intensity value
///         in the range from "outA" or "outB".
///
///         If thresh_01_a < thresh_01_b
///
/// @code
///
///           output
///           intensity
///    outB    /|\                              __________________\
/// (usually 1) |                           _.-'                  /
///             |                       _,-'                 
///             |                   _,-'            
///    outA     |________________,-'                     _________\ input
/// (usually 0)               thresh         thresh               / intensity
///                            01_a           01_b
/// @endcode
///
///
///     OR    if thresh_01_b < thresh_01_a
///
/// @code
///
///          output
///          intensity
///
///            /|\
///     outB    |________________
/// (usually 1) |                `-._
///             |                    `-._ 
///             |                        '-._            
///     outA    |                            `-.___________________\ input
/// (usually 0)               thresh         thresh                / intensity
///                            01_b           01_a 
///
/// @endcode

template<class Number>
Number Threshold2(Number intensity,
                  Number threshold_01_a,
                  Number threshold_01_b,
                  Number outA = 0,
                  Number outB = 1)
{

  Number g; // "g" is the new intensity after threshold filter.

  if (IsBetween(intensity,
                threshold_01_a,
                threshold_01_b))
    g =
      (intensity-threshold_01_a) /
      (threshold_01_b-threshold_01_a);
  else if ((intensity-threshold_01_a)
           *
           (threshold_01_b-threshold_01_a)
           > 0.0)
    g = 1.0;
  else
    g = 0.0;

  return outA + g*(outB-outA);
} // Threshold2()




/// @brief  Read an intensity value and return a new intensity value
///         in the range from "outA" or "outB".
///
/// If the caller selects thresh_10_b < thresh_01_a, then
/// the thresholding function is:
///
/// @code
///
/// output
/// intensity
///
///     /|\
/// outB |                 ________________                
///      |             _,-'                `-._
///      |         _,-'                        `-._
/// outA |______,-'                                `-._______\ input
///           thresh    thresh          thresh    thresh     / intensity
///            01_a      01_b            10_a      10_b
///
/// @endcode
///
///  Or, if the caller selects thresh_10_b < thresh_01_a, then use:
///
/// @code
///
///     /|\
/// outB |_____                                       _______\ input
///      |     `-._                               _.-'       / intensity
///      |         `-._                       _,-'
/// outA |             `-._________________,-'
///           thresh    thresh          thresh    thresh
///            10_b      10_a            01_b      01_a
///
/// @endcode

template<class Number>
Number Threshold4(Number intensity,
                  Number threshold_01_a,
                  Number threshold_01_b,
                  Number threshold_10_a,
                  Number threshold_10_b,
                  Number outA = 0,
                  Number outB = 1)
{
  Number g;
  g = Threshold2(intensity,
                 threshold_01_a,
                 threshold_01_b);

  if ((threshold_01_b == threshold_10_a) &&
      (threshold_01_b == threshold_10_b))
    return g; // Default: if 2 extra arguments don't make sense, ignore them

  if (IsBetween(intensity,
                threshold_01_a,
                threshold_01_b))
    g = Threshold2(intensity,
                   threshold_01_a,
                   threshold_01_b);
  else if (IsBetween(intensity,
                     threshold_10_a,
                     threshold_10_b))
    g = Threshold2(intensity,
                   threshold_10_a,
                   threshold_10_b);
  else if (threshold_01_b <= threshold_10_a) {
    if (IsBetween(intensity,
                  threshold_01_b,
                  threshold_10_a))
      g = 1.0;
    else
      g = 0.0;
  }
  else if (threshold_10_b <= threshold_01_a) {
    if (IsBetween(intensity,
                  threshold_10_b,
                  threshold_01_a))
      g = 0.0;
    else
      g = 1.0;
  }
  else
    assert(false); //then thresholds are not in increasing or decreasing order


  return outA + g*(outB-outA);

} // Threshold4()



/// @brief  Read an intensity value and return a new intensity value
///         which is either "outA" or "outB".
///
/// If the caller selects range_a < range_b, then
/// the thresholding function is:
///
/// @code
///
/// output
/// intensity
///
///     /|\
/// outA |              ___________________
///      |             |                   |
///      |             |                   |
/// outB |_____________|                   |________________\ input
///                  range_a             range_b            / intensity
///
/// @endcode
///
///  Or, if the caller selects range_b < range_a, then:
///
/// @code
///
///     /|\
/// outA |_____________                     ________________\
///      |             |                   |                /
///      |             |                   |
/// outB |             |___________________|                
///                  range_b             range_a
///
/// @endcode

template<class Number>
Number SelectIntensityRange(Number intensity,
                            Number range_a,
                            Number range_b,
                            Number outA = 0,
                            Number outB = 1)
{
  Number g;

  if (range_a < range_b) {
    if (IsBetween(intensity, range_a, range_b))
      g = 1.0;
    else
      g = 0.0;
  }
  else {
    if (IsBetween(intensity, range_b, range_a))
      g = 0.0;
    else
      g = 1.0;
  }

  return g;
} // SelectIntensityRange



/// @brief  Read an intensity value and return a new intensity value
///         in the range from "outA" or "outB".
///         In this case, the mapping function is an unnormalized Gaussian
///         of width "sigma", centered at x0, with base outA
///         and height outA-outB.
/// @code
///
/// outB /|\                      _---_
///       |                     ,'  :  `.
///       |                    /    :    \
///       |                _.-'     :sigma`-,_
/// outA  |------------''''         :         ````------------
///                                x0
/// @endcode

template<class Number>
Number SelectIntensityRangeGauss(Number intensity,
                                 Number x0,
                                 Number sigma,
                                 Number outA = 0,
                                 Number outB = 1)
{
  Number dx = intensity - x0;
  Number xr = (dx/sigma);
  return outA + (outB-outA) * exp(-0.5*xr*xr);
} // SelectIntensityRangeGauss


#endif //#ifndef _THRESHOLD_HPP
