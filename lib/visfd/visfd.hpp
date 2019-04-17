///   @file visfd.hpp
///   @brief a meta header file including all of the header files in the
///          VISFD library
///   @author Andrew Jewett
///   @date 2019-4-15

#ifndef _VISFD_HPP
#define _VISFD_HPP

// higher level functions and classes are defined here:

#include <feature.hpp>        // simple feature detectors (blob, surface,...)
#include <segmentation.hpp>   // common segmentation and cluster operations
#include <draw.hpp>           // functions for drawing and annotation
#include <filter3d.hpp>       // common 3D filters


// lower level functions and classes (upon which files above depend):


#include <err_visfd.hpp>      // defines the "VisfdErr" exception type
#include <eigen3_simple.hpp>  // defines matrix diagonalizer (DiagonalizeSym3())
#include <lin3_utils.hpp>     // defines DotProduct3(),CrossProduct(),quaternion
#include <visfd_utils.hpp>    // defines invert_permutation(), AveArray(), ...
#include <alloc2d.hpp>        // defines Alloc2D() and Deallox2D()
#include <alloc3d.hpp>        // defines Alloc3D() and Dealloc3D()
#include <filter1d.hpp>       // defines "Filter1D" (used in ApplySeparable())
#include <filter2d.hpp>       // defines "Filter2D"
#include <multichannel_image3d.hpp> // defines "CompactMultiChannelImage3D"


#endif //#ifndef _VISFD_HPP
