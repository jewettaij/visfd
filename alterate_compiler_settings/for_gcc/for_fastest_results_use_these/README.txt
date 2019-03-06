On 2019-3-06, I stopped using the "-O3" flag optimization flag with the gcc
compiler.  I noticed that using these flags caused strange behavior in
several different places in the code.
(Most recently, ConnectedClusters() and possibly also ApplySeparable3D().
 At the time, both functions were located in the "filter3d.hpp" file.)

This directory contains the original compiler settings used with
GCC at that time.  If I can figure out why these functions are in conflict
with GCC's -O3 flag, then I will move these files to a more visible location.
It would be nice to have code which is compiler agnostic.

-andrew 2019-3-06
