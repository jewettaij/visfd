[![Build Status](https://travis-ci.org/jewettaij/visfd.svg?branch=master)](https://travis-ci.org/jewettaij/visfd.svg?branch=master)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

visfd
===========

## Volumetric Image toolkit for Simple Feature Detection

This is a C++ template library for 3D image processing
("[visfd.hpp](./lib/visfd/)")
as well as reading & writing image files using the the MRC file format
("[mrc_simple.hpp](./lib/mrc_simple/mrc_simple.hpp)").
It is also a collection of stand-alone programs
which use this library
(including "[filter_mrc](./doc/doc_filter_mrc.md)",
 "[combine_mrc](./doc/doc_combine_mrc.md)",
 and
 "[pval_mrc](./doc/doc_pval_mrc.md)").
They are documented [here](./doc).
Multiprocessor support is implemented using
[OpenMP.](https://en.wikipedia.org/wiki/OpenMP)


#### WARNING: This project is changing.  Optimizations and features are being added rapidly.  Command-line syntax and the API could change.  Some features don't work well yet.  (-andrew 2019-4-11)


## programs included with this repository:

After compilation, all programs will be located in the "*bin/*" subdirectory.  Here is a brief description of some of them:


## filter_mrc
![example: a slice through a tomogram with a visible nucleoid](./doc/images/nucleoid_example_Hylemonella_gracilis.jpg)
![example: red: scale-free-blob-detection ("-blobr"), blue: fluctuation-filter ("-fluct")](./doc/images/nucleoid_example_Hylemonella_gracilis__red_blob_detection__blue_fluctuation_filter.jpg)  *<-- 2D slice through a masked, segmented 3D tomogram:*

**filter_mrc** is a stand-alone program which uses many of the
features of the **visfd** library.
It was intended to be used for
[membrane (surface) detection](https://www.ncbi.nlm.nih.gov/pubmed/24625523),
[surface closure](https://stackoverflow.com/questions/51149213/how-to-avoid-hole-filling-in-surface-reconstruction),
filament (curve) detection (*available soon*),
[scale-free blob-detection](https://en.wikipedia.org/wiki/Blob_detection),
and
[watershed segmentation](https://imagej.net/Classic_Watershed)
with non-max suppression.
A list of detected objects can be sorted, clustered, and saved to a text file.

*(filter_mrc can also be used to apply simple filters to images, including
low-pass, high-pass,
thresholding,
brightness inversions,
fluctuations,
[generalized Gaussian blur](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1),
[DoG](https://en.wikipedia.org/wiki/Difference_of_Gaussians),
[LoG](https://en.wikipedia.org/wiki/Blob_detection#The_Laplacian_of_Gaussian),
[Ridge-detection](https://en.wikipedia.org/wiki/Ridge_detection),
and
[3D tensor voting](https://www.cs.stevens.edu/~mordohai/public/TensorVotingTutorial_2007.pdf)
filters.)*

Documentation for this program is located
[here](./doc/doc_filter_mrc.md).
The source code for the filters used by this program
is located
[here](./lib/visfd/).
This program currently only supports the .mrc/.rec image file format.


## combine_mrc
**combine_mrc** is a program for combining two volumetric images (i.e. tomograms, both of identical size) into one image/tomogram, using a combination of addition, subtraction, multiplication, division, and thresholding operations.  These features can be used perform binary operations between two images (which are similar to "**and**", "**or**", and "**not**" operations.)
Documentation for this program is located
[here](./doc/doc_combine_mrc.md).

## histogram_mrc.py
**histogram_mrc.py** is a graphical python program which displays the
histogram of voxel intensities contained in an MRC file.
It can be useful when deciding what thresholds to use
with in the "**filter_mrc**" and "**combine_mrc**" programs.
Voxels and regions in the image can be excluded from consideration
by using the "-mask" and "-mask-select" arguments.
This software requires the *matplotlib* and *mrcfile* python modules
(both of which can be installed using pip).
Documentation for this program is located
[here](./doc/doc_histogram_mrc.md).

## sum_voxels
**sum_voxels** is a program for estimating volumes.
It is a simple program which
reads an MRC (.REC) file as an argument
and computes the sum of all the voxel intensities.
(Typically the voxel intensities are either 1 or 0.
 The resulting sums can be converted into volumes
 either by multiplying by the volume-per-voxel,
 or by specifying the voxel width using the "-w" argument,
 and including the "-volume" argument.)
For convenience, threshold operation can be applied
(using the "-thresh", "-thresh2", and "-thresh4" arguments)
so that the voxels intensities vary between 0 and 1
before the sum is calculated.
The sum can be restricted to certain regions
(by using the "-mask" and "-mask-select" arguments).
Documentation for this program is located
[here](./doc/doc_sum_voxels.md).


## pval_mrc
**pval_mrc** is a program for estimating the probability
that a cloud of points in an image is distributed randomly.
It looks for regions of high (or low) density in an image.
(The user can specify a *-mask* argument to perform the analysis
 in small, confined, irregularly-shaped subvolumes from a larger image.)
Documentation for this program is located
[here](./doc/doc_pval_mrc.md).


## Development Status: *pre-alpha*
Program names, command line arguments, file names, and function names
(in the API) may all change in the future.
Automated testing was added,
however as of 2019-2-27, some commits still (temporarilly) break everything.
(...because I'm too lazy to use branch & merge.
 This usually gets fixed within 24 hours.
 If the build is failing, choose a previous commit.)


## Compilation

## Linux:

    source setup_clang.sh
    make

(If you are not using the bash shell, enter "bash" into the terminal beforehand.)


Note: As of 2019-3-06, we have had good results using the **CLANG** compiler.
      ([Avoid using GCC](https://github.com/jewettaij/visfd/issues/2) for now.)


## Windows:

   It is recommended that you install the BASH shell environment on your computer, along with *clang* and *make*.  (If you decide to use a different compiler, modify the "setup_clang.sh" file accordingly.)  There are several ways to do that.

   Perhaps the easiest way is to install [virtualbox](https://www.virtualbox.org) in windows together with a linux distribution with a lightweight desktop, such as [xubuntu](https://xubuntu.org).  Alternatively, if you are using Windows 10 or later, you can try installing the "Windows Subsystem for Linux (WSL)", as explained
[here](https://solarianprogrammer.com/2017/04/15/install-wsl-windows-subsystem-for-linux/)
and
[here](https://msdn.microsoft.com/en-us/commandline/wsl/faq),
or
[Hyper-V](https://blogs.windows.com/buildingapps/2018/09/17/run-ubuntu-virtual-machines-made-even-easier-with-hyper-v-quick-create/).
Otherwise, if you are using an older version of windows, try installing
[CYGWIN](https://www.cygwin.com/) instead.

## Apple Mac:

*WARNING: The following proceedure has not been tested.
 If you have a mac, please test this and let me know if it worked
 (or if something else worked).
 -andrew 2019-3-08*

First follow the instructions
[here](https://iscinumpy.gitlab.io/post/omp-on-high-sierra/)
to install OpenMP support:

Then compile visfd using
```
    source setup_clang.sh
    make
```

If this doesn't work and you are desperate, you can try this instead:
```
    source alternate_compiler_settings/for_debugging_and_profiling/setup_gcc_dbg.sh
    make
```
*Unfortunately the resulting (already slow) program will run at least 20x slower.*


## Requirements:

**8GB** or higher is required.
(For membrane detection,
 your RAM must exceed 11x the size of the tomogram that you are analyzing.)

The CLANG compiler is strongly recommended.
([The clustering feature does not yet work with GCC.](https://github.com/jewettaij/visfd/issues/2))

Some programs
(such as "histogram_mrc.py" and "draw_filter_1D.py")
require python, along with the
"numpy", "matplotlib", and "mrcfile" modules.
(installable via "pip")

*Automatic membrane surface closure* (a.k.a. surface reconstruction) currently
requires "**PoissonRecon**" (an external program) which you can download
[here](https://github.com/mkazhdan/PoissonRecon).

## License

These programs are available under the terms of the open-source 3-clause BSD
license.  (See "[LICENSE.md](./LICENSE.md)")
