[![CircleCI](https://circleci.com/gh/jewettaij/visfd.svg?style=svg)](https://circleci.com/gh/jewettaij/visfd)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jewettaij/visfd)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()


visfd
===========

## Volumetric Image toolkit for Simple Feature Detection

This is a small C++ template library for general 3D image processing
("[visfd.hpp](./lib/visfd/visfd.hpp)"),
as well as reading & writing image files using the the MRC file format
("[mrc_simple.hpp](./lib/mrc_simple/mrc_simple.hpp)").
It is also a collection of stand-alone programs
which use this library
(including "[filter_mrc](./doc/doc_filter_mrc.md)",
 "[combine_mrc](./doc/doc_combine_mrc.md)",
 "[voxelize_mesh](./doc/doc_voxelize_mesh.md)",
 "[sum_voxels](./doc/doc_sum_voxels.md)",
 and
 "[pval_mrc](./doc/doc_pval_mrc.md)").
They are documented [here](./doc).
Multiprocessor support is implemented using
[OpenMP.](https://en.wikipedia.org/wiki/OpenMP)


## *Alternatives to VISFD*
Much more comprehensive libraries and software tools are available
for 3-D image processing, such as [scikit-image](https://scikit-image.org) and
[scipy.ndimage](https://docs.scipy.org/doc/scipy/reference/ndimage.html).
*(MRC files can be read into python arrays using the
[mrcfile](https://mrcfile.readthedocs.io/en/latest/readme.html#basic-usage)
module.)*
VISFD compliments these libraries by providing
more robust curve and surface detectors
(following [this paper](http://dx.doi.org/10.1016/j.jsb.2014.02.015)),
as well as tools for clustering and closed surface segmentation.


## programs included with this repository:

After compilation, all programs will be located in the "*bin/*" subdirectory.  Here is a brief description of some of them:


## filter_mrc
![example: a slice through a tomogram with a visible nucleoid](./doc/images/nucleoid_example_Hylemonella_gracilis.jpg)
![example: red: scale-free-blob-detection ("-blobr"), blue: fluctuation-filter ("-fluct")](./doc/images/nucleoid_example_Hylemonella_gracilis__red_blob_detection__blue_fluctuation_filter.jpg)
![example: membrane reconstruction using tensor voting and PoissonRecon (after cleaning up with meshlab)](./doc/images/nucleoid_example_Hylemonella_gracilis_inner_membrane.jpg)

**filter_mrc** is a stand-alone program which uses many of the
features of the **visfd** library.
This program was intended to be used for automatic
[membrane (surface) detection](https://www.ncbi.nlm.nih.gov/pubmed/24625523),
[surface closure](https://stackoverflow.com/questions/51149213/how-to-avoid-hole-filling-in-surface-reconstruction),
[edge detection](./doc/doc_filter_mrc.md#-edge-thickness), 
[filament (curve) detection](./doc/doc_filter_mrc.md#Detecting-curves), and
[scale-free blob-detection](https://en.wikipedia.org/wiki/Blob_detection).
Images can be segmented into distinct contiguous objects,
using a [variety](https://imagej.net/plugins/classic-watershed#introduction)
of [strategies](./doc/doc_filter_mrc.md#-connect-threshold),
and closed compartments can be hierarchically segmented using
[*voxelize_mesh.py*](doc/doc_voxelize_mesh.md).
A list of detected objects can be sorted, clustered
and saved to standard files for further analysis by 3rd party programs,
*(including [*PoissonRecon*](https://github.com/mkazhdan/PoissonRecon),
[*meshlab*](http://www.meshlab.net), and
[*voxelize_mesh.py*](doc/doc_voxelize_mesh.md)
to generate smooth, connected, closed surface meshes.)*

This program includes a manual (text-mode) 3D image editor,
as well as a variety of filters to clean up 3D images, including
low-pass, high-pass, dilation, erosion, opening ,closing,
thresholding, clipping, brightness inversions, fluctuation detectors,
[generalized Gaussian blur](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1),
[DoG](https://en.wikipedia.org/wiki/Difference_of_Gaussians),
[LoG](https://en.wikipedia.org/wiki/Blob_detection#The_Laplacian_of_Gaussian),
[Ridge-detection](https://en.wikipedia.org/wiki/Ridge_detection),
and
[3D tensor voting](https://www.ncbi.nlm.nih.gov/pubmed/24625523)
filters.
*(As of 2021-7-12, this program does not have a graphical user interface, so
examples explain how to use this software with
[IMOD](https://bio3d.colorado.edu/imod/) and
[meshlab](http://www.meshlab.net)
to display results and choose parameters.  Hopefully these tools can eventually
be integrated with visualizers like IMOD/3dmod, chimerax, or tomviz.)*


Documentation for this program is located
[**here**](./doc/doc_filter_mrc.md).
The source code for the filters used by this program
is located
[**here**](./lib/visfd/).
This program currently only supports the .mrc/.rec image file format.


## voxelize_mesh.py
**voxelize_mesh.py** is a program that finds the voxels in a volumetric
image that lie within the interior of a closed surface mesh.
It was intended for segmenting the interiors of membrane-bound
compartments in tomograms of cells.
The mesh files that this program reads are *typically* generated by
*filter_mrc* ([together with other tools](./doc/doc_filter_mrc.md#example-3)).
However it can read any standard PLY file containing a closed polyhedral mesh.
This program currently only supports the .mrc/.rec image file format.
Documentation for this program is located
[here](./doc/doc_voxelize_mesh.md).
**WARNING: This experimental program is very slow and currently requires
a very large amount of RAM.**


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


### Development Status: *alpha*

Program names, command line arguments, file names, and function names
(in the API) may change in the future.
Automated testing was added,
however as of 2020-12-15, some commits still (temporarily) break everything.
(...because I'm too lazy to use branch & merge.
 This usually gets fixed within 24 hours.
 If the build is failing, choose a previous commit.)

### Development Timeline: 2021

Work on this project was temporarily halted on 2019-7,
however I occasionally make small feature updates as I need them.
We hope to finish the remaining features
and submit a paper on this software in 2021.
Nevertheless, if you find a bug, please report it.  I will fix it.




## Requirements

- **16GB** of RAM or higher.
(64GB is recommended.  For membrane detection,
 your RAM must exceed 11x the size of the tomogram that you are analyzing.
 The *voxelize_mesh.py* program requires even more memory.
 You can reduce the memory needed and computation time dramatically
 by cropping or binning your tomogram.)
- A terminal (running BASH or CSH/TCSH) where you can enter commands.
  (*filter_mrc* is currently a text-only program.)
- A [C++ compiler](#Supported-compilers)
- [make](https://en.wikipedia.org/wiki/Make_(software))
- Software to visualize MRC/REC/MAP files
(such as [IMOD/3dmod](https://bio3d.colorado.edu/imod/))


### Recommended:

- A text editor.  (Word, Wordpad, and Notepad will not work.)
  Popular graphical text editors
  include **Atom**, **Sublime**, **Notepad++**, and **VSCode**.
  Older, non-graphical programs include **vim**, **emacs**,
  **nano**, **ne**, and **jove**.
  (Apple's TextEdit can be used if you save the file as *plain text*.)
- python
- The "numpy", "matplotlib", "mrcfile", and "pyvista" python modules
  (These are are installable via "pip3" or "pip".)
- [**PoissonRecon**](https://github.com/mkazhdan/PoissonRecon).
- Software to visualize mesh files (such as [meshlab](http://www.meshlab.net)).
- A computer with *at least* 4 CPU cores (8 threads).
- A **hard drive**.  (An HDD, not an SSD.)  The "filter_mrc"
  program creates many large temporary files.  Although they
  can be deleted, I worry that over time, this frequent writing and erasing
  could wear down an SSD.  Old-fashioned magnetic hard drives are
  better able tolerate this kind of usage.



### Supported compilers

This code has been tested with the **CLANG** and **GCC** (v9.3.0)
compilers.  I have run into
[difficulties in the past](https://github.com/jewettaij/visfd/issues/2).
with older versions of GCC (7.5.0).


## Compilation

## Linux:
```
    source setup_gcc.sh
    # (...or "source setup_gcc.sh", if you prefer gcc)
    make
```
(If you are not using the bash shell, enter "bash" into the terminal beforehand.)

Then copy the newly-compiled programs (and .PY files) from the
various bin/ subdirectories to a directory that is in your
[PATH](http://www.linfo.org/path_env_var.html)
(such as ~/bin/ or /usr/local/bin/).
*(This includes "filter_mrc", "combine_mrc", "voxelize_mesh.py",
and "sum_voxels".)*
Then log out and log in again for the changes to take effect.


#### Warning: Avoid the **-Ofast** compiler argument
As of 2020-2-11, we have had better results using
[**-O3 -ffast-math**](https://github.com/jewettaij/visfd/issues/6) instead.
(The problem with -Ofast seems to only appear during membrane detection.
 See [setup_gcc.sh](setup_gcc.sh) for a list of suggested compiler flags.)

*Note: This program has been carefully checked for memory errors using valgrind
and was found to be clean.  The issues mentioned above only occur when compiler
optimization flags like -Ofast are in use.*


## Windows:

It is recommended that you install the BASH shell environment on your computer,
along with *make* and either *gcc* or *clang*.  Once you have done that,
you can follow the instructions above for linux users.
There are several ways to to create a BASH environment,
but perhaps the easiest method is to install
[Windows Subsystem for Linux (WSL2)](https://docs.microsoft.com/en-us/windows/wsl/install-win10),
***or***
[virtualbox](https://www.virtualbox.org)
(In the later case, you will also need to install a linux distribution,
preferably with a lightweight
desktop such as [xubuntu](https://xubuntu.org).)
Alternatively, you can try 
[Hyper-V](https://www.nakivo.com/blog/run-linux-hyper-v/)
or (if you have an older version of windows)
[CYGWIN](https://www.cygwin.com/).

WSL and virtualbox are virtual machines that allow you to run an
alternate operating system from within windows.
In this case that operating system is linux.  The BASH shell and the
compiler tools that you need can be easily installed from within in linux.
Both WSL and virtualbox also create an alternate filesystem inside windows
where the linux operating system is stored.  Software (like *visfd/filter_mrc*)
that you download and install there can access the files in that filesystem.
So you may need to copy your tomograms and otherf files to this fileystem
beforehand.  If you **are using WSL or WSL2**, then you should
[use caution when using windows programs to edit files stored there](https://devblogs.microsoft.com/commandline/do-not-change-linux-files-using-windows-apps-and-tools/).
This includes text editors, image editors, and tomography processing software
(such as IMOD).
One possible way to avoid problems is to try to restrict yourself to using
programs which you downloaded and installed directly from within the
WSL or WSL2 environment *(if possible)*.  In particular, a couple of the
features of the *filter_mrc* program require you to learn and use a
(unix-style) text editor.  (Word, Wordpad, and Notepad will not work.)
Again, it is a good idea to install and run such programs from within WSL,
not windows.



## Apple Mac:

You must install OpenMP, and a compiler that supports it.
Otherwise, this software will run at an intolerably slow speed.
Last I checked OpenMP is disabled in MacOS by default.
You must figure out how to install it.
*(If it helps, there are some old instructions for doing that
[here](https://iscinumpy.gitlab.io/post/omp-on-high-sierra/).)*

I do not have access to a computer running MacOS.
If you figure out how to install OpenMP, please contact me
(or create a pull-request) so that we can update this section.


Once OpenMP is installed, compile the visfd tools using:
```
    source setup_gcc.sh
    # (...or "source setup_clang.sh", if you prefer clang)
    make
```

Then copy the newly-compiled programs (and .PY files) from the
various bin/ subdirectories to a directory that is in your
[PATH](http://www.linfo.org/path_env_var.html)
(such as ~/bin/ or /usr/local/bin/).
*(This includes "filter_mrc", "combine_mrc", "voxelize_mesh.py",
and "sum_voxels".)*


## License

All of the code in this repository
(except for code located in "lib/mrc_simple" and "lib/visfd/eigen_simple.hpp")
is available under the terms of the terms of the
MIT license.  (See "[LICENSE.md](./LICENSE.md)")


### Additional license dependencies

#### MPL-2.0 licensed code (eigen_simple)
The "lib/visfd/eigen3_simple.hpp" file contains code from
[Eigen](http://eigen.tuxfamily.org) which requires the 
[MPL-2.0 license](lib/eigen_simple/COPYING.MPL2).


#### GPLv2 licensed code (mrc_simple)
A small subset of the code in "lib/mrc_simple" was adapted from ***IMOD***.
The IMOD code uses the GPL license (version 2), which is more restrictive.
License details for the "mrc_simple" library can be found in the
[COPYRIGHT.txt](lib/mrc_simple/COPYRIGHT.txt)
file located in that directory.
*If you write your own code using the "visfd" library to analyze 3D images
(which you have loaded into memory by some other means),
then you can ignore this notice.*


## Funding

VISFD is funded by NIH grant R01GM120604.
