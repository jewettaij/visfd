mrc_simple_filter
===========

This is a collection of programs (mainly "filter_mrc" and "combine_mrc") for 3-D image processing using the MRC file format.

## programs included with this repository:

After compilation, all programs will be located in the "*bin/*" subdirectory.  Here is a brief description of each:


## filter_mrc

![A comparison of Difference-of-Gaussians Difference-of-Generalized-Gaussians filter weights](./doc/images/example_dogxy_w=2.516nm_a=13nm_b=20nm_m=2_n=2__vs__m=6_n=16.svg)

**filter_mrc** applies a filter to a tomogram in the X,Y,Z directions
and saves the result as a new .mrc/.rec file.
This program can be used to rescale or invert a 3-D image, remove its high or low spatial frequencies,
(smoothing, edge detection, band-pass filtering),
and perform 
[Scale-Free Blob-Detection](https://en.wikipedia.org/wiki/Blob_detection).
Currently, the program supports the following filters:
([generalized](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1))
[Gaussians](https://en.wikipedia.org/wiki/Gaussian_blur),
([generalized](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1))
[Difference-of-Gaussians](https://en.wikipedia.org/wiki/Difference_of_Gaussians),
inverse ("-invert", bright <--> dark),
and
several threshold filters ("-thresh", "-thresh2", "-thresh4".
(Threshold filters are useful for selecting voxels with intensities
 which lie within a certain range.)
Both isotropic and anisotropic versions of these filters are supported.
(Fast [separable](https://en.wikipedia.org/wiki/Separable_filter) filters are used when possible.)
This program can also perform
[template-matching](https://en.wikipedia.org/wiki/Template_matching)
with *error-estimation* 
(but only for simple spherical template shapes).
The user can exclude certain voxels or regions from consideration
(using the "-mask" argument).
Documentation for this program is located
[here](./doc/doc_filter_mrc.md).


## combine_mrc
**combine_mrc** is a program for combining two volumetric images (i.e. tomograms, both of identical size) into one image/tomogram, using a combination of addition, multiplication, and thresholding operations.  These features can be used perform binary operations between two images (which are similar to "**and**" and "**or**" operations.  "**not**" operations are also possible.)  As with the "*filter_mrc*" program, you can also use the "-mask" argument to restrict the operation to certain voxels from the image.  (See the documentation for that tool for details.)
Documentation for this program is located
[here](./doc/doc_combine_mrc.md).


## histogram_mrc.py
**histogram_mrc.py** is a graphical python program which displays the
histogram of voxel intensities contained in an MRC file.
It can be useful when deciding what thresholds to use
with in the "**filter_mrc**" and "**combine_mrc**" programs.
This software requires the *matplotlib* and *mrcfile* python modules
(both of which can be installed using pip).
Documentation for this program is located
[here](./doc/doc_histogram_mrc.md).

## sum_voxels
**sum_voxels** is a simple program which
reads an MRC (.REC) file as an argument
and computes the sum of all the voxel intensities.
(Typically the voxel intensities are either 1 or 0.
 The resulting sums can be converted into volumes by
 multiplying by the volume-per-voxel.)
For convenience, the voxel intensities can be clipped or scaled
between 0 and 1 beforehand
(by using the "-thresh", "-thresh2", and "-thresh4" arguments),
and the sum can be restricted to certain regions
(by using the "-mask" and "-mask-select" arguments). 


## Development Status: *pre-alpha*
As of 2018-2-12, this code has been lightly tested.
Program names and command line
arguments may change in the future.
Tagged *releases* of this repository should have some semblance
of functionality.  (More recent commits may not even compile.)


## Compilation

## Linux:

    cd src
    source setup_gcc_linux.sh
    make

(If you are not using the bash shell, enter "bash" into the terminal beforehand)

## Windows 10:

Install the Windows Subsystem for Linux (WSL) and run

    sudo apt-get install build-essential

and then follow the instructions above for linux.
(Older windows users can install Cygwin or MinGW.)

## Apple OSX:

    cd src
    source setup_gcc_OSX.sh
    make

## Requirements:

The optional "draw_filter_1D.py" script
(included in the "bin/filter_mrc" directory)
requires python, numpy, and matplotlib.
(It is useful only if you actually want to see
 the shape of the convolution filter that is currently in use.
 Most users can ignore this.)


## License

These programs are available under the terms of the open-source 3-clause BSD
license.  (See `LICENSE.md`.)
