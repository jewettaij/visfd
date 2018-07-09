mrc_simple_filter
===========

This is a collection of programs (mainly "filter_mrc" and "combine_mrc") for 3-D image processing using the MRC file format.

## programs included with this repository:

After compilation, all programs will be located in the "*bin/*" subdirectory.  Here is a brief description of some of them:


## filter_mrc
![example: a slice through a tomogram with a visible nucleoid](./doc/images/nucleoid_example_Hylemonella_gracilis.jpg)
![example: red: scale-free-blob-detection ("-blobr"), blue: fluctuation-filter ("-fluct")](./doc/images/nucleoid_example_Hylemonella_gracilis__red_blob_detection__blue_fluctuation_filter.jpg)

**filter_mrc** detects features and applies simple filters to 3-D images.
It supports 
[3D scale-free blob-detection](https://en.wikipedia.org/wiki/Blob_detection),
and local minima/maxima finding (with non-max suppresion).
Filters include
thresholding,
low-pass, high-pass, 
[generalized](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1)
[Gaussian](https://en.wikipedia.org/wiki/Gaussian_blur),
and
[DOG](https://en.wikipedia.org/wiki/Difference_of_Gaussians)
filters, and others.
A image *mask* can be used to exclude certain
voxels or regions from consideration.
(Typically these are voxels which have been characterized previously.)
A list of detected objects can be saved to a text file.
Annotated images can be saved to a new MRC/REC file.
Documentation for this program is located
[here](./doc/doc_filter_mrc.md).


## combine_mrc
**combine_mrc** is a program for combining two volumetric images (i.e. tomograms, both of identical size) into one image/tomogram, using a combination of addition, multiplication, and thresholding operations.  These features can be used perform binary operations between two images (which are similar to "**and**", "**or**", and "**not**" operations.)
Documentation for this program is located
[here](./doc/doc_combine_mrc.md).

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


## histogram_mrc.py
**histogram_mrc.py** is a graphical python program which displays the
histogram of voxel intensities contained in an MRC file.
It can be useful when deciding what thresholds to use
with in the "**filter_mrc**" and "**combine_mrc**" programs.
This software requires the *matplotlib* and *mrcfile* python modules
(both of which can be installed using pip).
Documentation for this program is located
[here](./doc/doc_histogram_mrc.md).


## Development Status: *alpha*
Program names and command line
arguments may change in the future.


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
