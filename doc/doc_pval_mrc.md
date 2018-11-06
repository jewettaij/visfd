pval_mrc
===========

**pval_mrc** performs
[P-value](https://en.wikipedia.org/wiki/P-value)
analysis on the densities of
objects detected in 3D images.

# Description

The **pval_mrc** program does the following:

1. Generate an image which is a map of the the
   local density of particles at a particular resolution.
2. Find the minimum density (or maximum density) in that image.
   (Note: Local minima or maxima located at edge of the image/mask boundaries
          will be ignored.)
3. Estimate the probability that such an image would contain *any*
   regions whose density is that low (or that high).
   This number is printed to terminal.
4. If the user wants this to be done at multiple resolutions
   (ie. using the "*-scan-bin*" or "*-scan-gauss* arguments)
   then steps 1-3 are repeated for each resolution requested by the user.


### This program requires:

* a multi-column text file containing a
   **list of 3D coordinates** (typically in units of Angstroms)
   corresponding the the centers of objects that were detected in an image file.
   (This is specified using the "*-coords*" argument.  See example below.)
   The first three numbers on each line of this file should contain
   the coordinates of an object detected in the image.
   ```
   2016.3 4012.8 1881.6
   6201.6 7852.8 1113.6
   3168.9 5894.4 1996.8
   4108.8 2745.6 2131.2
     :      :      :
   ```
   These files are typically generated using "*filter_mrc -blob ...*".
   (See the docs for "[filter_mrc](./doc_filter_mrc.md)" for details.
   *Note:* These coordinates are in physical units by default.
   If you use the "*-w 1*" argument, these numbers will be in units of voxels
   and will be positive integers.)
* an **image file** (containing the image from which the objects were detected,
  specified using the "*-in*" argument). Alternatively, you can specify the size
  of the image (using the "*-image-size*" argument) without specifying a file.
* a **resolution, δ**.
  This is the resolution at which you want to calculate the density.
  (It is specified with the
   "*-bin-width*", "*-sigma*", "*-scan-bin*", or "*-scan-gauss*" arguments
   and it typically has units of Angstroms.  We denote it using "δ" below.)
   
* Optionally, you can also supply
  a **mask file** countaining another image of the *same size*.
  Any voxels in *this* image with brightness 0 correspond to 
  regions of space that you wish to *ignore* from the original image.
  (These voxels typically lie outside the region of interest.
   You can specify this file using the "*-mask*" argument.)



# Return values

The program will print 6 numbers to the standard out:
1.   The probability that a random distribution of particles in this region
     would *not* contain *any* "bins" (regions of space of volume δ^3)
     whose density exceeded the observed minimum density (at this resolution)
     **(This is the P-value.)**
2.   the minimum (or maximum) density at that resolution (bin-size),
3.   x (location coordinate),
4.   y (location coordinate),
5.   z (location coordinate) of the global minima (or maxima) of the density
     (typically in Angstroms).
6.   the resolution (bin-size) considered (typically in Angstroms)

If the user wants to repeat the analysis at multiple resolutions
(using the "*-scan-bin*" or "*-scan-gauss*" arguments),
then it will print all these numbers on separate lines,
one line per resolution.


### Examples:
```
pval_mrc -w 15.6 \
         -coords detected_blobs.txt \
         -in orig.rec \
         -out density_at_500angstroms.rec \
         -bin-width 500.0
```
The particles listed in the file ("detected_blobs.txt") are
(effectively) binned into boxes of width 500.0 Angstroms,
and the minimum density bin is calculated.
Finally the P-value for that density is printed to the standard out
(along with other numbers).
An image of the density at that resolution is also saved in a file
("density_at_500angstroms.rec")

In this example, the "*-w*" argument specifies
the width of each voxel (in this case 15.6 Angstroms).
(If instead you wish to measure distances in units of voxels, use "*-w 1*".)

The next example demonstrates how to perform the analysis over a range
of resolutions:

```
pval_mrc -w 19.2 \
         -coords detected_blobs.txt \
         -in orig.rec \
         -mask segmented_image_periplasm1_cytoplasm2.rec \
         -mask-select 2 \
         -maxima \
         -scan-bin-width 500.0 1000.0 1.1
```
In this example, the program
searches for density *maxima* instead of *minima*.
It calculates their probability at various resolutions using 
bin widths ranging from 500.0 to 1000.0 (in physical units, Angstroms),
increasing the bin width by a factor of 1.1 (10%) each time.


# Required Arguments:

#### -coords xyz.txt
  Specify a file containing at least 3 numbers on each line.
  The numbers are assumed to be non negattive and are in the same units
  specified by the bin-width and -w arguments (typically Angstroms).
  (To override this and specify the coordinates in units of voxels, 
  use the "-w 1" argument.  See below.)

#### -bin-width δ
  Specify the size of the bins you wish to use when constructing the histogram
  used to estimate the local density of points.  (See details below.)
  This is the resolution of the density.
  (The bin width, δ, is typically much wider than the width of a voxel.)
  The units of δ are typically in Angstroms
  (or whatever physical units specified by the -w argument or 
   specifed by the header of the MRC file).
  To override this and specify the bin width (δ) in voxels, 
  use the "-w 1" argument.  (See below.)


#### -in image_file.rec   *(or -image-size Nx Ny Nz)*

  Specify an MRC/.REC file whose voxels (after multiplication by the voxel
  width) encompass the points of interest.  A new image of the same size will
  be created which displays the density cloud.  The minimum and maximum
  density of this cloud will be used later to calculate the P-value.

#### -image-size Nx Ny Nz

  Alternatively you can specify the *size* of that image (in voxels) 
  using the "*-image-size Nx Ny Nz*" argument.
  (If you are also supplying a *mask*, then the *Nx,Ny,Nz* numbers should match
   the size of the *mask* image you provided.
   In normal usage, the contents of the image contained in the file specified
   by the "*-in*" argument is ignored.  The only reason an image
   is needed is to determine the size of the image which will be used
   for the calculation of the density.  Rather than loading an image
   this argument allows you to simply specify the size of that image in voxels.)
  *Users can specify either the "-in" or "-image-size" arguments, but not both.*



# Optional Arguments:

#### -w  voxel_width
  Specify the width of each voxel (in physical units, typically Angstoms).
  If specified, then all distance units will be expressed 
  using these units.
  (Densities will be expressed using the inverted cube of these units.)
  To express all distance (and density) units in voxels, use "-w 1".

#### -maxima 
  Search for maxima instead of minima, 
  and compute the probability of
  such maxima occuring (instead of minima).

#### -mask  mask_file.rec
  Specify an image file indicating which voxels should be ignored.
  This image should be the same size as the image specified by
  the "*-in*" or "*-image-size*" arguments.
  Voxels in this image with brightness 0 will be ignored.
  Typically these are voxels which lie outside the region of interest
  (typically a cell, organelle, or compartment subvolume).
  Points of interest (from the text) file which lie within
  ignored voxels will not be counted.

#### -mask-select  intensity_value

If the "*-mask-select*" argument is specified, then instead of considering all
voxels with non-zero values from the "mask" image,
only voxels whose mask intensity equals the number following
this argument will belong to the mask.
All other voxels will be ignored during filtering.
*(Note: This disables "soft" masking.  In other words,
  during the process of filtering, all selected voxels will be weighted equally
  during filtering, and all others will be completely ignored.)*

#### -out  density.rec
  Optionally, you can specify name of a **new image file**
  (in MRC/.REC format) whose voxel brightnesses correspond
  to the local density of particles at that location.
  This is the image from which the minimum and maximum densities are calculated.
  (Typically the density will be expressed in units of 1/Angstroms^3.
   Use "-w 1" if you wish to express the density in units of voxels.  In that
   case, the numbers in the coordinate file will need to be scaled accordingly.)

#### -gauss σ

  Specify the width (σ) of the Gaussian used directly 
  (instead of specifying the effective "bin" width, δ=σ√(2π))

#### -scan-bin δ_min δ_max g_ratio

  Repeat the analysis using bin-widths ranging from δ_min to δ_max,
  increasing the size of the bin each time by a factor of *g_ratio*
  (a number > 1.0.  Recommended value 1.1)

#### -scan-gauss  σ_min  σ_max  g_ratio

  Repeat the analysis using Gaussian-widths ranging from σ_min to σ_max,
  increasing the size of the bin each time by a factor of *g_ratio*
  (a number > 1.0.  Recommended value 1.1)

#### -np  num_processors
  Specify the number of processors to use during the calculation.
  (By default this is the typically number of cores of your CPU,
   including hyperthreaded cores.)
  This argument is not available if you did not compile the program
  to enable support for OpenMP.

#### -pg
  If you have already created a density map image, then you can save time by
  omitting the step of blurring the image.  In that case, 
  the brightness of each voxel in the file specified by "-in" argument
  is assumed to contain the local density of the particles of interest.
  (*This density must be in units of 1/Angstroms^3* **not** *1/voxels^3*, 
   The number of particles in the blurred image will be estimated by 
   integrating the brightness of the voxels over the volume of the image.
   If the *-mask* argument was supplied, then the volume of the image will
   omit voxels ignored by the mask.)
   *The use of "-pg" is not necessary or recommended*.)



Note: The **-w**, **-in**, **-out**, **-mask**, and **-mask-select** arguments
are shared with the *filter_mrc* program
and are explained in more detail [here](./doc_filter_mrc.md).


## Details:

1. An image will be created of voxels whose brightness value is 1 if there is
an object centered at that location (as read from the text file), and 0
otherwise.  (If the user specified a mask, then ignore objects which are
located in voxels which lie outside the mask region.)

This image is subsequently blurred by a Gaussian
whose width (σ) equals δ/√(2π).  (Where "δ" is the argument following the
"*-bin-width*" argument.  The volume of this Gaussian is effectively δ^3.)
(Note: If the user specified a mask, then we must normalize the Guassian blur
 near the boundaries of the mask region so that the density there is not
 artificially lower than it should be.  Without this, the density would be
 depleted near the boundary because there are less objects nearby.)

This image should resemble the result of generating 
a histogram using a fixed-sized bin-width.
However the more complex proceedure we use here 
generates a smoother image that is more numerically stable
in tight confined spaces near the cell/mask boundaries
(compared to using a histogram).
 
2. Then all *local* density minima and maxima are found in this image.
   (Minima and maxima located right at the edge of the image boundaries are 
    ignored because spurious fluctuations in density are often located there.)
   From the list of local minima and maxima, 
   the global lowest and highest density is found.

3. The probability that a small ("bin"-sized) region of space (of volume δ^3)
   contains less than this number of particles is calculated using the
   Poisson distribution:
```
   P0=Σ λ^i exp(-λ) / i!
   where λ = the expected number of particles in this region = ρ0 δ^3
        ρ0 = the average density of particles in the entire image = N/V
   and the sum runs from i=0 to floor(ρ δ^3), where ρ is the local
   density of objects at that location (ie in that "bin"-sized region).
```
   The probability that *none* of the bins in the image have a density
   exceeding ρ is calculated using:
```
   P = 1 - (1 - P0)^n
   where n = V / δ^3 is the number of bins in the image
     and V = the volume of the cell being considered
```

