filter_mrc
===========

**filter_mrc** applies a filter to a 3D image, 
and saves the result as a new .mrc/.rec file.
It currently supports
low-pass, high-pass,
thresholding,
brightness inversions,
[generalized](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1)
[Gaussian](https://en.wikipedia.org/wiki/Gaussian_blur),
[Difference-of-Gaussian](https://en.wikipedia.org/wiki/Difference_of_Gaussians),
[Laplacian-of-Gaussian](https://en.wikipedia.org/wiki/Blob_detection#The_Laplacian_of_Gaussian),
3D planar [ridge detection](https://en.wikipedia.org/wiki/Ridge_detection),
and
[3D tensor voting]([example](http://scikit-image.org/docs/dev/auto_examples/features_detection/plot_blob.html)),
and brightness-fluctuation
filters.
Fast [separable](https://en.wikipedia.org/wiki/Separable_filter)
filters are used whenever possible.

**filter_mrc** can also be used for 3D
[scale-free blob-detection](https://en.wikipedia.org/wiki/Blob_detection) 
([example](http://scikit-image.org/docs/dev/auto_examples/features_detection/plot_blob.html)),
local minima-finding, and
[classic watershed segmentation](https://imagej.net/Classic_Watershed),
and the detection and segmentation of **1D curves** and **2D surfaces**
(including **membranes**).
A list of detected objects can be saved to a text file.
Annotated images can be saved to a new .mrc/.rec file.

All filters support "masking".

***Note: The masking feature is broken. Issue opened. -andrew 2019-3-01***

An image *mask* can be used to exclude certain
voxels or regions from consideration.
(Typically these are voxels which have been characterized previously.
The contributions from remaining voxels are normalized, so that objects located
within narrow confined spaces can be detected accurately and without penalty.)
Masks can also be used to give some voxels more consideration
than others during the bluring (filtering) process.  (A.K.A. "weighting".)
You can use a mask to to apply a filter to an image
whose boundaries are smooth and gradual as opposed to jagged and rectangular,



## Usage Examples:

### Example 1
```
# Detect membranes using tensor voting (target thickness ≈ 60.0 Angstroms)

filter_mrc -in tomogram.rec \
  -out membranes_tv.rec \
  -planar 60.0 -planar-best 0.06 \
  -planar-tv 5.0 -planar-tv-angle-exponent 4
```

Note: The computation time will be roughly proportional to the image size
      and the "-planar-best" argument (which ranges from 0 to 1).
      Try this on a small image first. (See example 3 below)

### Example 2

Detect all dark blobs ("minima") between 200 and 280 Angstroms in width.
This corresponds (approximately) to objects which are the size of ribosomes.


```
filter_mrc -in tomogram.rec \
  -blob minima tomogram_blobs.txt 200.0 280.0 1.01

# Now discard the faint, noisy, or overlapping blobs.

filter_mrc -discard-blobs tomogram_blobs.txt tomogram_ribosomes.txt \
  -minima-threshold -70 \
  -blob-separation 0.8

# Finally, display the remaining blobs we want to keep:

filter_mrc -in tomogram.rec \
  -out tomogram_ribosomes.rec \
  -spheres tomogram_ribosomes.txt
```

Note: All of these parameters make reasonable defaults for ribosome detection
      except the "*-minima-thresold*" parameter ("-70" in the example).
      It must be chosen carefully because it will vary from image to image.
      (Strategies for choosing this parameter are discussed below.)



### Example 3

Find the **largest membrane** in an image,
and **generate a closed surface** for that membrane.
(This example requires
 [*PoissonRecon*](https://github.com/mkazhdan/PoissonRecon) and
 [*meshlab*](http://www.meshlab.net).)
**WARNING: Experimental.**
*(This might not work, and these arguments could change in the future.
  -andrew 2019-3-04)*

```
filter_mrc -in tomogram.rec \
  -out membranes_tv_clusters.rec \
  -planar 60.0 -planar-best 0.08 \
  -planar-tv 5.0 -planar-tv-angle-exponent 4 \
  -connect 1.0e+09 -cts 0.707 -ctn 0.707 -cvn 0.707 -cvs 0.707 \
  -select-cluster 1 -planar-normals-file largest_membrane_pointcloud.ply

PoissonRecon --in largest_membrane_pointcloud.ply \
  --out largest_membrane.ply --depth 8
```

Note: This will generate a triangle-mesh file ("*largest_membrane.ply*")
      which can be imported into *meshlab* for smoothing and refinement.
      All of these parameters make reasonably good defaults for membrane
      detection except the "*-connect*" parameter ("1.0e+09" in the example).
      It must be chosen carefully because it will vary from image to image.
      (Strategies for choosing this parameter are discussed below.)



## Arguments:

### Input and Output files

The user must specify the name of the tomogram they wish to process using the
"-in" argument:
```
   -in SOURCE_FILE.mrc
```
(Note: files may also end in ".rec")

Users must specify the name of the new tomogram
(created by applying the filter to the original tomogram)
using the "-out" argument:
```
   -out DESTINATION_FILE.mrc
```


### Voxel Width

*By default, all of the parameters provided by the user
are assumed to be in **physical units**.  (Eg. *Angstroms*, not *voxels*.)
This makes it more convenient to apply the software
to images at different magnifications.
(This includes all of the width parameters used with the
 "-gauss", "-blob", "-dog", "-spheres", and "-template-gauss" filters.)

Consequently, this program needs to know the physical width of each voxel.
By default it will attempt to infer that from the MRC file
*(which are typically stored in units **Angstroms**, not **nm**)*.
If you have a more accurate estimate of the voxel width,
you can specify it using the "**-w**" argument:
```
   -w voxelwidth
```
The "-w" argument will override the voxel widths
specified in the MRC file.

Setting voxelwidth to 1 using "-w 1"
will allow you to enter distance parameters in units of voxels.

*(In the Jensen Lab tomography database, there is typically
 a ~10% difference between the voxel width stored in the
 MRC file, and the voxel width stored at the
 corresponding entry in the tomography database.
 The later is more accurate.)



## A Short Description of Each Filter Type:


The user can select the type of filter they want to use
using command line arguments.


The "**-invert**" filter exchanges bright and dark voxels (while keeping the average voxel intensity and standard deviation the same).


The "**-thresh**", "**-thresh2**", "**-thresh4**", and "**-clip**"
filters are used to clip and rescale voxel intensities,
or to select voxels whose intensities lie within in a certain range.


The "**-gauss**" filter uses a simple (low-pass)
[Gaussian blur](https://en.wikipedia.org/wiki/Gaussian_blur)
filter.
(The "**-gauss-aniso**" filter allows you to apply different amounts of blurring to the X,Y, and Z directions.)
[Generalized Gaussians](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1) are supported.


The "**-dog**"" filter uses a
[Difference-Of-Gaussians](https://en.wikipedia.org/wiki/Difference_of_Gaussians)
filter which can be used as a frequency band-pass filter (for high and/or low frequency removal).
(The "**-dog-aniso**" filter allows you to modify the properties of the filter in the X,Y,Z directions.  Filter sharpness can be customized using the "-exponents" argument.)
[Generalized Gaussians](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1) are supported.


The "**-watershed**" argument will generate a new image which has been
segmented using the
[classic watershed](https://imagej.net/Classic_Watershed)
algorithm.


The "**-fluct**" and "**-fluct-aniso**" filters calculate the magnitude of
the voxel intensity fluctuations relative to their local surroundings.
They can be used to find locations in the image where the brightness
remains relatively constant or fluctuates wildly.  It can also be useful
for characterizing regions within the image that have poor contrast.


The "**-find-minima**" and "**-find-maxima**" will find all of the local
intensity minima and maxima respectively in the image.
The user has the option to discard poor scoring minima or maxima,
or minima and maxima which are too close together.


The "**-blob**", "**-blobd**", and "**-blobr**" filters can be used for [scale-free blob detection](https://en.wikipedia.org/wiki/Blob_detection#The_Laplacian_of_Gaussian).  Objects in the image of various sizes can be detected using this filter.  Depending upon how many sizes are considered, the computation can be slow.  (For a fast and sloppy alternative, you can use the "**-log-d**" filter.)


The "**-planar**" and "**-planar-tv**" filters are used to detect thin, membrane-like structures using
[3D ridge detection](https://en.wikipedia.org/wiki/Ridge_detection)
and 
[3D tensor voting](https://www.cs.stevens.edu/~mordohai/public/TensorVotingTutorial_2007.pdf), respectively.
Voxels belonging to different membranes can be grouped into different clusters
using the "**-connect**" argument.
Voxels belonging to the same membrane can be analyzed and their orientations
can be saved to a file in order to repair holes and perform further analysis
using the "**-planar-orientations-file**" argument.






# Feature detection

## Detecting membranes

### -planar  thickness

When the "**-planar**" filter is selected, a
[3D ridge detection](https://en.wikipedia.org/wiki/Ridge_detection)
filter will be used.
By default, it detects thin membrane-like
structures which are dark on a white background.
The **-planar** argument must be followed by a number indicating
the approximate *thickness* of the membrane-like feature 
you are interested in detecting.
(in physical units, not voxels.  
 Membranes whose thickness within a factor of 2 of this target
 are also likely to be detected.
 Technically, the *thickness* parameter controls the width Gaussian that will
 be initially convolved with the original image, according to σ=thickness/√3.
 For details, see Steger IEEE Trans. Pattern Anal. Machine Int. 1998.)

The *output* of this filter will be bright wherever the derivative of
the brightness varies much more rapidly in one direction than it does
in the two orthogonal directions.
Because this filter depends on second derivatives,
it is prone to detect a large number of spurious fluctuations in brightness
(in addition to the membrane-like structures you are interested in).
Tensor-voting (using the *-planar-tv* argument)
can be used to remove these spurious detected peaks
and improve the signal-to-noise ratio.

*(Techincal details: This will generate an output image whose brightness equals
Ngamma_norm = λ1^2 - λ2^2,
where λ1, abd λ2 are the largest and second-largest eigenvalues of the
hessian matrix, multiplied by σ^2 [in units of voxels].
See Lindeberg Int.J.Comput.Vis.(1998) for details.)*

**WARNING** 
*This filter requires a very large amount of memory
(enough to store at least 8 copies of the original image in 32bit-float mode).*


### -planar-tv  σ_ratio
The "**-planar-tv**" argument enables refinement of (*-planar*) results using
[3D tensor voting](https://www.cs.stevens.edu/~mordohai/public/TensorVotingTutorial_2007.pdf).
It performs a kind of directional blurring which reinforces regions in the image
where detected ridges in the image point in the same or similar directions.
(for a nice visual demonstration, see 
[these slides](http://www.sci.utah.edu/~gerig/CS7960-S2010/handouts/Slides-tensorVoting-Zhe-Leng.pdf).)
The σ_ratio argument is a number larger than 1 which controls the distance
over which this blurring occurs.
(It is typically 5.0.  See technical details below.)
Tensor voting is a very effective method to improve the signal-to-noise
ratio when detecting curves and surfaces.
Tensor-voting refinement is not done by default 
because it can be a very expensive computation.

*(Note: -planar-tv is not an independent filter.
        It enables refinement of existing results from the -planar filter.
        This argument is ignored unless the -planar argument is also used.)*

*(Techincal details: The width of the Gaussian used in the radial-decay
                     function used in tensor voting has a width of
                     σ_tv = σ_ratio ⨉ σ,
                     where "σ" is the bluring used in the ridge detection
                     and it has a value of σ=thickness/√3)*


### -planar-tv-angle-exponent  *n*
The "**-planar-tv-angle-exponent**" parameter (*n*) 
controls the penalty applied to membrane-like regions which are pointing
in incompatible directions.  It is 4 be default.


### -planar-threshold  *threshold*

This will discard voxels whose saliency (after ridge detection)
falls below *threshold*.
This will make subsequent steps in the calculation
**(such as tensor voting)** faster.
Users can choose this threshold by running **filter_mrc** using only
the "-planar" filter argument, and then visualizing the resulting file.
*(In IMOD, you can click on the image and
 select the "Edit"->"Point"->"Value" menu option to
 see the saliency of voxels of interest.)*

In practice, it is often easier to use the **-planar-best** argument
because no intermediate visualization step is required.

### -planar-best  *fraction*

This will discard voxels whose saliency (after ridge detection)
is not in the top *fraction* of voxels in the image
(excluding masked voxels, if applicable).
This will reduce the computation needed for any subsequent steps
in the calculation 
**(such as tensor voting)** 
faster by a factor which is roughly proportional to this number.
The *fraction* parameter should lie in the range from 0 to 1.
(Using *0.1* is a conservative choice, but you can often get away 
 with using lower values.)

If the resulting membranes stored in the "membrane_tv.rec" file
are missing or are incomplete, then increase this number and try again.

### -connect  *threshold*

After membrane detection is performed, 
voxels from different membranes can be grouped into different clusters
using the "**-connect**" argument.
Nearby ridge-like (ie membrane-like) voxels of similar orientation
can be grouped into connected islands.
If the *threshold* parameter is chosen carefully, then these
different islands will hopefully correspond to different objects
(eg. membranes) in the original image.
The "*threshold*" parameter determines how "*membrane-like*"
a voxel must be in order for it to be classsified as a membrane.
It will vary from image to image it must be chosen carefully.

#### Strategies for determining the -connect *threshold* parameter

      To do choose the *threshold* parameter, 
      run membrane-detection (for example using "-planar" and "-planar-tv")
      once in advance without the "-connect" argument
      (as we did in the membrane-detection example).
      Open the file created during that step
      (eg. "membranes_tv.rec") in IMOD.
      Find a place in the image where the saliency (brightnees)
      of the membrane you are interested in is weak.
      Click on voxels located near the weakest point (a.k.a. "junction point",
      or "saddle point") between two different bright blobs
      corresponding to the *same* surface you are interested in.
      These two islands will not be joined unless the *-connect* argument
      is less than the weakest link connecting these two islands.
      (and even then, they might not be joined 
       if the voxel orientations are dissimilar.)
      Select "Edit"->"Point"->"Value" menu option in IMOD to
      see the "saliency" (brightness) of that voxel.
      Do this several times in different places near the junction
      write down the largest "saliency" number.
      Then reduce this number by 20% (ie. multiply it by 0.8).
      This makes a good first guess for the "*-connect*" parameter.

      After using the "*-connect*" argument you can can 
      open the REC/MRC file we created
      (eg "*membranes_tv_clusters.rec*")
      in IMOD, and click on different voxels (using "Edit"->"Point"->"Value")
      to see whether the voxels were clustered correctly into the same object.
      The voxel brightness values in that image should be integers
      indicating which cluster they belong to
      (reverse-sorted by cluster size, starting at 1).

      If some clusters are too big, you can either increase the *threshold*
      value, *or* you can alter increase angular sensitivity by increasing 
      the *-cts*,*-ctn*,*-cvn*, and *-cvs* parameters from 0.707 to, say 0.9.
      (See below.)

      Because it might take several tries, and membrane detection is slow,
      it is a bad idea to try this on a full-sized tomogram image.
      Instead try this on one or more small cropped versions of the image
      near the junction points of interest.
      (You can crop images either by using IMOD, 
       or by using the "crop_mrc" program distributed with *visfd*.)

      Make sure clustering was successful *before* attempting to
      close holes in the surface using *PoissonRecon*.


*(Note: The "-connect" and *-cts*, *-ctn*, *-cvs*, *-cvn* arguments 
  currently have no effect unless the "-planar" argument was also supplied.
  -andrew 2019-3-04)*

### -cts *threshold*
### -ctn *threshold*
### -cvs *threshold*
### -cvn *threshold*

The *-cts*, *-ctn*, *-cvs*, *-cvn* arguments determine how similar the 
orientations of each voxel must be in order for a pair of neighboring voxels
to be merged into the same cluster (ie. the same membrane or same filament).
Threshold values in the range from 0 to 1 are supported.
A *-threshold* value of 0.707 ≈ cos(45°) and corresonds to a
45 degree difference between orientations of neighboring voxels.
In that case, neighboring voxels pointing in directions which 
differ by more than 45°
will be assigned to different clusters (membranes).
This is the default behavior.  (-andrew 2019-3-04)

Most of the time, *-cts*, *-ctn*, *-cvs*, *-cvn* arguments can be omitted.
The difference in meaning between these 4 arguments will be ellaborated on
in the future.  For now it is safe to set them all equal to the same value,
or omit them entirely since the default value of 0.707 works well in most cases
(-andrew 2019-3-04).


### -select-cluster  *cluster-ID*
### -planar-orientations-file  *file_name*

Once clustering is working, you can select one of the clusters using
the "**-select-cluster**" argument.
(Cluster-IDs are assigned in reverse order according to their size, 
 beginning with the largest cluster, which has ID $1$.)
You can create a file which contains a list of voxels locations
(for the voxels belonging to that cluster),
as well as their planar orientations
using the "**-planar-orientations-file**".
The resulting file will be in .PLY format.
This oriented point-cloud file can be used for further processing
(for example for hole-repair using the "PoissonRecon" program).

Note: Voxel coordinates in the point-cloud file are expressed
in physical units (ie Angstroms), not voxels,
Consequently they are not integers (unless the "-w 1" argument was used).



### **-find-minima**, "**-find-maxima**"

Usage:
```
  -find-minima  filename
```
or
```
  -find-maxima  filename
```
The **-find-minima** and **-find-maxima** arguments will
create a file containing the locations of the local minima or maxima
of the voxel brightnesses within the image, respectively.
A "minima" is defined as one or more neighboring voxels
of identical brightness, surrounded by neighbors of higher brightness.
(See below for details about which voxels are considered *"neighbors"*.)

This will generate a text file indicating the location of the minima.
This file is a 5-column ascii text file.
Each file is a 5-column ascii text file with the following format:
```
x1 y1 z1 num_voxels1 brightness1
x2 y2 z2 num_voxels2 brightness2
 :  :  :   :     :
xM yM zM num_voxelsM brightnessM
```
The x,y,z coordinates of each minima or maxima are in the first 3 columns,
followed by the number of voxels in the maxima or minimaa (which is usually 1), and finally
it's "score" (which is the brightness of the voxel at that location).
*(This format is nearly identical to the format used by the
  "-blob",  "-blobr",  "-blobd" and "-spheres" arguments.)*
When a minima or maxima contains multiple voxels of identical brightness,
the position of one of the voxels among them is chosen arbitrarily.

You can use the "**-out filename.rec**" to create 
an image showing which voxels belong to each minima or maxima.
(In that image, the voxel brightness will be an integer indicating 
 the minima or maxima to which the voxel belongs, if any, starting at 1.  
 When both minima and maxima are
 displayed, the minima are represented by negative integers.)


*Note:* When searching neighboring voxels, all 26 neighbors (3x3x3-1)
are considered by default.
(You can use the **-neighbor-connectivity nc** argument to skip
 3D corner voxels by setting nc=2, and also 2D corners by setting nc=1.
 This may increase the number of spurious local minima and maxima discovered.)

*Note:* By default, local minima and maxima which lie on the boundaries
of the image (or on the boundary of the mask region), are considered.
To ignore these extrema (ie, to consider only voxels which are surrounded
by other voxels) use the "**-ignore-boundary-extrema-**" argument.

#### Non-max suppression

The "**-find-minima**" and "**-find-maxima**" arguments 
are sometimes used together with the following arguments:
```
  -minima-threshold   threshold
```
or
```
  -maxima-threshold   threshold
```
This allows users to discard poor scoring minima (or maxima), 
whose "scores" (brightnesses) do not fall below (or above) "threshold".

It may also be useful to discard minima which are too close together.
This can be done using a combination of these two arguments:
```
  -sphere-radius    radius
```
and
```
  -blobr-separation  ratio
```
The "**-sphere-radius** argument allows you to assign a radius to each minima
(or maxima).  When used together with "**-blobr-separation**", it means that
if a pair of minima (or maxima) lie within the sum of their
radii times ratio (*2\*radius\*ratio*), 
then the poorer scoring minima will be discarded.




## Blob detection


### -blob,  -blob-r,  -blob-s

The "**-blob**", "**-blob-r**", and "**-blob-s**" arguments are used for
[Scale-Free Blob-Detection](https://en.wikipedia.org/wiki/Blob_detection).
When this is used, the program will apply a LoG filter to the image
multiple times using a range of Gaussian widths (σ) (specified by the user)
in an attempt to determine the correct size (scale) for the relevant objects
in the image automatically.  (As such, this operation is comparatively slow.)
A list of local minima and maxima in *X,Y,Z* space (and scale-space)
will generated and saved in a file, using the method described in:
Lindeberg,T., Int. J. Comput. Vision., 30(2):77-116, (1998)

The "**-blob**", "**-blob-r**" and "**-blob-s**" filters are followed 
by 5 arguments (whose meaning depends on the filter selected):
```
  -blob   type  file_name  d_min  d_max  gratio
  -blob-s type  file_name  r_min  r_max  gratio
  -blob-r type  file_name  σ_min  σ_max  gratio
```
The "**type**" argument must be either "*minima*", "*maxima*", or "*all*".
*(It's meaning is explained in detail below.
  For nearly all CryoEM images, it is safe to use "minima".)*

The "**file_name**" argument is the name of a text file which will store
a list of detected blobs.  The format of this file is explained below.

If "**-blob**" is selected, then the **d_min** and **d_max** parameters
(3rd and 4th parameters), specify a range of diameters
of the objects that you wish to detect.
(Simlarly, 
"--blob-r" allows the user to specify blob sizes in terms of their radii.)
Either way, a LoG filter will be applied to the image using a series of
different Gaussians whose widths (σ) vary between "**σ_min**", and "**σ_max**".
In this series, each successive Gaussian width (σ) will be larger than the
previous one by a factor of (no larger than) "**gratio**",
a number which should be > 1.
(**1.01** is a safe choice for **gratio**,
 but you can speed up the calculation by increasing this parameter.
 Values as high as 1.1 are not uncommon.)
(Note that: *σ_min*, and *σ_max* are equal to *d_min/(2√3)* and *d_max/(2√3)*,
 or *r_min/√3* and *r_max/√3*, respectively.
 If you prefer, you can use the "*-blob-s*" argument to directly specify the
 range of Gaussian widths you wish to use"*σ_min*", and "*σ_max*".)
After applying the LoG filter, if a voxel in the image is a (scale-adjusted)
local maxima (or minima) in *x,y,z,* and *σ*,
*then* its location and size will be recorded in a file.
"**file_name**" is the name of a file which will store these
locations of all the blobs detected in the image
as well as their sizes and scores (see below).
These detected blobs are either local minima or maxima in
X,Y,Z,[scale-space](https://en.wikipedia.org/wiki/Blob_detection#The_difference_of_Gaussians_approach).
By default, two files will be created, named
*file_name.minima.txt* (for storing local minima), and
*file_name.minima.txt* (for storing local maxima).
*(Note: If the "-minima-threshold" or "-maxima-threshold"
 argument is specified, then only one file is generated.)*
Each file is a 5-column ascii text file with the following format:
```
x1 y1 z1 diameter1 score1
x2 y2 z2 diameter2 score2
 :  :  :   :     :
xM yM zM diameterM scoreM
```
...where
"**M**" is the number of blobs (maxima or minima) which were detected,
On each line of that file,
**x,y,z** are the coordinates for the blob's center,
**diameter** is the diameter of the blob (which equals (2√3)σ),
and **score** is the intensity of that voxel after
a (scale-adjusted) LoG filter of that size was applied to it.
For *minima* type blobs, the list is ordered
from low score (most negative) to high score score
For *maxima* type blobs, the list is ordered from high score to low score.
These blobs can be visualized using the "**-spheres**" argument (see below).

#### Blob types

The first argument indicates the **type** blob that the user wants to detect.
(The **type** is the 1st argument following the "-blob" arument.)
It must be one of the following choices:

|   type   | meaning |
| -------- | ------- |
| *minima* | detect dark blobs on a bright background |
| *maxima* | detect bright blobs on a darker background |
|  *all*   | detect both dark and bright blobs |

Note: If **type** is set to *all*, then two different files will be generated
whose names will end in ".minima.txt" and ".maxima.txt", respectively.


#### Automatic disposal of blobs (non-max suppression)

***By default, all minima and maxima are reported during blob detection ***
(no matter how insignificant).
Unfortunately, the blob-detector will typically an enormous number of 
blobs in a typical image.
(Running a blob detector on some Electron-Cryotomography images can result
 in tens of millions of spurious blobs.)
The *vast majority of these detected blobs* are either
not relevant due to noise in the source image, or are overlapping.
You can discard these blobs during the blob detection process
using the arguments below.
However, because blob-detection is slow,
it is recommended that you save all of these blobs to a file.
Later on you can decide which of these blobs are meaningful
by running "filter_mrc" again with the "**-spheres-nonmax**" argument.
This will discard the meaningless blobs.
Then you can run "filter_mrc" a third time using the "**-spheres**" 
argument to visualize the remaining blobs you did not throw away.
Sometimes several iterations of non-max suppression followed by visualization
are needed before you get visually pleasing results.

Details will be provided elsewhere below explaining how visualize these blobs
(using the "**-spheres**" argument),
as well as how to discard low-quality and overlapping blobs
(using the "**-spheres-nonmax**" argument).


#### Recommendation:

Blob-detection is computationally expensive,
but it only has to be performed once.
Typically, users will perform blob detection with the default (permissive)
score thresholds, so that they keep most of the minima and maxima.
Then later, they will run **filter_mrc** with the
"**-spheres**" and "**-minima-threshold**" (or "**-maxima-threshold")
arguments several times with different thresholds to decide which of
these blobs are worth keeping.


## Morphology

### -watershed-minima

If the "**-watershed**" argument is selected, the image will be segmented
using the
[classic watershed](https://imagej.net/Classic_Watershed)
algorithm.
Afterwards, each voxel in the image will be assigned to an integer 
indicating the local-minima basin to which it belongs.
By default voxels which lie at the "ridges" 
(ie., at the boundary of multiple different drainage basins),
are assigned a value of 0.
(You can prevent this using the
 "**-watershed-hide-boundaries**" argument.)

*Note:* These minima are sorted according to depth and begin at 1, not zero.
(Hence voxels in the first minima's basin, ie, the global minima basin, will be
 assigned to intensities of 1.  Voxels in the next deepest basin, 2, etc...)

*Note:* The locations of the corresponding local minima can be found by
invoking the program on the same image using the "**-find-minima**" argument.
(They are also sorted according to minima depth.)

*Note:* When searching neighboring voxels, all 26 neighbors (3x3x3-1)
are considered by default.
(You can use the **-neighbor-connectivity nc** argument to skip
 3D corner voxels by setting nc=2, and also 2D corners by setting nc=1.
 This may increase the number of spurious basins discovered.)

*Note:* Oversegmentation can be reduced by performing a Gaussian blur 
on the image to remove small, insignificant local minima beforehand.
(This algorithm is usually only applied to images that have been smoothed
or filtered in advance.)


### -watershed-threshold  threshold

If the "**-watershed-threshold**" argument is also supplied, then voxels
whose intensity exceeds *threshold*, will be assigned to
the *number of basins (whose depth lies below this value) + 1*.
(Since this is an impossible value,
 these will be the brightest voxels in the image.
 IE. They will have the highest intensity.)

### -watershed-maxima

This performs watershed segmentation starting from maxima instead of minima.


### -watershed-hide-boundaries

By default, voxels which lie on the border between two or more different basins
will be assigned a value of 0.
When this argument is selected, the basin to which a boundary voxel is assigned
will be chosen (arbitrarily) from among the neighboring basins.


### -watershed-boundary  label

By default, voxels which lie on the border between two or more different basins
will be assigned a value of 0.
When this argument is included, these voxels will have intensities
which are assigned to *label* instead.  (This parameter is a number.)




## General filters:


### -gauss and -ggauss
The **-gauss** and **-gauss-aniso** arguments must be followed by one or more numbers specifying the width of the Gaussian filter to apply to your image:
```
   -gauss  σ
```
or
```
   -gauss-aniso  σ_x  σ_y  σ_z
```
In either case, the original image is convolved with:
```
   h(x,y,z) = A*exp(-0.5 * ((x/σ_x)^2 + (y/σ_y)^2 + (z/σ_z)^2))
```
When *-gauss* is used, *σ_x = σ_y = σ_z = σ*

If **-ggauss** or **-ggauss-aniso** is used,
then the image is instead convolved with
[a generalized gaussian function](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1):
```
   h(x,y,z) = A*exp(-r^m)
    where r = √((x/σ_x)^2 + (y/σ_y)^2 + (z/σ_z)^2))
```
The width of the Gaussian (ie, the σ_x, σ_y, σ_z arguments) should be specified in units of physical distance, not in voxels.
(The A coefficient will be determined automatically by normalization.)
By default **m** *= 2*,
however this can be overridden using the optional "**-exponent m**"
*Note:* that these *σ*, *σ_x*, *σ_y*, *σ_z*, parameters are a factor
of *√2 times larger* than the corresponding *σ* ("sigma") parameter
traditionally used in the
[Gaussian function](https://en.wikipedia.org/wiki/Normal_distribution).
(For your convenience, you can plot these functions for various parameters using
the "draw_filter1D.py -ggauss A s m" python script which is located in the same directory.
In this case, the "s" parameter is specified in voxels, not physical distance.)
*Note:* The calculation is
[fast](https://en.wikipedia.org/wiki/Separable_filter)
if you use the default exponent of 2.
*Changing the exponent will slow down the filter considerably.*

The filter is truncated far away from the central peak at a point which is chosen automatically according the σ, σ_x, σ_y, σ_z parameters selected by the user.  However this can be customized using the "-truncate-threshold" and "-truncate" arguments if necessary (see below).

### -dog
When the "**-dog**" or "**-dog-aniso**" filter is selected, a
[Difference-of-Gaussians](https://en.wikipedia.org/wiki/Difference_of_Gaussians)
filter will be used.
The **-dog** and **-dog-aniso** arguments must be followed by numbers indicating the width of the Gaussians you wish to convolve with your image.
(Equivalently, the numbers indicate the band of frequencies you wish to keep from your original image.)
```
  -dog  a  b
```
or:
```
  -dog-aniso  a_x  a_y  a_z  b_x  b_y  b_z
```
In either case,
The original image will be convolved with the following function:
```
  h(x,y,z) =   A*exp(-0.5 * ((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))
             - B*exp(-0.5 * ((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))
```
When *-dog* is used,
*a_x = a_y = a_z = a*
and
*b_x = b_y = b_z = b*

### -log,  -log-r, -log-d,  -log-aniso


When the "**-log**", "**-log-r**", and "**-log-d**", or "**-log-aniso**"
arguments are selected, a
[Laplacian-of-a-Gaussian (LoG)](https://en.wikipedia.org/wiki/Blob_detection)
filter is applied on the source image.
(See implementation details below.)
The Laplacian-of-a-Gaussian filter applies a Gaussian blur to the image
(whose width is specified by the user),
and then calculates the Laplacian of the result.
Features in the image of various sizes can be emphasized
using this kind of filter, and the result can be saved as a new image using
"**-out** filename.mrc".
The "**-log**", "**-log-aniso**", "**-log-r**", and "**-log-d**"
arguments described here are typically only chosen when the user already 
knows the approximate size of the features they want to emphasize
in the image.
```
    -log  σ
```
When "**-log  σ**" is used, a LoG filter will be applied
to the image using a Gaussian width (σ).
Alternatively, if the user wishes to specify the actual size 
(the radius or diameter) of the objects they want to emphasize, 
then they can (equivalently) use 
```
    -log-r radius
```
or
```
    -log-d diameter
```
instead.
There is also an anisotropic version of this filter as well:
```
    -log-aniso  σ_x  σ_y  σ_z**
```
The new image which is generated will *tend* to have have bright and dark spots 
wherever bright and dark objects of the corresponding size can be found.
This image can be searched for (local) minima and maxima
by running this program again on the new image using the 
"**-find-minima**", "**-find-maxima**" and "-blobr-separation" arguments.
These locations can be saved to a file and then a new anotated image can be
created using the "**-spheres**" argument.  This provides a fast, sloppy, 
and unselective way to search for features in an image.
However scale-free 
[**blob detection**](https://en.wikipedia.org/wiki/Blob_detection)
is a more robust (and computationally expensive) way to detect objects 
of a given size within an image.
"Blob detection" is discussed elsewhere in this document, and should not be
confused with the discussion here.

*(Implementation note: The LoG filter described here is approximated internally 
 using a DoG filter. You can control the accuracy of this approximation using
 the "-dog-delta δ" argument, which is 0.02 by default.  Using smaller values
 of δ can improve the approximation, but could lead to spurious artifacts.)*



### -fluct,   -fluct-aniso
Usage:
```
  fluct  r
```
The "**-fluct**" filter calculates the magnitude of
the voxel intensity fluctuations relative to their local surroundings.
This is useful for finding regions in an image with nearly uniform brightness
and regions whose brightness fluctuates wildly.
*r* is the effective search radius.

To compute this, near every voxel, a Gaussian-weighted average
intensity is computed with the same effective radius, *r*.
A Gaussian-weighted squared-difference between all the nearby voxel intensities
and this local average intensity is also computed.
A new image is generated whose voxel intensities equal this
local average squared difference intensity of nearby voxels.

The "**-fluct-aniso**" variant allows the user to control the width
of the Gaussian independently in the x,y,z directions:
```
  fluct-aniso   r_x  r_y  r_z
```
#### Implementation detail:
 The width of the Gaussian, σ, is chosen so that the effective
 number of voxels under the 3D Gaussian,
 (which can be interpreted as the 3D integral of the unnormalized 3D Gaussian,
 (2\*π\*σ^2)^(3/2))
 equals the volume of a sphere of radius r
 (4/3)\*π\*r^3.
 This yields σ=(9π/2)^(-1/6)r.

 It might seem more straightforward to simply consider the fluctuations in
 brightnesses of all of the voxels which lie within a fixed radius from
 each target voxel.  However using Gaussian weighted averages instead
 accomplishes essentially the same goal and is much more
 [computationally efficient](https://en.wikipedia.org/wiki/Separable_filter).

 If you prefer, instead of using a Gaussian, you can instead perform the
 fluctuation calculation over a hard spherical region by using the
 "-exponent n" argument.
 (This parameter is 2 by default.  Changing it will slow down the calculation.)
 As n is increased, the region over which the fluctuations in brightness
 are calculated becomes more and more like a uniform sphere with radius σ.
 (So if you do this, be sure to multiply your radius argument, r,
  by (9π/2)^(1/6)≈1.5549880806696572 to compensate.)




## Thresholding, Clipping and Rescaling (Intensity Map Filters):

Thresholding filters provide a way to rescale (map) the
intensity of all of the voxels in the image to new intensities.
For example, you could multiply these intensities by a constant,
and/or clip them within a specified range.
Thresholding operations can be used to insure that the intensities
(ie. brightnesses) for every voxel in the image lie between 0 and 1.
They can also be used to select only voxels whose intensities lie within
a narrow range of interest.

*Note:* All threshold operations are performed *after* normal filtering operations have been applied (such as -gauss, -dog, -dogg, -fluct, or -log filters).

### -invert

This filter replaces bright voxels with dark ones, and visa versa,
while keeping the average voxel brightness the same.
(To calculate the new voxel intensity,
 the *difference* between the original voxel intensity
 and the average intensity is *subtracted* from the average intensity.)
Inverting an image can be useful to prepare files which you plan
to read with other programs (such as ChimeraX and IMOD).



### -rescale m b

Rescale the voxel intensities (multiply the brightnesses by *m*)
and add an offset (*b*).
Afterwards *new_intensity=m\*old_intensity+b*.
    
```
         output
        intensity
           /|\        /
            |        /
            |       /     /|
            |      /   m / |
            |     /     /__|
            |    /                      
  <---------|---/------------------->  input
        /|\ |  /                       intensity
      b  |  | /
        \|/ |/
            /
           /|
          / |
         /  |/
```


### -rescale-min-max outA outB

If you use the **-rescale-min-max** argument, then
image voxel intensities in the final image
(after all other processing has been applied)
will be shifted and rescaled so that the
minimum and maximum voxel intensities will be *outA* and *outB*.

*(Note: As of 2018-6-30, the "Black" and "White" controls in *IMOD*
  report values between 0 and 255, even if the actual voxel brightnesses
  in the file are floating point numbers (between 0 and 1, for example).
  Ignore those numbers.
  You can find the brightness of a particular voxel by left-clicking
  on that voxel (to move the yellow-cursor)
  and then by selecting the "Edit"->"Point"->"Value" menu option within IMOD.
  You can also use the "histogram_mrc.py" script to
  see if your intensity values in the image lie in the range you expect.)*


### -clip a b

The "**--clip a b**" argument will
restrict the range of voxel intensities to a range of your choosing.
This can be useful to eliminate extremely dark or bright voxels from
your image.
(Such spurious voxels can make it difficult to
 detect faint objects in the original image.)

```
-clip 0.48 0.52
```
The resulting voxels intensities will lie in the range from 0.48 to 0.52.
```
 output
 intensity
    /|\                              _________________
0.52 |                           _.-'                 
     |                       _,-'                 
     |                   _,-'            
0.48 |________________,-'                     ________\ input
                      ^              ^                / intensity
                     0.48          0.52
                  (thresh_a)    (thresh_b)
```

*Note: If instead, you want the output intensities to range from*
      **0 to 1**,
      (for example), then use the more general*
      "**-thresh2**"
      *argument. (See below.)*

*Note: Unfortunately, intensity values can be difficult to guess.
       When choosing thresholds, the
       [histogram_mrc.py program](doc_histogram_mrc.md)
       can be useful.
       It displays the range of voxel intensities in an image.
       Alternatively you can use the "**-cl a b**" argument.*

### -cl
```
-cl  a  b
```
The "**-cl**" argument is similar to the "**-clip**" argument, however it
allows you to specify the minimimum and maximum intensity parameters
in units of σ, where σ is the standard deviation of the brightness values
present in the image.  Typical usage:
```
-cl  -2.5 2.5
```
This will clip voxel intensities which are either 2.5 standard-deviations
below or above the average intensity value, respectively
**Note:** You can use the "**-mask**" argument to exclude certain voxels
       from the clipping and thresholding operations.
       When using "**-mask**" together with "**-cl**", then
       these voxels will also be excluded from consideration when
       calculating the average and spread(σ) of voxel intensities.


### -thresh  threshold

If the "-thresh" argument is passed, it must be followed by a number
("thresh01").
Voxels with intensities *below* this number will be replaced 0,
and voxels *above* this number will be replaced with 1 (by default, see below).
For example:
```
-thresh 0.5
```
results in:
```
          output
         intensity
            /|\
             |
   outB   -->|                       __________________________\
(usually 1)  |                      |                          /
             |                      |                 
             |                      |            
   outA   -->|______________________|  ________________________\ input
(usually 0)                         ^                          / intensity
                                   0.5
                               (threshold)

```

*Note: When choosing thresholds, the
       [histogram_mrc.py program](doc_histogram_mrc.md)
       can be useful.
       It displays the range of voxel intensities in an image.*


*Note: You can use the "-rescale-min-max outA outB" argument to scale the
       resulting voxel intensities from outA to outB (instead of from 0 to 1).*




### -thresh-range  range_a  range_b
    Select a range of voxels whose intensities fall within
    the range from *range_a* to *range_b*.
    Each voxel's new intensity will be a function of its former intensity.
    If the range_a < range_b, then the thresholding function will be

```
         output
        intensity
           /|\
            |
            |
   outB  -->|              ___________________
(usually 1) |             |                   |
            |             |                   |
   outA  -->|_____________|                   |________________\ input
(usually 0)             range_a             range_b            / intensity
```

*Note: You can use the "-rescale-min-max outA outB" argument to scale the
       resulting voxel intensities from outA to outB (instead of from 0 to 1).*

*Note: This command is equivalent to "-thresh4 range_a range_a range_b range_b"*



### -thresh2  thresh_a  thresh_b

If the "-thresh2" argument is passed, then it must be followed by 2 numbers.
For example:
```
-thresh2 0.48 0.52
```
In this case, the resulting voxel intensities in the range
from *thresh_a* to *thresh_b* will be scaled between 0 and 1
and clipped above and below.
Graphically, the threshold filter resembles a step function
with a smooth ramp between *thresh_a* and *thresh_b*:

```
           output
         intensity
            /|\
             |
   outB   -->|                               _________________
(usually 1)  |                           _.-'                 
             |                       _,-'                 
             |                   _,-'            
   outB   -->|________________,-'                     ________\ input
(usually 0)                   ^              ^                / intensity
                             0.48          0.52
                          (thresh_a)    (thresh_b)
```

*Note: When choosing thresholds, the
       [histogram_mrc.py program](doc_histogram_mrc.md)
       can be useful.
       It displays the range of voxel intensities in an image.*

*Note: You can use the "-rescale-min-max outA outB" argument to scale the
       resulting voxel intensities from outA to outB (instead of from 0 to 1).*


### -thresh4  thresh01_a  thresh01_b  thresh10_a  thresh10_b

Sometimes it is useful to select a ***narrow range of voxel intensities***.
*filter_mrc* allows you to do this using the **"-thresh4"** argument.
The **"-thresh4"** must be followed by 4 numbers in order, for example:
```
 -thresh4 0.3 0.4 0.5 0.6
```
The resulting image voxels will be scaled
between 0 and 1 according to the following function:

```
         output
        intensity
           /|\
            |
   outB  -->|                 ________________                
(usually 1) |             _,-'                `-._
            |         _,-'                        `-._
   outA  -->|______,-'                                `-._______\ input
(usually 0)        ^         ^                ^         ^       / intensity
                  0.3       0.4              0.5       0.6
             thresh_01_a  thresh_01_b    thresh_10_a  thresh_10_b
```

*Note: When choosing thresholds, the
       [histogram_mrc.py program](doc_histogram_mrc.md)
       can be useful.
       It displays the range of voxel intensities in an image.*


*Note: You can use the "-rescale-min-max outA outB" argument to scale the
       resulting voxel intensities from outA to outB (instead of from 0 to 1).*




### -thresh-gauss x0 σ
Select a range of voxels whose intensities fall within a Gaussian
centered around *x0* with standard deviation σ.
    
```
         output
        intensity
           /|\
            |
            |                       _---_
(usually 1) |                     ,'  :  `.
            |                    /    :<-->\
            |                _.-'     :  σ  `-,_
    outA -->'------------''''         :         ````------------> input
(usually 0)                          x0                          intensity
```

*Note: You can use the "-rescale-min-max outA outB" argument to scale the
       resulting voxel intensities from outA to outB (instead of from 0 to 1).*



## Annotation of detected objects in an image

## -spheres filename

The "**-spheres**" argument does not perform a filtering operation
on the original image.
Instead it reads a text file containing between 3 and 5 columns (eg "filename")
indicating the location, size, and brightness of a series of points in space.
After reading that file, a new image will be created
(the same size as the input image)
with each blob represented (by default) by a hollow spherical shell
that surrounds the point of interest.

Each line in the file corresponds to a different sphere.
The first three numbers on each line are assumed to store the x,y, and z
coordinates of the center of that sphere.
If the file contains only 3 columns, the "spheres" will be only 1-voxel wide
by default.
If the file contains a 4th column, then it is assumed to store the diameter of
the sphere (in physical units, not voxels).
(Either way, the sphere size can be overridden using
 the "**-sphere-diameter d**" argument.)
If the file contains a 5th column, it is assumed to represent the brightness
of that sphere.
(This can be overridden using the "**-sphere-foreground brightness**" argument.)
Incidentally, the format of this file matches the formmat of the text file
generated during blob detection (using the "**-blob**" argument).

The thickness of the shell can be controlled using the
"**-spheres-shell-ratio ratio**" argument.
Setting the "ratio" to 1 results in a solid rather than hollow sphere.
(The default shell thickness is 0.08.
 The minimum width of the shell is 1 voxel wide.)

By default, these spherical shells will be superimposed upon the
original image (whose voxel's brightness will be shifted and scaled
so that the voxels remain easy to see in the presence of the spherical shells).
The average *brightness* and *contrast* of these background voxels
is controlled by the "**-sphere-background brightness**", and the
"**-sphere-background-scale ratio**" arguments, respectively.
(Using "**-sphere-background-scale 0**" makes this background image invisible.
 The default scale ratio is 0.333.)

When displayed in IMOD,
blobs with good scores typically appear as black or white spheres,
and grey spheres correspond to blobs with poor scores.
The score of a given blob
can be queried in IMOD by clicking on a voxel somewhere on the spherical shell
and selecting the "Edit"->"Point"->"Value" menu option.)




## -spheres-nonmax  
```
   -spheres-nonmax  orig_blobs.txt  selected_blobs.txt
```
When *filter_mrc* is run using the "**-spheres-nonmax**" argument,
no image processing is done on the source image.
Instead, *filter_mrc* will read the list of blobs reported in a file
(created by running *filter_mrc* with the *-blob* argument),
and discard blobs with poor scores, or blobs which overlap with other blobs.
The new list of blobs can then be visualized later by running
*filter_mrc* again using the "**-spheres**" argument.

Note that the "**-spheres-nonmax**" argument by itself does nothing.
You must run it with other arguments telling it which blobs you want to discard.
Typically you would run *filter_mrc* together with the "**-spheres-nonmax**",
**-maxima-threshold** (OR **-minima-threshold**), AND the
**-max-volume-overlap** arguments, simultaneously.  (See below for details.)

#### Automatic disposal of poor scoring blobs
```
   -minima-ratio  ratio_min
   -maxima-ratio  ratio_max
```
The *-maxima-ratio* argument allows you to discard maxima whose score 
is less than *ratio_max* times the blob with the highest score,
(where *ratio_max* is a number between 0 and 1.
 *-minima-ratio* behaves in a similar way.)
Alternatively, if you prefer, you can specify the exact threshold used
to discard local minima and maxima (instead of specifying ratios),
using the following arguments:
```
   -minima-threshold  thresh_min
   -maxima-threshold  thresh_max
```
...where "thresh_min" and "thresh_max" are the thresholds below which,
and above which, the blob score must lie in order not to be discarded,
respectively.
*(Note: You can either use thresolds or ratios, but you cannot mix both.)*


#### Automatic disposal of overlapping blobs (non-max suppression)
There are several arguments you can use to discard overlapping blobs
which are typically invoked together with "**-spheres-nonmax**"
```
   -max-volume-overlap  fraction
   -max-volume-overlap-small  fraction
   -radial-separation  fraction
```

Whenever two blobs overlap, the one with the better score
(ie. higher score for maxima, and lower score for minima)
is superimposed upon the one with a poorer score.
If they overlap too much, the user can request that the 
one with a poorer score is discarded.
Again, this does not happen by default.

Disposal of overlapping blobs can be accomplished using the
"**-max-volume-overlap fraction**" argument.
If this argument is specified, then the poorer scoring blob is discarded only
when the volume overlapping region between the two spheres exceeds *fraction*
multiplied by the volume of the *larger* sphere.
(A larger *fraction* means that more overlap is permitted.  
Setting it to 1.0 disables overlap prevention.)
Additionally, you can specify the maximum amount of overlap with the *smaller*
of the two spheres using the "**-max-volume-overlap-small fraction**" argument.
(See example below.)

Alternatively, you can use the
"**-radial-sepration fraction**" argument
to discard blobs if the distance is less than **fraction**
multiplied by the sum of their radii.
The "**fraction**" parameter is a number between 0 and 1.
(A larger number means less overlap is permitted.
 Setting *fraction=0.0* allows complete overlap and disables this argument.)
*Note:* You can *either* discard blobs based on *overlapping volume*
or *overlapping radii*, but not both.

##### Example: Hierarchical blob detection
Sometimes a small spherical blob may reside *within* a significantly larger 
sphere. If the small sphere's volume is less than half of the larger sphere
(for example), you may want to keep both spheres.
To prevent the small sphere from being discarded, you can use
**-max-volume-overlap 0.5** together with
**-max-volume-overlap-small 1**.





#### Specifying the radius or Gaussian-sigma parameters for the objects of interest:

The "**-blob-r**", "**-blob-s**", (and "**-log-r**" "**-log**") arguments
are variants of the "**-blob**" (and "**-log-d**") argument.
Their parameters are specified by the approximate radius(r≈σ√3)
or Gaussian width (σ) of the objects
that you wish to detect within the image (instead of the object's diameter).



####  Visualizing blobs using the "**-spheres**" argument

Again, after blob-detection AND nonmax-suppression *(discussed elsewhere)*,
you can visualize the resulting blobs using the "**-spheres**"
and "**-out**" arguments.







## Masking

The optional "-mask", "-mask-select", and "-mask-out" arguments allow you to
ignore certain voxels from the source image (tomogram).

Using "masks" allows the user to perform the filtering considering only a
subset of voxels in the tomogram (ie, the "mask").

Masks can also be used to give some voxels more consideration
than others during the bluring (filtering) process.  (A.K.A. "weighting".)
This effectively makes it possible to apply a filter to an image
whose boundaries are smooth and gradual as opposed to jagged and rectangular,
slowly fading from 1 to 0 (near a curve-shaped boundary, for example).


```
   -mask  file.mrc
```
The argument following the
"-mask" argument should be the name of a tomogram file (MRC format) of the
same size as the input tomogram.
During the filtering process, only voxels belonging to the tomogram
with non-zero voxel values in the "mask" tomogram will be considered,
and these remaining voxels will be multiplied by the mask's voxel value
before filtering.
(Usually every voxel in the mask tomogram is a number between 0 and 1.)
 This allows you to use soft masks with gradual boundaries.)
(This can be overridden using the "-mask-select" argument.)
*(Note:
 Certain filters like "-gauss", which have unit norm, will rescale the result of
 the filter by the sum of the mask values in the region considered.)*
```
   -mask-out  intensity_value
```
If the "-mask-out" argument is specified, then voxels outside the mask will
be assigned to the intensity following this argument.  Otherwise they are
assigned to 0 by default.  
*(Note that by default, voxel brightnesses
are rescaled between 0 and 1 before and after filtering.  
This means that "-mask-out 0" assigns these voxels same brightness as the lowest brightness voxels in the image.
If you are using "-rescale-min-max", then
"-mask-out 1" assigns them to the brightest voxels in the image.)


```
   -mask-select  intensity_value
```
If the "-mask-select" argument is specified, then instead of considering all
voxels with non-zero values from the "mask" image,
only voxels whose mask intensity equals the number following
this argument will belong to the mask.
All other voxels will be ignored during filtering.
*(Note: This disables "soft" masking.  In other words,
  during the process of filtering, all selected voxels will be weighted equally
  during filtering, and all others will be completely ignored.)*




## Miscellaneous

### -np  num_processors
  Specify the number of processors to use during the calculation.
  (By default this is the typically number of cores of your CPU,
   including hyperthreaded cores.)
  This argument is not available if you did not compile the program
  to enable support for OpenMP.


### Filter Size
```
   -truncate-threshold threshold
```
This specifies the number of voxels in the filter which will be convolved
with the image.
The filter window width is extended in the X, Y, and Z directions until
the filter intensity decays to less than **threshold** (a fraction),
compared to its peak value.
(For Difference-of-Gaussian filters, both Gaussians must decay to this value.)
If unspecified, this parameter is set to 0.02.
(This overrides the "-truncate-ratio" argument.)

```
   -truncate Wx Wy Wz
```
If you prefer to specify the size of the filter window manually, you can
use the "-window-ratio" argument.
This argument specifies the size of the 3-D filter box in x,y, directions
(in units of the **a** and **b** distance parameters, or **σ** ("sigma")
 parameters which appear in the formulas above).
This overrides the "-truncate-threshold threshold" argument.

*(Note: Incidentally, the product of these 3 numbers, Wx*Wy*Wz, is proportional to the running time of the filter.
However when using the "-gauss" or "-gdog" filters with default exponent settings,
the running time is proportional to the sum of these numbers, Wx+Wy+Wz.
Keep this in mind when specifying filter window widths.)*


### Distance Units: Angstroms or Nanometers
```
   -a2nm
```
When reading MRC files, this argument will
multiply the voxel width
specified in the MRC file by 0.1.
This allows you to specify your distance parameters in units of nm
instead of Angstroms.
(The voxel width in most tomograms is in Angstroms by default.)

This argument has no effect if the "-w" argument is used.



### -distance-points  coordinate_file

The **-distance-points** argument reads a file containing a list of 3D 
coordinates and generates an image whose voxel intensities equal
the *distance* to the nearest point.
(This should not be confused with the
 [distance transform](https://en.wikipedia.org/wiki/Distance_transform).
 The size of the image matches the size of the image
 specified with the the "**-in**" argument.
 The points should lie within the boundaries of the image.)
Example:

```
   -distance-points coordinates.txt
```
The file should contain at least 3 numbers per line
(the x,y,z coordinates of that point, one point per line).
The x,y,z coordinates are assumed to be in physical units
(ie Angstroms) not voxels, *unless* the "**-w 1**" argument was also used.
*(WARNING: As of 2019-1-11, the calculation is quite slow.
  I've been too lazy to bother optimizing this for speed.)*


### -dogg
*Depreciation warning: This type of filter is unlikely to be useful
                       and may be removed in the future.*

If the **-dogg** or **-dogg-aniso* arguments are selected,
then the image will instead be convolved with the following function:
```
  h(x,y,z) = A*exp(-r_a^m) - B*exp(-r_b^n)
   where r_a = √((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))
     and r_b = √((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))
```
The width of the Gaussian (the a_x, a_y, a_z, b_x, b_y, b_z arguments) should be specified in units of physical distance, not in voxels.
The A and B coefficients will be automatically chosen
so that each Gaussian is normalized.
(The sum of h(x,y,z) over x,y,z is 0.
 The sum of the A and B terms over x,y,z is 1.)
The filter is truncated far away from the central peak at a point which is chosen automatically according the shape of the filter selected by the user.  However this can be customized using the "**-truncate-threshold**" and "**-truncate**" arguments if necessary (see below).
![A comparison of Difference-of-Generalized-Gaussians filter weights](./images/example_dogxy_w=2.516nm_a=13nm_b=20nm_m=2_n=2__vs__m=6_n=16.svg)
```
   filter_mrc -w 25.16 -dog 130 200 -exponents 2 2
   filter_mrc -w 25.16 -dog 130 200 -exponents 6 16
```
By default **m** *=* **n** *= 2*,
however this can be overridden using the optional "**-exponents m n**"
(For your convenience, you can plot these functions for various parameters
using the "draw_filter1D.py -dogg A B a b m n" python script which is
located in the same directory.
In this case, the "a" and "b" parameters should be specified in voxels, not physical distance.)
*Note:* that these *a*, *a_x*, *a_y*, *a_z*,
and *b*, *b_x*, *b_y*, and *b_z* parameters are a factor
of *√2 times larger* than the corresponding *σ* ("sigma") parameter
traditionally used in the
[Gaussian function](https://en.wikipedia.org/wiki/Normal_distribution).
*Note:* The calculation is
[fast](https://en.wikipedia.org/wiki/Separable_filter)
if you use the default exponent of 2.
*Changing the exponents will slow down the filter considerably.*




####  Implementation of LoG filters and blob-detectors:

To speed up the calculation,
the Difference-of-Gaussian approximation to the Laplacian-of-Gaussian
filter us used.
Specifically, the original image is convolved with a
[Difference-of-Gaussians (DoG)](https://en.wikipedia.org/wiki/Blob_detection#The_difference_of_Gaussians_approach)
filter, *h(x,y,z)*
```
   h(x,y,z) = scale * ( A*exp(-0.5*r_a^2) - B*exp(-0.5*r_b^2) )
  where r_a = √((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))
    and r_b = √((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))
      a_x = σ_x*(1-0.5*δ), a_y = σ_y*(1-0.5*δ), a_z = σ_z*(1-0.5*δ)
      b_x = σ_x*(1+0.5*δ), b_y = σ_y*(1+0.5*δ), b_z = σ_z*(1+0.5*δ)
    scale = (1.0 / δ^2)
```
The A and B parameters are determined automatically by normalization.
The "*δ*" parameter is *0.02* by default.
(This can be overridden using the "-dog-delta δ" argument.
A smaller "*δ*" value may improve the approximation,
but may also result in a noisier filtered image.)
As always, the width of the Gaussian (the *σ_x*, *σ_y*, *σ_z* arguments) should be specified in units of physical distance, not in voxels.
The *A* and *B* coefficients will be automatically chosen using normalization
(so that the discrete integral sums of *exp(-0.5*r_a^2)* and *exp(-0.5*r_b^2)*
functions are both 1).
The filter is multiplied by (1.0 / δ^2) to achieve *scale invariance*.
*("Scale invariance" means an object of width W filtered with a Gaussian
of width σ, receives the same score as an
object of width 2W filtered with a Gaussian of width *2σ*, for example.)
This also means that you can use the
"**-log**",  "**-log-r**", and "**-log-d**" arguments
to perform scale-free-blob-detection, *one Gaussian-width at a time,*
in such a way that that you can directly compare
the results of using different Gaussian-widths
(the same way they are compared when using
"**-blob**", "**-blobr**", and "**-blob**").

In all versions, the filter is truncated far away from the central peak at a point which is chosen automatically according the shape of the filter selected by the user.  However this can be customized using the "-truncate-threshold" and "-truncate" arguments if necessary (see below).

*Note:* The Gaussian "σ" arguments (*σ*, *σ_x*, *σ_y*, and *σ_z*)
and the *t* ("time") scaling parameter (traditionally used
in the scale-space literature), are related according to:
```
   σ^2 = t
```
