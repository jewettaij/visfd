filter_mrc
===========

**filter_mrc** applies a filter to a tomogram in the X,Y,Z directions
and saves the result as a new .mrc/.rec file.
This program can be used to rescale or invert a 3-D image,
remove high or low frequencies (smoothing, edge detection, band pass filter).
perform 3-D blob detection.
Currently, the program supports the following list of filters:
(generalized) Gaussians,
(generalized) Difference-of-Gaussians,
Laplacian-of-Gaussians, and others.
Both isotropic and anisotropic filters are supported.


This program supports masks and thresholding using
the "*-mask*" and "*-thresh*" (and "*-thresh2*" and "*-thresh4*") arguments.


## Usage Example:

```
filter_mrc \
   -in Bdellovibrio_1K.rec \
   -out Bdellovibrio_1K_ribosome_peaks.mrc \
   -w 19.2 \
   -template-gauss 100 130 \
   -exponents 4 8 \
   -cutoff 0.0001 \
   -mask Bdellovibrio_1K_mask_water=0_periplasm=1_cytoplasm=2.mrc \
   -mask-select 2 -mask-out 0.0
```

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
(This includes all of the width parameters used with the "-gauss", "-blob", "-dog", and "-template-gauss" filters.)

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


The "**-thresh**", "**-thresh2**", "**-thresh4**",
filters are used to clip and rescale voxel intensities,
or to select voxels whose intensities lie within in a certain range.


The "**-gauss**" filter uses a simple (low-pass)
[Gaussian blur](https://en.wikipedia.org/wiki/Gaussian_blur)
filter.
(The "**-gauss-aniso**" filter allows you to apply different amounts of blurring to the X,Y, and Z directions.)
[Generalized Gaussians](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1) are supported.


The "**-blob**", "**-blobd**", and "**-blobr**" filters can be used for [scale-free blob detection](https://en.wikipedia.org/wiki/Blob_detection#The_Laplacian_of_Gaussian).  Objects in the image of various sizes can be emphasized or selected using this filter.  Unfortunately, the computation is slow compared to filters like -gauss and -dog.  (However, if you can estimate the size of the objects you want to detect in advance, then you can use "**-blob1**", "**-blobr1**", or "**-blobd1**" which are much faster.)

The "**-dog**"" filter uses a (band-pass)
[Difference-Of-Gaussians](https://en.wikipedia.org/wiki/Difference_of_Gaussians)
filter which can be used as a frequency band-pass filter (for both high and low frequency removal).
(The "**-dog-aniso**" filter allows you to modify the properties of the filter in the X,Y,Z directions.  Filter sharpness can be customized using the "-exponents" argument.)
[Generalized Gaussians](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1) are supported.

The "**-template-gauss**" filter
is also used for blob detection.
It performs
[template-matching](https://en.wikipedia.org/wiki/Template_matching)
on 
[(generalized) Gaussians](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1)
Specifically, it calculates the overlap 
(cross-correlation) 
of the image with the template
(a generalizec Gaussian),
as well as the RMSE (root-mean-squared-error) between the
original image and the Gaussian template
after optimal overlap.
(IE., after optimal scaling of the template voxel intensities and background subtraction).
This allows the user to insure that the blobs in the image
actually resemble the shape of the template (in this case a Gaussian),
and cannot be easily fooled by a particularly bright or dark voxel
in the source image.  (See below for details.)




## Filter Details:

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

The filter is truncated far away from the central peak at a point which is chosen automatically according the σ, σ_x, σ_y, σ_z parameters selected by the user.  However this can be customized using the "-cutoff" and "-window-ratio" arguments if necessary (see below).


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

Alternatively, if the **-dogg** or **-dogg-aniso* arguments are selected, 
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
The filter is truncated far away from the central peak at a point which is chosen automatically according the shape of the filter selected by the user.  However this can be customized using the "**-cutoff**" and "**-window-ratio**" arguments if necessary (see below).
![A comparison of Difference-of-Gaussians Difference-of-Generalized-Gaussians filter weights](./images/example_dogxy_w=2.516nm_a=13nm_b=20nm_m=2_n=2__vs__m=6_n=16.svg)
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


### -doggxy
The **-doggxy** argument must be followed by 3 numbers:
```
  -doggxy  a  b  c
```
When the "**-dogxy**" filter is selected,
the original image is convolved with the following function:
```
   h(x,y,z) = h_xy(x,y) * h_z(z)
```
 In the XY plane, the filter used is:
```
   h_xy(x,y) = A*exp(-(|r|/a)^m) - B*exp(-(|r|/b)^n)
           r = √(x^2 + y^2)
 and A,B coefficients are chosen so that each Gaussian is normalized

```
 The "m" and "n" parameters are exponents.
 They are both set to 2 by default, however you can override them
 using the "**-exponents m n**" argument.

 Along the Z direction, the filter used is a simple Gaussian:
```
   h_z(r) = C * exp(-0.5*(z/c)^2)
```
(The "C" constant is determined by normalization.)

The computational cost of "**-gdogxy**" lies in between the ordinary and *generalized* difference-of-Gaussian (DOG) filters discussed above.  (Features in electron tomography are typically blurred more in the Z direction due to the effect of the missing wedge.  So it may be pointless and impossible to use the computationally more expensive generalized ("**-gdog**", "**-exponents**") filter in an effort to find the precise boundaries of objects in the Z direction.  In these cases, the "**-gdogxy**" (and "**-dog**") filters may work just as well and are significantly faster.)


### -blob,  -blobr,  -blobd,  -blob1,  -blobr1,  -blobd1,  -blob1-aniso

When the "**-blob1**", "**-blobr1**", and "**-blobd1**" arguments are selected a
[Laplacian-of-a-Gaussian (LOG)](https://en.wikipedia.org/wiki/Blob_detection)
filter is applied on the source image.
(See implementation details below.)
The Laplacian-of-a-Gaussian filter applies a Gaussian blur to the image
(whose width is specified by the user), 
and then calculates the Laplacian of the result.
Features in the image of various sizes can be emphasized 
using this kind of filter.
Detected blobs can be displayed to the user using the 
"**-o** filename.mrc" (and "**-spheres**") arguments.

The "**-blob1**", "**-blob1-aniso**", "**-blobr1**", and "**-blobd1**"
arguments is typically chosen when the user already knows the approximate
size of the objects or features they want to detect in the image.
When these arguments are selected, an LOG filter will be applied
to the image only once using a Gaussian width (σ) specified by the user,
and the resulting filtered image will be saved as a file.
(The local minima and maxima of this image can then be extracted later using
 the "**-find-minima**", "**-find-maxima**" and "**-max-overlap**" arguments.)
```
  -blob1  σ
```
There is also an anisotropic version of this filter.
```
  -blob1-aniso  σ_x  σ_y  σ_z
```

####  -blob,  -blobr,  -blobd

The "**-blob**", "**-blobr**", and "**-blobd**" arguments are used for
[Scale-Free Blob-Detection](https://en.wikipedia.org/wiki/Blob_detection).
When this is used, the program will apply a LOG filter to the image 
multiple times using a range of Gaussian widths (σ) (specified by the user)
in an attempt to determine the correct size (scale) for the relevant objects
in the image automatically.  (As such, this operation is comparatively slow.)
A list of local minima and maxima in *X,Y,Z* space (and scale-space)
will generated and saved in a file, using the method described in:
Lindeberg,T., Int. J. Comput. Vision., 30(2):77-116, (1998)

The "**-blob**", "**-blobr**" and "**-blobd**" arguments
filters are followed by 4 arguments:
```
  -blob   file_name.txt  σ_min  σ_max  gratio
  -blobr  file_name.txt  r_min  r_max  gratio
  -blobd  file_name.txt  d_min  d_max  gratio
```
A LOG filter will be applied to the image using different Gaussians
whose widths (σ) vary between "**σ_min**", and "**σ_max**"
(or, equivalently whose radii vary between "**r_min**" and "**r_max**").
Each Gaussian will be wider than the previous Gaussian by a fraction
(no larger than) "**gratio**", a number which should be > 1.
(**1.01** is safe, recommended choice,
 but you can speed up the calculation by increasing this parameter.)
"**file_name.txt**" is the name of a file which will store the 
locations of all the blobs detected in the image
as well as their sizes and scores (see below).
These detected blobs are either local minima or maxima in
X,Y,Z,[σ space](https://en.wikipedia.org/wiki/Scale_space_representation).
By default, two files will be created, named
*file_name_minima.txt* (for storing local minima), and 
*file_name_maxima.txt* (for storing local maxima).
Each file is a 5-column ascii text file with the following format:
```
x1 y1 z1 size1 score1
x2 y2 z2 size2 score2
 :  :  :   :     :
xM yM zM sizeM scoreM
```
...where
"**M**" is the number of blobs (maxima or minima) which were detected,
On each line of that file,
**x,y,z** are the coordinates for the blob's center,
**size** is the size of the blob, 
(characterized either using either σ, radius, or diameter, as explained below),
and **score** is the intensity of that voxel after 
a LOG filter of that size was applied to it.
The list is ordered from high score to low score
(for maxima, or low score to high score for minima).

#### Automatic disposal of poor scoring blobs

There are an enormous number of local minima in X,Y,Z,scale space.
By default, local maxima whose score is less than 50% of the global maxima
are discarded.  A similar rule is used for discarding local minima.
You can override these rules by supplying the following arguments:
```
   -blob-minima-ratio  ratio_min
   -blob-maxima-ratio  ratio_max
```
...where "ratio_min" and "ratio_max" are fractions between 0 and 1.
If you prefer, you can alternately specify the exact threshold used
to discard local minima and maxima (instead of specifying ratios).
You can do this using the following arguments:
```
   -blob-minima-threshold  thresh_min
   -blob-maxima-threshold  thresh_max
```
...where "thresh_min" and "thresh_max" are the thresholds below which,
and above which, the blob score must lie in order not to be discarded,
respectively.
*(Note: You can either use thresolds or ratios, but you cannot mix both.)*


#### Automatic disposal of overlapping blobs (non-max suppression)

Whenever two blobs overlap, the one with the better score
(ie. higher score for maxima, and lower score for minima)
is superimposed upon the one with a poorer score.
By default, if the distance between two blobs is less than 1.0 times
the sum of their radii, the poorer-scoring blob will be discarded.
This can be overridden using the "**-blobr-sepration fraction**" argument
where "**fraction**" is typically a number between 0 and 1.
(A larger number means less overlap is permitted.
 Setting *fraction=0.0* disables this feature.
 You can also specify the overlap in terms of the Gaussian width of each blob, 
 σ, using the "**-blob-sepration fraction**" argument.  In that case,
 setting *fraction* to 1 will discard blobs if they lies closer than 
 the sum of their Gaussian widths, σ1+σ2 = (r1+r2)/√3.
 *Be careful not to use "-blob-sepration" when you really mean
 "-blobr-sepration"*)


#### Specifying the radius or diameter of the objects of interest:

The "**-blobr**", "**-blobd**", (and "**-blobr1**" "**-blobd1**") arguments
are all variants of the
"**-blob**" (and "**-blob1**") argument whose parameters are (more conveniently)
specified by the approximate radius(r≈σ√3) or diameter(d≈σ2√3) of the objects 
that you wish to detect within the image (instead of the Gaussian width, σ).


####  Visualizing blobs using the "**-o**" argument

*If* the "**-o**" argument is specified during blob-detection,
then a new tomogram will be created with each blob represented
by a hollow spherical shell that surrounds the point of 
interest, whose *radius* and *brightness* indicate *size* and *score* of 
that blob, respectively.
*(The intensity of voxels belonging to the spherical
shell should exactly match the score of that blob.  
In IMOD, blobs with good scores typically appear as black or white spheres,
and grey spheres correspond to blobs with poor scores.
The score of a given blob
can be queried in IMOD by clicking on a voxel somewhere on the spherical shell
and selecting the "Edit"->"Point"->"Value" menu option.)

The radii and brightnesses of the spheres can be overidden using the
"**-sphere-radius r**" 
and "**-sphere-foreground intensity**" arguments, respectively.
(Note: The *r* parameter should be in physical units/Angstroms, not voxels.
 If you prefer to scale all the radii up or down by the same ratio,
 use the "**-sphere-scale ratio**" argument instead.)

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
 The default scale ratio is 0.3.)


####  Visualizing blobs using the "**-spheres**" argument

If you did not generate an image showing the blob locations at the time
you performed the blob detection, you can generate this image later using:
```
  filter_mrc -spheres blob_filename.txt -o new_image.mrc
```
This can be useful because, quite often, 
these default score thresholds too permissive.
(Far too many blobs are detected.)
As before, 
the user can discard poor scoring blobs by using any one of these arguments
*(explained elsewhere)*:
"**-blob-minima-threshold**", 
"**-blob-maxima-threshold**", 
"**-blob-minima-ratio**", 
"**-blob-maxima-ratio**"

#### Recommendation:

Blob-detection is computationally expensive, 
but it only has to be performed once.
Typically, users will perform blob detection with the default (permissive)
score thresholds, so that they keep most of the minima and maxima.
Then later, they will run **filter_mrc** with the 
"**-spheres**" and "**-blob-minima-threshold**" (or "**-blob-maxima-threshold")
arguments several times with different thresholds to decide which of 
these blobs are worth keeping.





####  Implementation:

To speed up the calculation, 
the Difference-of-Gaussian approximation to the Laplacian-of-Gaussian 
filter us used.
Specifically, the original image is convolved with a 
[Difference-of-Gaussians (DOG)](https://en.wikipedia.org/wiki/Blob_detection#The_difference_of_Gaussians_approach)
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
A smaller "*δ*" value may improve the approximation.)
The width of the Gaussian (the *σ_x*, *σ_y*, *σ_z* arguments) should be specified in units of physical distance, not in voxels.
The *A* and *B* coefficients will be automatically chosen using normalization
(so that the discrete integral sums of *exp(-0.5*r_a^2)* and *exp(-0.5*r_b^2)*
functions are both 1).
The filter is multiplied by (1.0 / δ^2) to achieve *scale invariance*.
*("Scale invariance" means an object of width W filtered with a Gaussian 
of width σ, receives the same score as an
object of width 2W filtered with a Gaussian of width *2σ*, for example.)
This also means that you can use the
"**-blob1**",  "**-blobr1**", and "**-blobd1**" arguments 
to perform scale-free-blob-detection, *one Gaussian-width at a time,*
in such a way that that you can directly compare
the results of using different Gaussian-widths 
(the same way they are compared when using
"**-blob**", "**-blobr**", and "**-blob**").

In all versions, the filter is truncated far away from the central peak at a point which is chosen automatically according the shape of the filter selected by the user.  However this can be customized using the "-cutoff" and "-window-ratio" arguments if necessary (see below).

*Note:* The Gaussian "σ" arguments (*σ*, *σ_x*, *σ_y*, and *σ_z*)
are the related to the *t* ("time") scaling parameter traditionally used 
in the scale-space literature according to:
```
   σ^2 = t
```
(See implementation details above)


*Note:* As a reminder, 
a more general version of the the DOG filter used above can be 
applied to the image by using the "-dog" and "-dogg" arguments,
(although it's not clear whether this would be useful for blob detection).
These alternate versions of the DOG filter allow the user 
to specify the *a_x, a_y, a_z, b_x, b_y, b_z* parameters  directly,
and use
[generalized Gaussians](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1)
instead of
[ordinary Gaussians](https://en.wikipedia.org/wiki/Normal_distribution)
(by using the "-exponents m n" argument).




### -template-gauss
The **-template-gauss** and **-template-gauss-aniso** arguments must be followed by one or more numbers specifying the width of the Gaussian filter to apply to your image:
```
   -template-gauss-aniso  a_x  a_y  a_z  b_x  b_y  b_z
```
or:
```
   -template-gauss  a  b   (this means a = a_x = a_y = a_z, b = b_x = b_y = b_z)
```
When the "-template-gauss" or "-template-gauss-aniso" filter is selected, the
original image compared with
[a generalized Gaussian function](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1):
```
   h(x,y,z) = A*exp(-r^m)
    where r = √((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))
```
The width of the Gaussian (ie, the *a_x*, *a_y*, *a_z* parameters)
should be specified in units of physical distance, not in voxels.
(The A coefficient will be determined automatically by normalization.)

When using this filter,
the *relative* intensity of the voxels (compared to voxels nearby)
is compared with the (Gaussian) template.
The *relative* voxel intensity is defined as the
original voxel intensity minus the average intensity of nearby voxels.
The *average* voxel intensity of
nearby voxels is calculated (by convolvution)
using the following weights:
```
   w(x,y,z) = B*exp(-r^n)
    where r = √((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))
```
The *b_x*, *b_y*, *b_z* arguments define an ellipse over which this average
is calculated. (They should be specified in physical units, not voxels.
The **n** exponent defines the sharpness of the ellipsoidal boundary.)
(The B coefficient is determined automatically by normalization.)

The "**m**" and "**n**" parameters are both *2* by default,
however you can override them using the "**-exponents m n**" argument.
*(Doing so will slow down the calculation considerably.)*

#### Two different ways to assess goodness of fit:
After the *relative* voxel intensity is calculated
(original - average_of_nearby),
the (similarly weighted) average of the *template* image
is also subtracted from the template image.
Finally, the intensities of the shifted template are multiplied by a constant
("*c*") so that the template intensities optimally overlap with the
*relative* voxel intensities in the original image.
(I.E., with lowest possible root-mean-squared-error, RMSE).
This fitting is performed everywhere, sliding the template function
(a Gaussian whose peak height=1) along the image, centered it on every voxel
and comparing the fit with the voxels at that location in the original image.
The value of the constant, **c**, is calculated everywhere, and the new
image file (specified by the "*-o*" argument)
will store the value of **c** calculated at all of these locations.
The **c** value turns out to equal the 
["cross-correlation"](https://en.wikipedia.org/wiki/Template_matching#Template-based_matching_explained_using_cross_correlation_or_sum_of_absolute_differences),
between the relative voxel intensities of the source image and the template.
The cross-correlation frequently used in template matching,
however it is often an inadequate measure of similarity on its own.
A high **c** value indicates that the Gaussian that fits as well as possible
must have a high peak, but it does not indicate whether or not 
that Gaussian is actually a good fit.
A single bright or dark voxel could cause a large **c** value.



Consequently
the RMSE (root-mean-squared-error) of comparison is also calculated
everywhere (with the template Gaussian centered on every voxel).
The RMSE is the root sum-squared error between the source image
and the template (after optimal intensity shifting and scaling).
These RMSE values are saved in a new image whose filename
ends in "rmse.rec".
*(The user can decide to exclude voxels where either
the c value is too low, or
the RMSE value is too high,
by combining the two
files together using the "combine_mrc" program
which is documented [here](doc_combine_mrc.md).)*


(Details: The filter is truncated far away from the central peak at a point which is chosen automatically according the a_x, a_y, a_z, b_x, b_y, b_z parameters selected by the user.  However this can be customized using the "-cutoff" and "-window-ratio" arguments if necessary.)




## Thresholding, Clipping and Rescaling (Intensity Map Filters):

Thresholding filters provide a way to rescale (map) the
intensity of all of the voxels in the image to new intensities.
For example, you could multiply these intensities by a constant,
and/or clip them within a specified range.
Thresholding operations can be used to insure that the intensities
(ie. brightnesses) for every voxel in the image lie between 0 and 1.
They can also be used to select only voxels whose intensities lie within
a narrow range of interest.

*Note:* All threshold operations are performed *after* normal filtering operations have been applied (such as Gaussian blurs, -gdog filters, or -dog filters).

### -invert

This filter replaces bright voxels with dark ones, and visa versa,
while keeping the average voxel brightness the same.
(To calculate the new voxel intensity,
 the *difference* between the original voxel intensity
 and the average intensity is *subtracted* from the average intensity.)
Inverting an image can be useful to prepare files which you plan
to read with other programs (such as ChimeraX and IMOD).

### -rescale


If you use the **-rescale** argument, then
image voxel intensities in the final image
(after all other processing has been applied)
will be shifted and rescaled so that the
minimum and maximum voxel intensities are 0 and 1.

*(Note: As of 2018-1-31, IMOD reports voxel brightnesses between 0 and 255
  even if the actual voxel brightnesses are floating point numbers between 
  0 and 1.  Keep this in mind when using that software.
  Instead, you can use the "histogram_mrc.py" script to 
  see if your intensity values lie in the range you expect.)*

### -clip a b
The "**--clip a b**" argument will
restrict the range of voxel intensities to a range of your choosing.
This can be useful to eliminate extremely dark or bright voxels from
your image.
(Such spurious voxels can make template matching less reliable.)

```
-clip 0.48 0.52
```
The resulting voxels will have intensities in the range from 0.48 to 0.52.
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
*Note: You can use the "**-mask**" argument to exclude certain voxels 
       from the clipping and thresholding operations.
       When using "**-mask**" together with "**-cl**", then 
       these voxels will also be excluded from consideration when 
       calculating the average and spread(σ) of voxel intensities.


### -thresh  threshold

If the "-thresh" argument is passed, it must be followed by a number ("thresh01").  Voxels with intensities *below* this number will be replaced 0, and voxels *above* this number will be replaced with 1.  For example:
```
-thresh 0.5
```
results in:
```
 output
 intensity
  /|\                      _____________________________\
 1 |                      |                             /
   |                      |                 
   |                      |            
 0 |______________________|  __________________________\ input
                          ^                            / intensity
                         0.5
                     (threshold)

```


*Note: When choosing thresholds, the 
       [histogram_mrc.py program](doc_histogram_mrc.md)
       can be useful.
       It displays the range of voxel intensities in an image.*


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
    /|\                              _________________
 1   |                           _.-'                 
     |                       _,-'                 
     |                   _,-'            
 0   |________________,-'                     ________\ input
                      ^              ^                / intensity
                     0.48          0.52
                  (thresh_a)    (thresh_b)
```
***Alternatively***, if the threshold parameters are reversed
(if ***thresh_b < thresh_a***), as in this example:
```
-thresh2 0.52 0.48
```
...then the output intensity will be inverted
(i.e., bright voxels become dark, and dark voxels become bright):

```
   output
  intensity
     /|\
      |
  1   |________________
      |                `-._
      |                    `-._
      |                        '-._         
  0   |                            `-.___________________\ input
                       ^              ^                  / intensity
                      0.48          0.52
                   (thresh_b)     (thresh_a)
```

*Note: When choosing thresholds, the 
       [histogram_mrc.py program](doc_histogram_mrc.md)
       can be useful.
       It displays the range of voxel intensities in an image.*


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
 1 |                 ________________                
   |             _,-'                `-._
   |         _,-'                        `-._
 0 |______,-'                                `-._______\ input
         ^         ^                 ^         ^       / intensity
        0.3       0.4               0.5       0.6
   thresh_01_a  thresh_01_b    thresh_10_a  thresh_10_b
```
This assumes the numbers were listed in *increasing* order.
***Alternatively***, if the parameters are listed in *decreasing* order,
for example:
```
 -thresh4 0.6 0.5 0.4 0.3
```
...then the output is inverted:
```
 output
 intensity
  /|\                                                   
 1 |_____                                       _________
   |     `-._                               _.-'       
   |         `-._                       _,-'             
 0 |             `-._________________,-'          _______\ input
         ^         ^                 ^         ^         / intensity
        0.3       0.4               0.5       0.6
   thresh_10_b  thresh_10_a    thresh_10_b  thresh_10_a
```


*Note: When choosing thresholds, the 
       [histogram_mrc.py program](doc_histogram_mrc.md)
       can be useful.
       It displays the range of voxel intensities in an image.*


## Masking

The optional "-mask", "-mask-select", and "-mask-out" arguments allow you to
ignore certain voxels from the source image (tomogram).

Using "masks" allows the user to perform the filtering considering only a
subset of voxels in the tomogram (ie, the "mask").

```
   -mask  file.mrc
```
The argument following the
"-mask" argument should be the name of a tomogram file (MRC format) of the
same size as the input tomogram.
During the filtering process, only voxels belonging to the tomogram
with non-zero voxel values in the "mask" tomogram will be considered,
and these remaining voxels will be multiplied by the mask value there
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
If you are using "-rescale", then
"-mask-out 1" assigns them to the brightest voxels in the image.)


```
   -mask-select  intensity_value
```
If the "-mask-select" argument is specified, then instead of considering all
voxels with non-zero values from the "mask" image/tomogram,
only voxels whose mask intensity equals the number following
this argument will belong to the mask.
All other voxels will be ignored during filtering.
*(Note: This disables "soft" masking.  In other words,
  during the process of filtering, all selected voxels will be weighted equally
  during filtering, and all others will be completely ignored.)*




## Additional Arguments:

### Filter Size
```
   -cutoff threshold
```
This specifies the number of voxels in the filter which will be convolved
with the image.
The filter window width is extended in the X, Y, and Z directions until
the filter intensity decays to less than **threshold** (a fraction),
compared to its peak value.
(For Difference-of-Gaussian filters, both Gaussians must decay to this value.)
If unspecified, this parameter is set to 0.02.
(This overrides the "-window-ratio" argument.)

```
   -window-ratio Wx Wy Wz
```
If you prefer to specify the size of the filter window manually, you can
use the "-window-ratio" argument.
This argument specifies the size of the 3-D filter box in x,y, directions
(in units of the **a** and **b** distance parameters, or **σ** ("sigma")
 parameters which appear in the formulas above).
This overrides the "-cutoff threshold" argument.

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
