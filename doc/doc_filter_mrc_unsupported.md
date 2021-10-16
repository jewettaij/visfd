Unsupported features:
===========
* DEPRECIATION WARNING:
This file contains documentation for features of the "filter_mrc" program
which are no longer maintained or supported (and were probably never useful).
This file may be deleted in the future.*

When using "filter_mrc", the user must select the kind of filter that will 
be used by supplying a list of command-line arguments to the program.
*The complete list of arguments and filters for the "filter_mrc" program
 is described [here](./doc_filter_mrc.md).*
The small list of optional arguments shown below correspond to 
a filters that will probably be deleted in the future.

### -doggxy
The **-doggxy** argument must be followed by 3 numbers:
```
  -doggxy  a  b  c
```
When the "**-doggxy**" filter is selected,
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

The computational cost of "**-doggxy**" lies in between the ordinary and *generalized* difference-of-Gaussian (DOG) filters discussed above.  (Features in electron tomography are typically blurred more in the Z direction due to the effect of the missing wedge.  So it may be pointless and impossible to use the computationally more expensive generalized ("**-dogg**", "**-exponents**") filter in an effort to find the precise boundaries of objects in the Z direction.  In these cases, the "**-doggxy**" filter may work just as well and is faster.)



### -template-gauss
The "**-template-gauss**" filter can also be used for blob detection.
It performs
[template-matching](https://en.wikipedia.org/wiki/Template_matching)
on
[(generalized) Gaussians](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1)
Specifically, it calculates the overlap
(cross-correlation)
of the image with the template
(a generalized Gaussian),
as well as the RMSE (root-mean-squared-error) between the
original image and the Gaussian template
after optimal overlap.
(IE., after optimal scaling of the template voxel intensities and background subtraction).
This allows the user to insure that the blobs in the image
actually resemble the shape of the template (in this case a Gaussian),
and cannot be easily fooled by a particularly bright or dark voxel
in the source image.

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
image file (specified by the "*-out*" argument)
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




### -blob-intensity-vs-radius center_type blobs_file.txt base_name

The "**-blob-intensity-vs-radius**" argument creates a list of txt
files (one per blob) storing the brightness as a function of
distance from the center of the blob.
The distance is expressed in *physical* units
(such as Angstroms, **not voxels**).
The *first* argument **center_type** must be one of these choices:
***center***, ***max***, ***min***.
This indicates whether we should measure distance to the *center*
of the blob, or the location of its intensity *maxima* or *minima*,
respectively.
The *second argument* is the file containing a list of blob
locations and radii.
It should be in the same format as the file created by the
"**-blob**" argument  (See documentation above.)
The *third argument* is the name we will attach to the beginning
of each file that we create.  The program will create a series
of files with names lie "base_name_1.txt", "base_name_2.txt", ...
"base_name_N.txt", where *N* is the number of blobs in the
"blobs_file.txt" file.
