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






### -auto-thresh score -supervised file_accept.txt file_reject.txt

When discarding poor-scoring blobs using the
["-discard-blobs"](./doc_filter_mrc.md#-discard-blobs)
argument, it can be difficult to choose the correct threshold
for differentiating between real blobs and noise.
As an alternative to specifying the threshold(s) manually,
(using the "-minima-threshold" or "-maxima-threshold" arguments),
you can instead supply examples of blobs that you want to keep,
and blobs you want to discard 
(using the "-auto-thresh" and "-supervised" arguments).

*(Note: To use these arguments, the image must be associated 
 with a list of blobs that was detected by previously running 
 this program with the "-blob" argument on that image.
 You must specify that list of blobs using the "-discard-blobs" 
 argument, which must also appear in the argument list,
 along with "-auto-thresh" and "-supervised".)*


Blobs will be discarded if their "score" is not sufficiently high.
(Equivalently, they will be discarded if they are not sufficiently
 bright or dark compared to their surroundings.)
These blob locations are listed
on separate lines in the "file_accept.txt" and "file_reject.txt" arguments.
The file format for both of these files is described
[here.](#-must-link-FILE)
*(Note: As explained in that link, you can obtain the contents of this file
        by clicking on objects of interest in IMOD.)*
This range of allowed scores is determined automatically from these examples.
The (upper-bound or lower-bound) threshold will be chosen
which minimize the number of incorrectly classified blobs.
Here is an example of the output:
```
  threshold lower bound: -inf
  threshold upper bound: -37.921
```
Equal weight is given to false-positives and false-negatives,
(so choose your examples accordingly).
Choosing blobs which are "edge-cases" is recommended.
(IE. blobs that would difficult to classify or are barely visible.)


*(Note: These arguments must be supplied together as shown.
The "-auto-thresh" argument will have no effect unless
you are also using the "-supervised" and the
"-discard-blobs" arguments as well.
It is also common to include the "-radial-separation"
or "-max-volume-overlap" arguments as well, but these arguments are optional.
The "score" argument is not a number, but literally the word "score".
Eventually, it will be possible to supply a list of other criteria used
for classifying blobs, but as of 2019-5-07, only "score" is available.)*

*(Note: You must provide examples of both
blobs that you want to keep and
blobs that you want to discard.
These example blobs must already be present
within the list of blobs that you have provided
to the "-discard-blobs" argument.
If you are using a mask, then blobs whose centers lie outside
the mask will be ignored.)*

*(Note: Usually one of these bounds is infinite, and the other is finite.
The infinite bound can be be ignored.
Occasionally an upper bound AND a lower bound are both not infiinte.
This usually happens when you did not supply enough training data.
Be sure to include some positive training examples of obvious blobs with good
(large magnitude) scores, as well as
negative training examples with faint blobs and poor scores.
If you have much more positive data than negative data (or visa-versa),
then include more data so that the number of positive and negative
training examples is not so lopsided.
Alternatively, you can ignore one of these thresholds
(whichever threshold which does not make sense).  For example, if your blobs
are light on a dark background, then only the lower bound is relevant.
You can ignore the upper bound in that case. Similarly, if your blobs are dark
on a white background then you can ignore the lower bound.


### -auto-thresh score -supervised-multi list_of_files.txt

*(WARNING: This is an experimental feature as of 2019-6-19.
           It has not been tested carefully, and may be difficult to use.)*

The "*-supervised-multi*" argument is a variant of "*-supervised*"
which allows you to use multiple training sets 
associated with multiple different images.

*(Note: Each of those images must be associated with a list of blobs 
  that were detected by previously running this program
  with the "-blob" argument on that image.)*

The program will consider the blobs at these locations, 
and calculate a range of scores which encloses
the positive training data with as few errors as possible.
The result is printed to the standard-error.
Here is an example of the output:
```
  threshold lower bound: -inf
  threshold upper bound: -131.016
```

*(Note: Unlike the "-supervised" argument, 
you cannot use the "-supervised-multi" argument simultaneously
with the "-discard-blobs" or "-mask" arguments.
Discarding overlapping or masked blobs must be done in advance.
See below.)*

#### File format details:

The file supplied to this argument (eg, "list_of_files.txt"),
should be a 3-column text file containing file names.
Each line of this file contains the names of 3 other files containing 
information about a particular image.

Suppose you have collected training data from 5 different images.
Example ("list_of_files.txt"):

```
training_pos_1.txt  training_crds_neg_1.txt  blob_info_1.txt
training_pos_2.txt  training_crds_neg_2.txt  blob_info_2.txt
training_pos_3.txt  training_crds_neg_3.txt  blob_info_3.txt
training_pos_4.txt  training_crds_neg_4.txt  blob_info_4.txt
training_pos_5.txt  training_crds_neg_5.txt  blob_info_5.txt
```
The first two files, *("training_pos_1.txt", "training_crds_neg_1.txt")*,
contain positive and negative training data corresponding to objects of
interest from that image (image#1).
(The format of these files is explained earlier in the documentation 
for the "-supervised" argument.)
The third file listed on that line contains a list of blobs that were 
detected in the image.)
(Those files are generated by running this program on that image
 separately using the ["**-blob**"](#Blob-detection) argument.
 See the docmentation for the ["**-blob**"](#Blob-detection) argument
 for details.)

*NOTE: This list of blobs (eg, contained in "blob_info_1.txt")
must not contain overlapping blobs, or blobs that lie
outside the mask.*
So before you use the "-supervised-multi" argument,
you should delete the overlapping blobs in each of these lists
by running this program again separately on the 
blobs you detected earlier
(by using the
["-discard-blobs"](#-discard-blobs),
["-radial-separation"](#Automatic-disposal-of-overlapping-blobs),
or
["-max-volume-overlap"](#Automatic-disposal-of-overlapping-blobs),
and/or
["-max-volume-overlap-small"](#Automatic-disposal-of-overlapping-blobs)
and
["-mask"](#-mask-MRC_FILE)
arguments).
***All*** of this must be done separately for each of the N images you want to 
simultaneously analyze (5 in this example), ***before*** you attempt to use
this program with the "-supervised-multi" argument.

