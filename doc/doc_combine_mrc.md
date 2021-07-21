
combine_mrc
===========
**combine_mrc** is a program for combining two volumetric images (i.e. tomograms, both of identical size) into one image/tomogram, using a combination of addition, multiplication, and thresholding operations.  These features can be used perform binary operations between two images which are similar to "**and**" and "**or**" operations.  ("**not**" operations are also possible.  See below.) As with the "filter_mrc" program, you can also use the "-mask" argument to restrict the operation to certain voxels from the image.
*(A detailed description of the threshold and rescaling functions used in this step are provided at the end of the the "doc_filter_mrc.md" file.)*


### combine_mrc examples:

#### + binary operation
To **add** the brightness values between to tomograms and save the result in a new file ("out_file.mrc"), use:
```
   combine_mrc file1.mrc  +  file2.mrc  out_file.mrc
```
#### * binary operation
To **multiply** the brightness values use:
```
   combine_mrc file1.mrc "*" file2.mrc  out_file.mrc
```
(In the last example, the quotes around the "\*" character are necessary to prevent the shell interpreting "\*" incorrectly as a wildcard.)

For completeness, division and subtraction are also supported:
```
   combine_mrc file1.mrc - file2.mrc  out_file.mrc
   combine_mrc file1.mrc / file2.mrc  out_file.mrc
```
The following examples apply 1, 2 or 4 thresholds to the input tomograms before performing the **+** or **\*** operation.  You specify the thresholds you want to use with commas placed after the file name (no spaces).  The following command will replace all of the brightness of all the voxels whos brightess is below 0.5 with 0, and all of brightnesses above 0.5 with 1, and *then* multiply them together:
```
   combine_mrc file1.mrc,0.5 "*" file2.mrc,0.5 out_file.mrc
```
(Note that 8-bit and 16-bit integer brightnesses are replaced with floating point numbers in the range from 0 to 1 beforehand, and the resulting tomogram is saved in 32bit float format.)

### -rescale

If you use the **-rescale** argument, then
image voxel intensities will be shifted and rescaled so that the
minimum and maximum voxel intensities are 0 and 1.
This will be done for both the input images and well as the output image.
*(Note: Unlike the thresholding operations described below,
 this does not erase (clip) any voxel intensities.)*


### Applying threshold filters to the (input and output) images:

It is often convenient to rescale or clip the brightnesses of the voxels from either image (tomogram) *before* adding or multiplying them together.  This way, you can insure that most of the voxels in the image are either ***0*** or ***1*** beforehand.  (This way, multiplying and adding the resulting voxel brightnesses is equivalent to performing an "and" and "or" gate operations on these 0,1 values.)  The resulting image (tomgram) created by *combine_mrc* can be useful for segmentation (or masking).

The following command will replace all of the brightness of all the voxels whose brightess falls below *0.48* with ***0***, and all of brightnesses above *0.52* with ***1*** (linearly scaling any voxels with brightnesses between *0.48* and *0.52* to fill the range from *0* to *1*).  *Then* it will add them together
```
   combine_mrc file1.mrc,0.48,0.52 + file2.mrc,0.48,0.52 out_file.mrc
```
Graphically, the relationship between the voxel's input intensity ("density")
and its output intensity is the following:
```
  output
  intensity
    /|\                              _________________
   1 |                           _.-'                 
     |                       _,-'                 
     |                   _,-'            
   0 |________________,-'                     ________\ input
                      ^             ^                 / intensity
                     0.48          0.52
                  (thresh_a)    (thresh_b)
```
*Note:* If the user also includes the "-rescale" argument,
then the resulting intensities will be vertically scaled between 0 and 1.
*Note:* If the order of thresholds is reversed, the inverse images is generated.
For example:
```
   combine_mrc file1.mrc,0.52,0.48 + file2.mrc,0.52,0.48 out_file.mrc
```
...runs the input image intensities through an inverted threshold filter
(i.e., bright voxels become dark, and dark voxels become bright):
```
 output
 intensity
  /|\
   |________________
 1 |                `-._
   |                    `-._
   |                        '-._            
 0 |                            `-.___________________\ input
                    ^              ^                  / intensity
                   0.48          0.52
                (thresh_b)     (thresh_a)
```
*(Again, a detailed description of the threshold and rescaling functions used in this step are provided at the end of the the "doc_filter_mrc.md" file.)*


### NOT-gates
This can be used to invert the image intensities so that bright voxels
are replaced by 0, and dark voxels are replaced by 1.

You can also apply thresholding to the *output image*, as demonstrated
[below](#NAND-gate example:)


#### OR-gate example:

For example, the following command will clip the output intensities
between 0 and 1.
```
   combine_mrc file1.mrc,0,1 + file2.mrc,0,1 out_file.mrc,0,1
```
This operation is similar to an *OR-gate*.
If either of the voxel intensities are above 0.5,
then the resulting output voxel intensity will be 1.
If both voxel intensities are below 0.4999, the resulting
output voxel intensity will be 0.

#### AND-gate example:
```
   combine_mrc file1.mrc,0,1 "*" file2.mrc,0,1 out_file.mrc,0,1
```

#### NAND-gate example:
This version will invert the output of the AND-gate by applying a "1,0"
threshold filter to the output voxel brightnesses.
*(Note that "0.52,0.48" would also work.)*
```
   combine_mrc file1.mrc,0,1 "*" file2.mrc,0,1 out_file.mrc,1,0
```

#### NOR-gate example:
```
   combine_mrc file1.mrc,0,1 + file2.mrc,0,1 out_file.mrc,1,0
```

Sometimes it is useful to select a ***narrow range of voxel intensities***.
When 4 numbers follow an input file name, the brightness of all the voxels from that tomogram in that file will be run through a **double-threshold** filter.
You might want to replace all voxels whose intensities lie between *0.4* and *0.5* with ***1***, and all of the brightness values *outside the wider range*
from *0.3* to *0.6* with 0.0:

```
   combine_mrc file1.mrc,0.3,0.4,0.5,0.6 + file2.mrc,0.3,0.4,0.5,0.6 out_file.mrc
```
*(Voxels whose intensities fall in the ranges between [0.3,0.4] or [0.5,0.6] will be scaled linearly).*

Graphically, the relationship between input voxel intensity and output will be:
```
 output
 intensity
  /|\
 1 |                 ________________                
   |             _,-'                `-._
   |         _,-'                        `-._
 0 |______,-'                                `-._________\ input
        0.3       0.4               0.5       0.6        / intensity
```

*Note:* If the order of thresholds are listed in ***decreasing order***, then the inverse images is generated:
For example:
```
   combine_mrc file1.mrc,0.6,0.5,0.4,0.3 + file2.mrc,0.6,0.5,0.4,0.3 out_file.mrc
```
This replaces all voxels whose intensities lie between *0.4* and *0.5* with ***0***, and all of the brightness values *outside the wider range*
from *0.3* to *0.6* with 1.
This will invert the intensity of each voxel compared to the previous example:
```
 output
 intensity
  /|\                                                   
 1 |_____                                       _________
   |     `-._                               _.-'       
   |         `-._                       _,-'             
 0 |             `-._________________,-'         ________\ input
        0.3       0.4               0.5       0.6        / intensity
```

### "not" operations
*In addition to "**and**" and "**or**" operations, you can also perform "not" operations.  You can perform the "**not**" operation on each **input** by reversing the order of the thresholds following a file name.  The example above ("file1.mrc,0.52,0.48) demonstrates how to do that.*


