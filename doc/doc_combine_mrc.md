
combine_mrc
===========
**combine_mrc** is a program for combining two volumetric images (i.e. tomograms, both of identical size) into one image/tomogram, using a combination of addition, multiplication, and thresholding operations.  These features can be used perform binary operations between two images which are similar to "**and**" and "**or**" operations.  ("**not**" operations are also possible.  See below.) As with the "filter_mrc" program, you can also use the "-mask" argument to restrict the operation to certain voxels from the image.
*(A detailed description of the threshold and rescaling functions used in this step are provided at the end of the the "doc_filter_mrc.md" file.)*

### -rescale

If you use the **-rescale** argument, then
image voxel intensities will be shifted and rescaled so that the
minimum and maximum voxel intensities are 0 and 1.
This will be done for both the input images and well as the output image.

### Clipping the output range
*Note:*  After adding or multiplying the voxel brightnesses, the brightness of each resulting voxel is automatically clipped between 0 and 1, *before* saving the result to a file.
(In other words, the resulting brightnesses below 0 are replaced with 0, and the resulting brightnesses above 1 are replaced with 1.
*Perhaps I will provide a more general way to rescale and clip the output brightnesses later...*)


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


### Applying threshold filters to the input images:

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
0.52 |                           _.-'                 
     |                       _,-'                 
     |                   _,-'            
0.48 |________________,-'                     ________\ input
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
*In addition to "**and**" and "**or**" operations, you can also perform "not" operations.  You can perform the "**not**" operation on each **input** by reversing the order of the thresholds following a file name.  The example above ("file1.mrc,0.6,0.4) demonstrates how to do that.*


*Incidentally, "**not**" operations can also be performed on the **output** image by saving it to a file, and later using the "filter_mrc" program
to perform a threshold operation using the "-thresh2" argument
with the first threshold number greater than the second.  For example:*
```
   filter_mrc -thresh2 0.51 0.49 -in file.mrc -out not_file.mrc
```
