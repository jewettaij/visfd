sum_voxels
===========
**sum_voxels** computes the sum of the voxels in a 3D image
(in MRC/.REC format).  The syntax is:
```
   sum_voxels INPUT_FILE [-mask MASK_FILE] [-w VOXELWIDTH] [-thresh THRESHOLD] [-volume]
```

## Arguments:


### -ave
The **-ave** argument will print the average voxel brightness.


### -stddev
The **-stddev** argument will print the standard deviation
of voxel brightnesses.


### -volume
The **-volume** argument will multiply the resulting sum of the voxel
brightness by the volume of each voxel.
The width of each voxel is either inferred from
the input file or specified manually using the **-w** argument.
(This argument is called "-volume" because this sum happens to be a volume
 when the voxel brightnesses are allways either 0 or 1.
 To insure that this is so, you can use the
 **-thresh**, **-thresh-range**, **-thresh2**, or **-thresh4**, arguments.)


### -w VOXELWIDTH
The **-w VOXELWIDTH** argument allows the user to specify the width
of each voxel in physical units (ie in Angstroms or nm).
(This argument is ignored unless the **-volume** argument is also specified.)


### -mask  mask_file.rec
  Specify an image file indicating which voxels should be ignored.
  This image should be the same size as the image INPUT_FILE.
  Voxels in this image with brightness 0 will be ignored.
  Typically these are voxels which lie outside the region of interest
  (typically a cell, organelle, or compartment subvolume).


### -mask-select  INTENSITY_VALUE

If the **-mask-select** argument is specified, then instead of considering all
voxels with non-zero values from the "mask" image,
only voxels whose mask intensity equals the number following
this argument will belong to the mask.
All other voxels will be ignored during filtering.
*(Note: This disables "soft" masking.  In other words,
  during the process of filtering, all selected voxels will be weighted equally
  during filtering, and all others will be completely ignored.)*

***NOTE FOR IMOD USERS:***
Some MRC files use signed bytes.
This is a problem because for these files
because IMOD reports all voxel brightnesses
in the range from 0-255 instead of -128-127.
This means that if IMOD reports that a voxel
has brightness of "1", (for example),
the true brightness stored in the file might be -127.
In that case using "-mask-select 1" will fail.
No other software I know does this.
You can tell if your MASK image file potentially has this problem
by using IMOD's "header" program, for example:
```
   header mask_file.rec
```
If the "Map mode" reported by "header" equals 0
then I suggest that you should convert the file
to a floating-point MRC format before using it with this software.
One way to do this is to use the
"convert_to_float" program (distributed with visfd):
```
   convert_to_float mask_file.rec new_file.rec
```
(If you prefer, you can also use IMOD's "newstack" program:
 "newstack -input mask_file.rec -output new_file.rec -mode 2")
Either way, this insures that the brightnesses in the "new_file.rec"
will be interpreted the same way in IMOD and other programs.
*Afterwards, if IMOD says that a particular voxel in this file
("new_file.rec") has value "1", then you can safely use "-mask-select 1"*



### -thresh threshold

For convenience, threshold operation can be applied
(using the **-thresh**, "**-thresh-range**",
 **-thresh2**, or **-thresh4**, arguments)
so that the voxels intensities are changed to either 0 or 1
*before* the sum is calculated.
(This will convert all the voxels above/below the threshold(s) to 1
 and the others to 0.)


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

### -thresh2 thresh_a thresh_b

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

*Note: You can use the "-thresh-range outA outB" argument to scale the
       resulting voxel intensities from outA to outB (instead of from 0 to 1).*


### -thresh4 thresh01_a thresh01_b thresh10_a thresh10_b

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


*Note: You can use the "-thresh-range outA outB" argument to scale the
       resulting voxel intensities from outA to outB (instead of from 0 to 1).*


*Note:* The **-w**, **-mask**, **-mask-select**,
**-thresh**, **-thresh-range**, **-thresh2**, **-thresh4**, **-clip**
arguments are shared with the *filter_mrc* program
and are explained in more detail [here](./doc_filter_mrc.md).


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
