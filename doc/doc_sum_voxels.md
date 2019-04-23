sum_voxels
===========
**sum_voxels** computes the sum of the voxels in a 3D image
(in MRC/.REC format).  The syntax is:
```
   sum_voxels INPUT_FILE [-mask MASK_FILE] [-w VOXELWIDTH] [-thresh THRESHOLD] [-volume]
```

## Arguments:

#### -mask  mask_file.rec
  Specify an image file indicating which voxels should be ignored.
  This image should be the same size as the image INPUT_FILE.
  Voxels in this image with brightness 0 will be ignored.
  Typically these are voxels which lie outside the region of interest
  (typically a cell, organelle, or compartment subvolume).

#### -mask-select  INTENSITY_VALUE

If the **-mask-select** argument is specified, then instead of considering all
voxels with non-zero values from the "mask" image,
only voxels whose mask intensity equals the number following
this argument will belong to the mask.
All other voxels will be ignored during filtering.
*(Note: This disables "soft" masking.  In other words,
  during the process of filtering, all selected voxels will be weighted equally
  during filtering, and all others will be completely ignored.)*

#### -ave
The **-ave** argument will print the average voxel brightness.

#### -stddev
The **-stddev** argument will print the standard deviation
of voxel brightnesses.

#### -volume
The **-volume** argument will multiply the resulting sum of the voxel
brightness by the volume of each voxel.
The width of each voxel is either inferred from
the input file or specified manually using the **-w** argument.
(This argument is called "-volume" because this sum happens to be a volume
 when the voxel brightnesses are allways either 0 or 1.
 To insure that this is so, you can use the
 **-thresh**, **-thresh-range**, **-thresh2**, or **-thresh4**, arguments.)

#### -w VOXELWIDTH
The **-w VOXELWIDTH** argument allows the user to specify the width
of each voxel in physical units (ie in Angstroms or nm).
(This argument is ignored unless the **-volume** argument is also specified.)

#### -thresh
For convenience, threshold operation can be applied
(using the **-thresh**, "**-thresh-range**",
 **-thresh2**, or **-thresh4**, arguments)
so that the voxels intensities are changed to either 0 or 1
*before* the sum is calculated.
(This will convert all the voxels above/below the threshold(s) to 1
 and the others to 0.)

*Note:* The **-w**, **-mask**, **-mask-select**,
**-thresh**, **-thresh-range**, **-thresh2**, **-thresh4**, **-clip**
arguments are shared with the *filter_mrc* program
and are explained in more detail [here](./doc_filter_mrc.md).
