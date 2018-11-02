sum_voxels
===========
**sum_voxels** computes the sum of the voxels in a 3D image
(in MRC/.REC format).  The syntax is:
```
   sum_voxels INPUT_FILE [-mask MASK_FILE] [-w VOXELWIDTH] [-thresh THRESHOLD] [-volume]
```

The sum can be restricted to certain regions
(by using the "-mask" and "-mask-select" arguments).

For convenience, threshold operation can be applied
(using the "-thresh", "-thresh2", and "-thresh4" arguments)
so that the voxels intensities are changed to either 0 or 1
before the sum is calculated.
(This will convert all the voxels above/below the threshold(s) to 1
 and the others to 0.  In some cases, this can make it easier to
 use the program to estimate volumes.)

*Note:* The **-w**, **-mask**, **-mask-select**,
**-thresh**, **-thresh2**, **-thresh4**, **-clip**
arguments are shared with the *filter_mrc* program
and are explained in more detail [here](./doc_filter_mrc.md).
