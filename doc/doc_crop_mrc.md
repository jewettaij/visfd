crop_mrc
===========
**crop_mrc** is a simple program for cropping 3D images in MRC (.REC)
file format.  It can also be used to shift (translate) the contents
of an image or extendding (padding) its size.
The syntax is:
```
   crop_mrc INPUT_FILE OUTPUT_FILE xmin xmax ymin ymax zmin zmax
```
The INPUT_FILE is cropped and a new file (OUTPUT_FILE) is created.
The last 6 arguments are nonnegative integers indicating the boundaries
of the rectangle you wish to extract from the original image.

The program can also be used to shift the position of an image and/or increase
the size of an existing image by supplying 7 additional parameters.
These 7 additional arguments allow you to add "padding" to the outside of the
cropped image.  The syntax is demonstrated below:
```
   crop_mrc INPUT_FILE OUTPUT_FILE xmin xmax ymin ymax zmin zmax \
            xpad Xpad ypad Ypad zpad Zpad B
```
This will allow you to add extra voxels (of brightness *B*)
before and after the image in the X,Y,Z directions.

As of 2022-2-03, the resulting image will be saved in
an MRC file in mode 2 (32-bit floats).


### Example 1
```
   crop_mrc orig_file.rec new_file.mrc 290 450 80 300 38 137
```
This will extract a rectangular region containing voxels
290-450 (inclusive) in the X direction,
80-300 in the Y direction, and
38 to 137 in the Z direction.


### Example 2
```
   crop_mrc orig_file.rec new_file.mrc 290 450 80 300 38 137 0 30 40 0 0 0 0
```
As before, this will extract a rectangular region from the "orig_file.rec" file.
However this time, *30* voxels of brightness *0* will be appended to the
end of the image in the X direction, and *40* voxels will be inserted
at the beginning of the image in the Y direction.
The resulting image will be size 191 x 261 x 100 voxels wide.


### NOTE: Voxel coordinates begin at 0

Unlike some programs (such as IMOD), the X,Y,Z coordinate
of the voxels in the image are counted beginning at 0 not 1.
So to extract the first 100 voxels from the image in the X,Y,Z directions,
use this:
```
   crop_mrc orig_file.rec new_file.mrc 0 99 0 99 0 99
```
*not* this:
```
   crop_mrc orig_file.rec new_file.mrc 1 100 1 100 1 100
```
Both commands will extract a 100x100x100 rectangular region from the image,
but the second command excludes the voxels located at X=0, Y=0, or Z=0.
