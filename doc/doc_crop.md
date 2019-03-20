crop_mrc
===========
**crop_mrc** is a crude, simple program for cropping 3D images in MRC (.REC)
file format.  The syntax is:
```
   crop_mrc INPUT_FILE OUTPUT_FILE xmin xmax ymin ymax zmin zmax
```
The INPUT_FILE is cropped and a new file (OUTPUT_FILE) is created.
The last 6 arguments are nonnegative integers indicating the boundaries
of the rectangle you wish to extract from the original image.

### Example:
```
   crop_mrc orig_file.rec new_file.mrc 290 450 80 300 38 137
```
This will extract a rectangular region containing voxels
290-450 (inclusive) in the X direction,
80-300 in the Y direction, and
38 to 137 in the Z direction.

Note: If you are using IMOD to locate voxel coordinates,
       keep in mind that IMOD uses
      "1-indexing".  This means that it 
      prints voxel coordinates starting at (1,1,1).
      In this program, voxel coordinates begin at (0,0,0).
