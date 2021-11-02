voxelize_mesh.py
===========
**voxelize_mesh.py** is a command-line python script which voxelizes closed
meshes.

### Warning: This is experimental software

If you run into bugs, please email me or post an issue with the github
issue tracker. -andrew 2021-6-08


### Explanation

A "mesh" is a collection of connected polygons which form a surface.  Unlike
image files, meshes contain vector geometry and are resolution independent.
The [*filter_mrc*](doc_filter_mrc.md) program (included with this repository),
together with [*PoissonRecon*](https://github.com/mkazhdan/PoissonRecon)
and [*meshlab*](http://www.meshlab.net)
can be used to create closed meshes which define the boundary surface of an
object detected from a 3D image (such as the boundary of a cell or organelle).
If you're goal is to segment the original image, you must determine which
voxels lie within this closed surface, and which voxels lie outside.
This is called "voxelization".

## Requirements

*voxelize_mesh.py* requires a python environment containing the
[pyvista](https://docs.pyvista.org),
and
[mrcfile](https://mrcfile.readthedocs.io)
packages.
These modules can be installed using
```
pip install pyvista
pip install mrcfile
```

# WARNING: Extremely memory inneficient. Use *ulimit*
This program was cobelled together using 3rd-party open-source tools:
(pyvista, vtk) which are extremely slow and memory inneficient for this task.
Currently, this program requires 160 bytes of RAM *per voxel*,
This means you may need RAM which is up to 160 times larger
than the size of the tomogram file you are segmenting.
A typical tomogram (1000x1000x600 voxels) would require 90Gb of RAM
to voxelize.
*(Unfortunately, this is beyond my control,
unless I rewrite the voxelizer from scratch.)*
For this reason, this program is barely useable on a typical desktop computer.

Machines with enough RAM can be rented from cloud vendors like
[Amazon](https://aws.amazon.com/ec2/instance-types/#Memory_Optimized)
[Google](https://cloud.google.com/compute/docs/machine-types), and
[Microsoft](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes-memory).

If you exceed the memory available on a rented computer,
the virtual machine instance can be terminated safely.
But exceeding the memory on your computer can make it unresponsive.
**Consequently, if you are running this program on your own desktop
or a shared cluster, it is strongly recommended that you use the 
[ulimit -v SIZE_IN_KB](https://ss64.com/bash/ulimit.html)
command to prevent this**
(If my understanding is correct, running "ulimit -v 14000000" beforehand
should prevent voxelize_mesh.py from consuming more than 14Gb of RAM.)

## To reduce computation time, reduce the number of vertices

A dense mesh can take a couple hours to voxelize on a full-size tomogram.
However automatically generated meshes (using *filter_mrc* and *PoissonRecon*,
for example) can be extremely dense, containing multiple triangles per voxel.
You can smooth, and simplify ("resample") the mesh using 3rd-party tools
like *meshlab*. (This can be done using the *"Filters"->"Remeshing, Simplification, and Reconstruction"->"Quadratic Edge Collapse Decimation"* menu option,
and reducing the number of faces.)  Later on, when you run *voxelize_mesh.py*,
this reduction will significantly speed up the calculation (although it is
still slow).


## Usage

Typical usage:
```
voxelize_mesh.py -m mesh.ply -i orig_image.rec -o segmented_image.rec -w 19.6
```
If you are running out of memory,
you can use the "-c" argument to reduce the image size:
```
voxelize_mesh.py -m mesh.ply -o segmented_image.rec -w 19.6 -c 330 430 80 140 50 100
```


## Details
1) You must provide the name of a file containing a closed
mesh (typically in .ply format) using the
**-m** argument.
2) You must provide the name of the volumetric image file (in MRC/REC format)
that you wish to create using the **-o** argument.
3) You can provide a 3D image file which is the same size as the
image you want to create using the **-i** argument.
(If you used *filter_mrc* to create the mesh, then you would
provide the same image from which the mesh was detected.)
Alternatively, you can directly specify the size of the image you want to
create using *either* the **-c** or **-b** arguments.
*(Note: If the object you are segmenting does not fill up the entire image
then you can reduce the memory and time needed to run the program significantly
by focusing on the region containing the feature of interest.)*
4) If the coordinates of the points in the mesh are in physical units
(such as Angstroms, or nm) then you must also specify the physical
size of each voxel using the **-w** argument.  (If you used *filter_mrc*
then make sure that you use the same **-w** argument that it uses.)


## Arguments

### -w voxel_width
The physical size of each voxel. If the mesh was created by *filter_mrc*,
this should match the **-w** argument that you used with that program.

### -m MESH_FILE
Specify the name of a file containing a closed mesh (typically in PLY format).

### -o NEW_IMAGE_FILE
Specify the name of the 3D image file (MRC or REC format)
that you wish to create.

### -i ORIG_IMAGE_FILE
Specify the name of a 3D image file (MRC or REC format)
which is the same size as the image you want to create.
(Alternatively, you can specify the size of the image file manually
using the "-c" or "-b" arguments.)

### -b x<sub>min</sub> x<sub>max</sub> y<sub>min</sub> y<sub>max</sub> z<sub>min</sub> z<sub>max</sub>
Specify the physical size of the 3-D image.
This will crop the voxelized image to a rectangular region specified
by the 6 floating point arguments.  These are in units of physical
distance, not voxels. (If you specify this argument, you do not need
to use the **-i** or **-c** arguments.)

### -c ix<sub>min</sub> ix<sub>max</sub> iy<sub>min</sub> iy<sub>max</sub> iz<sub>min</sub> iz<sub>max</sub>
Specify the size of the 3-D image (in units of voxels).
It will crop the voxelized image to the size indicated by the 6 integer
arguments.  (If you specify this argument, you do not need to use the
**-i** or **-b** arguments.)

### -s dx dy dz
Shift the coordinates of the mesh in the dx dy dz direction before voxelization.
(The dx, dy, and dz numbers are in units of voxels, not physical distance.)
This is useful if you want to move the mesh by subvoxel amounts.
(Typical usage: "-s 0.5 0.5 0.5".  *Note: Shifting the image will not effect
any of the numbers in the header of the MRC file created by this script.*)


## Installation

The *voxelize_mesh.py* program is not currently installable using pip.
To use *voxelize_mesh.py*, copy the "voxelize_mesh.py"
file to somewhere in your path (such as /usr/local/bin/)
```
sudo cp voxelize_mesh.py /usr/local/bin
```

Alternatively, you can edit your PATH variable manually to include
the subdirectory where the "voxelize_mesh.py" script is located.
Suppose the directory with this README file is named ``visfd/doc''
and it is located in your home directory:

If you use the bash shell, typically you would edit your 
`~/.profile`, `~/.bash_profile` or `~/.bashrc` files 
to contain the following lines:

```
    export PATH="$PATH:$HOME/visfd/bin/voxelize_mesh"
```

## Possible issues

1) As of 2020-12-15, I have not yet tested to make sure that the images
generated by this program match the size of the images passed to the **-i**
argument.  If not, please email me or post an issue with the github
issue tracker.  Meanwhile, you can use the **-c** or **-b** arguments
to tune the size of the output image.
2) As of 2020-12-15, I suspect that membranes detected with *filter_mrc*
are shifted by half a voxel in the x,y,z directions, but I have not tested
to see if this is true.  If you eventually plan to search for objects that lie
within 1 voxel of the boundary of a cell, a shift of half a voxel could effect
your ability to detect them.
As a workaround, you can use the **-s** argument to shift
the image created by this program.  (See explanation above.)
