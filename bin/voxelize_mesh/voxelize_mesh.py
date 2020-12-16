#!/usr/bin/env python

import sys
import argparse
import numpy as np
import pyvista as pv
from vtk.util import numpy_support as vtknp
import mrcfile

g_filename = __file__.split('/')[-1]
g_module_name = g_filename
if g_filename.rfind('.py') != -1:
    g_module_name = g_filename[:g_filename.rfind('.py')]
g_date_str = '2020-12-13'
g_version_str = '0.0.1'
g_program_name = g_filename
#sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+' ')


class InputError(Exception):
    """ A generic exception object containing a string for error reporting.
        (Raising this exception implies that the caller has provided
         a faulty input file or argument.)

    """
    def __init__(self, err_msg):
        self.err_msg = err_msg
    def __str__(self):
        return self.err_msg
    def __repr__(self):
        return str(self)



def voxelize_numpy(mesh,
                   density=None,
                   check_surface=True,
                   bounds=None,
                   shift=(0.0,0.0,0,0)):
    """Voxelize mesh to create a 3D numpy array of bools (True, False),
       indicating whether the corresponding voxel is within the closed surface
       formed by the mesh.

    Parameters
    ----------
    density : float
        The uniform size of the voxels. Defaults to 1/100th of the mesh length.

    check_surface : bool
        Specify whether to check the surface for closure. If on, then the
        algorithm first checks to see if the surface is closed and
        manifold. If the surface is not closed and manifold, a runtime
        error is raised.

    bounds : Size of image in units of physical distance:
             (x_min, x_max, y_min, y_max, z_min, z_max)
             By default, mesh.bounds is used.
    """
    if not pv.is_pyvista_dataset(mesh):
        mesh = pv.wrap(mesh)
    if density is None:
        density = mesh.length / 100
    if bounds == None:
        bounds = mesh.bounds
    x_min, x_max, y_min, y_max, z_min, z_max = bounds
    x = np.arange(x_min, x_max, density)
    y = np.arange(y_min, y_max, density)
    z = np.arange(z_min, z_max, density)
    x, y, z = np.meshgrid(x, y, z)
    # Create unstructured grid from the structured grid
    grid = pv.StructuredGrid(x, y, z)
    ugrid = pv.UnstructuredGrid(grid)
    # Get part of the mesh within the mesh's bounding surface.
    selection = ugrid.select_enclosed_points(mesh.extract_surface(),
                                             tolerance=0.0,
                                             check_surface=check_surface)
    mask = selection.point_arrays['SelectedPoints'].view(np.bool_)
    # Mask contains an array of True, False values indicating whether
    # the corresponding voxel is inside the closed mesh.
    # However it is a 1-dimensional VtkArray, not a 3D numpy array.
    # It must be converted to a 1D numpy array and converted to 3D with reshape.
    # It also must be transposed (x,y,z axes swapped).  (I have no idea why.)
    data = np.transpose(vtknp.vtk_to_numpy(mask).reshape(grid.dimensions[2],
                                                         grid.dimensions[1],
                                                         grid.dimensions[0]),
                        (0,2,1))
    return data



def main():
    try:
        ap = argparse.ArgumentParser()
        ap.add_argument('-m', '--mesh', dest='fname_mesh', required=True,
                        help='file containing closed mesh (eg. a ".ply" file)')
        ap.add_argument('-o', '--out', dest='fname_out', required=True,
                        help='name of the output file that will contain the voxelized mesh (mrc/rec format)')
        ap.add_argument('-i', '--in', dest='fname_mrc_orig', required=False,
                        help='file name of an MRC (or REC) file with the same size as the target. (Typically it is the original image in which the mesh surface was detected.)')
        ap.add_argument('-w', '--width', dest='voxel_width', required=False, type=float,
                        help='-w (or --width) should be followed by the voxel width (default 1)')
        ap.add_argument('-c', '--crop', dest='ibounds', required=False, type=float, nargs=6,
                        help='6 numbers indicating desired boundaries of the resulting cropped image: xmin xmax ymin ymax zmin zmax.  (These numbers are in units of voxels. Note: This will override the image size determined from the "-i" or "--in" argument.)')
        ap.add_argument('-b', '--bounds', dest='bounds', required=False, type=float, nargs=6,
                        help='6 numbers indicating desired image size: xmin xmax ymin ymax zmin zmax.  (If the voxel width is known, these numbers are in units of distance, not voxels. Note: This will override the image size determined from the "-i" or "--in" argument.)')
        ap.add_argument('-s', '--shift', dest='shift', required=False, type=float, nargs=3,
                        help='3 numbers indicating a shift in the x,y,z coordinates of the mesh before voxelization.  (These numbers are in units of voxels, not physical distance.)')
        args = ap.parse_args()


        # Now process the argument list

        # Read the mesh file
        try:
            mesh = pv.read(args.fname_mesh)
        except IOError:
            raise InputError('Error: Unable to open file "'+
                             args.fname_mesh+'" for reading.\n')

        # Determine the voxel width
        voxel_width = args.voxel_width
        if voxel_width == 0.0:
            raise InputError('Error: voxel width cannot be non-zero\n')

        # Determine how big the image is
        bounds = args.bounds

        # Make sure the user does not accidentally erase their original tomogram
        if args.fname_out == args.fname_mrc_orig:
            raise InputError('Error: Input and output image files cannot have the same name.\n')

        # Did the user specify an input image?
        # If so, use it to determine the output image size (and voxel width)
        if args.fname_mrc_orig:
            try:
                with mrcfile.open(args.fname_mrc_orig, 'r') as mrcdata:
                    if voxel_width == None:
                        # mrcdata.voxel_size contains the width of the voxel
                        # (Sometimes it is a numpy array with 3 elements.)
                        if hasattr(mrcdata.voxel_size, 'x'):
                            voxel_width = float(mrcdata.voxel_size.x)
                            # Due to roundoff error,the following is always true
                            # COMMENTING OUT:
                            #if ((mrcdata.voxel_size.x != mrcdata.voxel_size.y)
                            #    or
                            #    (mrcdata.voxel_size.y != mrcdata.voxel_size.z)):
                            #    sys.stderr.write('Warning: The voxels in file "'+args.fname_mrc_orig+'"\n'
                            #                     '         have a different width in the X,Y,Z directions.\n'
                            #                     '         Using the voxel width in the X direction\n')
                        else:
                            voxel_width = voxel_size
                    bounds = (0.0, mrcdata.header.nx*voxel_width,
                              0.0, mrcdata.header.ny*voxel_width,
                              0.0, mrcdata.header.nz*voxel_width)
                    mrcdata.close()
            except IOError:
                raise InputError('Error: Unable to open file "'+
                                 args.fname_mrc_orig+'" for reading.\n')
        elif voxel_width == None:
            voxel_width = 1
        assert(voxel_width != None)

        # Alternatively, did the user specify the bounds in units of voxels?
        if args.ibounds:
            bounds = (args.ibounds[0]*voxel_width,
                      (args.ibounds[1]+0.99)*voxel_width,
                      args.ibounds[2]*voxel_width,
                      (args.ibounds[3]+0.99)*voxel_width,
                      args.ibounds[4]*voxel_width,
                      (args.ibounds[5]+0.99)*voxel_width)

        # Did the user want us to shift the x,y,z coordinates of the mesh?
        if args.shift:
            bounds[0] -= args.shift[0]*voxel_width
            bounds[1] -= args.shift[0]*voxel_width
            bounds[2] -= args.shift[1]*voxel_width
            bounds[3] -= args.shift[1]*voxel_width
            bounds[4] -= args.shift[2]*voxel_width
            bounds[5] -= args.shift[2]*voxel_width

        # Now convert the mesh into an image whose (physical) size is "bounds"
        voxels = voxelize_numpy(mesh,
                                density=voxel_width,
                                check_surface=True,
                                bounds=bounds)

        # Now save the resulting numpy array as an MRC file
        mrcdata = mrcfile.new(args.fname_out, overwrite=True)
        mrcdata.voxel_size = voxel_width
        mrcdata.set_data(voxels)
        mrcdata.close()


        # Note to self: Here's a way to visualize the results using pyvista:
        #voxels["density"] = np.full(voxels.n_cells, 3.65) # 3.65 is arbitrary
        #voxels.plot(scalars="density")
        # Alternative method:
        #voxels.compute_implicit_distance(mesh, inplace=True)
        #contours = voxels.contour(20, scalars="implicit_distance")
        #p.add_mesh(contours, opacity=0.5, scalars="implicit_distance")
        #p.show()
        #p = pv.Plotter()

    except (InputError, ValueError) as err:
        sys.stderr.write('\n' + str(err) + '\n')
        sys.exit(-1)

    return


if __name__ == '__main__':
    main()
