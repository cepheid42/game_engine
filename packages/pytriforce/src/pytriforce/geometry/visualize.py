import numpy as np
import pyvista as pv

from .geometry import ExtrudedShape, Composite3D


def visualize_voxels(geom_object, grid, scale = None, mode = "cube", show_lines = False, show=False, show_domain = False, show_full_domain = False, cut_plane = (), pl = pv.Plotter(), color = None, opacity = 1.0, plot_origin = False):
    pv.global_theme.allow_empty_mesh = True

    boolean_map = geom_object

    if isinstance(geom_object, ExtrudedShape) or isinstance(geom_object, Composite3D):
        geom_object.voxelize(grid)
        boolean_map = geom_object.boolean_map.copy()

    # cut planes to see internal structure
    if "x" in cut_plane:
        x_index = np.argwhere(grid["x"][:,0,0] >= (grid["x"][0,0,0] + grid["x"][-1,0,0])/2)[0][0]
        boolean_map[x_index:,:,:] = False

    if "y" in cut_plane:
        y_index = np.argwhere(grid["y"][0,:,0] >= (grid["y"][0,0,0] + grid["y"][0,-1,0])/2)[0][0]
        boolean_map[:,y_index:,:] = False

    if "z" in cut_plane:
        z_index = np.argwhere(grid["z"][0,0,:] >= (grid["z"][0,0,0] + grid["z"][0,0,-1])/2)[0][0]
        boolean_map[:,:,z_index:] = False

    dx, dy, dz = grid["deltas"]
    offsets = grid["offsets"]
    points = np.vstack([grid["x"][boolean_map] + 0.5 * dx * offsets[0,0], grid["y"][boolean_map] + 0.5 * dy * offsets[0,1], grid["z"][boolean_map] + 0.5 * dz * offsets[0,2]]).T

    pdata = pv.PolyData(points)

    if scale is None:
        scale = (show_lines * 0.95 + (not show_lines))

    sphere = pv.Sphere(radius=scale/2, phi_resolution=10, theta_resolution=10)
    cube = pv.Cube(x_length = scale * dx, y_length = scale * dy, z_length = scale * dz)

    pc=pdata.glyph(geom=cube, scale = False, orient = False)

    pl.camera_position = "zx"
    pl.camera.roll = 90

    pl.add_mesh(pc, color=color, opacity = opacity)
    # pl.enable_eye_dome_lighting()

    if plot_origin:
        pl.add_mesh(pv.PolyData([0.0,0.0,0.0]).glyph(geom=cube, scale=False, orient = False),color="red")

    if show:
        pl.add_axes()
        pl.show()

    return pl


    # if show_full_domain:
    #     mlab.mesh(grid["x"][0,:,:], grid["y"][0,:,:], grid["z"][0,:,:], representation="wireframe", color=(1,0,0), opacity = 0.05)
    #     mlab.mesh(grid["x"][-1,:,:], grid["y"][-1,:,:], grid["z"][-1,:,:], representation="wireframe", color=(1,0,0), opacity = 0.05)
    #     mlab.mesh(grid["x"][:,0,:], grid["y"][:,0,:], grid["z"][:,0,:], representation="wireframe", color=(0,1,0), opacity = 0.05)
    #     mlab.mesh(grid["x"][:,-1,:], grid["y"][:,-1,:], grid["z"][:,-1,:], representation="wireframe", color=(0,1,0), opacity = 0.05)
    #     mlab.mesh(grid["x"][:,:,0], grid["y"][:,:,0], grid["z"][:,:,0], representation="wireframe", color=(0,0,1), opacity = 0.05)
    #     mlab.mesh(grid["x"][:,:,-1], grid["y"][:,:,-1], grid["z"][:,:,-1], representation="wireframe", color=(0,0,1), opacity = 0.05)