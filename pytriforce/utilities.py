import os
import subprocess
import numpy as np
from dataclasses import dataclass, field
from scipy import constants


def calculate_cfl(dt, dx, dy, dz):
    return constants.c * dt * np.sqrt(1.0/dx**2 + 1.0/dy**2 + 1.0/dz**2)


def calculate_dt(cfl, dx, dy, dz):
    return cfl / (constants.c * np.sqrt(1.0/dx**2 + 1.0/dy**2 + 1.0/dz**2))


def velocity_from_gamma(gamma):
    return constants.c * np.sqrt(1.0 - 1.0 / gamma**2)


def gamma_from_velocity(v):
    return 1.0 / np.sqrt(1.0 - (v / constants.c)**2)


def get_rotation_matrix(angle, axis = (0.0, 0.0, 1.0)):
    norm = np.linalg.norm(axis)
    u = axis / (norm + (norm==0))
    cos_th = np.cos(angle)
    cos_thm = 1 - cos_th
    sin_th = np.sin(angle)

    return np.array(
        [[u[0] * u[0] * cos_thm +		 cos_th, u[0] * u[1] * cos_thm - u[2] * sin_th, u[0] * u[2] * cos_thm + u[1] * sin_th],
         [u[0] * u[1] * cos_thm + u[2] * sin_th, u[1] * u[1] * cos_thm +		cos_th, u[1] * u[2] * cos_thm - u[0] * sin_th],
         [u[0] * u[2] * cos_thm - u[1] * sin_th, u[2] * u[1] * cos_thm + u[0] * sin_th, u[2] * u[2] * cos_thm +		   cos_th]]
    )


# def get_field_points_component(component, grid, grid_shape, field_object):
#     # get field component index
#     index = 0
#     if "y" in component:
#         index = 1
#     elif "z" in component:
#         index = 2
#
#     # set the points
#     field_object.set_points(grid)
#
#     # get the field at the points
#     if "B" in component:
#         # gets magnetic field
#         field_points = field_object.B()
#
#     else:
#         # gets magnetic vector potential
#         field_points = field_object.A()
#
#     # reshape the field points, select components, and make contiguous
#     return np.ascontiguousarray(field_points.T.reshape([3, *grid_shape])[index])


def get_field_points_components(component, grid, grid_shape, field_object):
    # set the points
    field_object.set_points(grid)

    # get the field at the points
    if "B" in component:
        # gets magnetic field
        field_points = field_object.B()

    else:
        # gets magnetic vector potential
        field_points = field_object.A()

    # reshape the field points, select components, and make contiguous

    fields = [
        np.ascontiguousarray(field_points.T.reshape([3, *grid_shape])[0]),
        np.ascontiguousarray(field_points.T.reshape([3, *grid_shape])[1]),
        np.ascontiguousarray(field_points.T.reshape([3, *grid_shape])[2])
    ]

    # todo: condense this a bit
    # check for nans and replace them with an average of the adjacent values
    for i, field in enumerate(fields):
        if np.isnan(field).any():
            for j, nanid in enumerate(np.argwhere(np.isnan(field))):
                # get adjacent points
                points = []
                if nanid[0] - 1 >= 0:
                    if not np.isnan(field[nanid[0] - 1, nanid[1], nanid[2]]):
                        points.append(field[nanid[0] - 1, nanid[1], nanid[2]])
                if nanid[0] + 1 < field.shape[0]:
                    if not np.isnan(field[nanid[0] + 1, nanid[1], nanid[2]]):
                        points.append(field[nanid[0] + 1, nanid[1], nanid[2]])
                if nanid[1] - 1 >= 0:
                    if not np.isnan(field[nanid[0], nanid[1] - 1, nanid[2]]):
                        points.append(field[nanid[0], nanid[1] - 1, nanid[2]])
                if nanid[1] + 1 < field.shape[1]:
                    if not np.isnan(field[nanid[0], nanid[1] + 1, nanid[2]]):
                        points.append(field[nanid[0], nanid[1] + 1, nanid[2]])
                if nanid[2] - 1 >= 0:
                    if not np.isnan(field[nanid[0], nanid[1], nanid[2] - 1]):
                        points.append(field[nanid[0], nanid[1], nanid[2] - 1])
                if nanid[2] + 1 < field.shape[2]:
                    if not np.isnan(field[nanid[0], nanid[1], nanid[2] + 1]):
                        points.append(field[nanid[0], nanid[1], nanid[2] + 1])

                field[*nanid] = np.array(points).mean()

    return fields

def generate_grid(x_range_full, y_range_full, z_range_full, shape, deltas, field_type = "B", staggered_fields = False):
    # unpack deltas
    dx,dy,dz = deltas

    # create the x, y, z ranges
    xs = np.linspace(*x_range_full, shape[0])
    ys = np.linspace(*y_range_full, shape[1])
    zs = np.linspace(*z_range_full, shape[2])

    # check that all deltas are consistent
    # in principle this should work for spatially-varying deltas too
    # this will incidentally catch errors with the shape too
    assert((np.abs(dx - np.diff(xs)) < 1e-15).all())
    assert((np.abs(dy - np.diff(ys)) < 1e-15).all())
    assert((np.abs(dz - np.diff(zs)) < 1e-15).all())

    offsets = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=int)
    if staggered_fields:
        if "E" in field_type or "A" in field_type or "eps" in field_type:
            offsets = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
        elif "B" in field_type:
            offsets = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=int)
        elif "geom" in field_type:
            offsets = np.array([[1, 1, 1]], dtype=int)

    x, y, z = np.meshgrid(xs, ys, zs, indexing="ij")

    grids = {
        "x": x,
        "y": y,
        "z": z,
        "deltas": deltas,
        "offsets": offsets
    }

    return grids

def generate_field_points(x_range, y_range, z_range, shape, deltas, field_object, field_type = "B", staggered_fields = False):
    # xs = np.linspace(*x_range, shape[0])
    # ys = np.linspace(*y_range, shape[1])
    # zs = np.linspace(*z_range, shape[2])
    # offsets = np.zeros((3, 3), dtype=int)
    # if staggered_fields:
    #     if field_type in ['E', 'A', 'eps']:
    #         offsets = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
    #     elif field_type == 'B':
    #         offsets = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=int)
    #     elif field_type == 'geom':
    #         offsets = np.array([[1, 1, 1]], dtype=int)
    #
    # x, y, z = np.meshgrid(xs, ys, zs, indexing='ij')
    # flat_grid = np.hstack([x.flatten(), y.flatten(), z.flatten()])

    grid = generate_grid(x_range, y_range, z_range, shape, deltas, field_type, staggered_fields)

    # save the shape of the grids for the different components
    grid_shape = grid["x"].shape

    # create the n x 3 coordinate arrays that simsopt's field objects require
    grid_flat = np.ascontiguousarray(np.vstack([grid["x"].flatten(), grid["y"].flatten(), grid["z"].flatten()]).T)

    # obtain field points grids (properly reshaped) and return
    field_x, field_y, field_z = get_field_points_components(field_type, grid_flat, grid_shape, field_object)

    components = {
        field_type : {
            "x": field_x,
            "y": field_y,
            "z": field_z,
        }
    }

    return components, grid


def create_data_dir(data_path):
    if not os.path.exists(data_path):
        print(f'Creating simulation data directory "{data_path}"...')
        os.makedirs(data_path)


def compile_project(build_path, check=True, output=False):
    command = ['meson', 'compile', '-C', build_path, '-j4']
    if output:
        subprocess.run(command, check=check)
    else:
        subprocess.run(
            command,
            check=check,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT
        )


def run_project(exec_path, check=True, output=False):
    if output:
        subprocess.run(exec_path, check=check)
    else:
        subprocess.run(
            exec_path,
            check=check,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT
        )


@dataclass
class ParticlePlotData:
    velocities : np.ndarray = field(default_factory=lambda: np.zeros(shape=(1, 3), dtype=np.float64))
    positions : np.ndarray = field(default_factory=lambda: np.zeros(shape=(1, 3), dtype=np.float64))
    weights : np.ndarray = field(default_factory=lambda: np.zeros(shape=(1,), dtype=np.float64))
    gammas : np.ndarray = field(default_factory=lambda: np.zeros(shape=(1,), dtype=np.float64))
    times : np.ndarray = field(default_factory=lambda: np.zeros(shape=(1,), dtype=np.float64))