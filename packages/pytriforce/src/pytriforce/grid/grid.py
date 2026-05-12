from dataclasses import dataclass
import numpy as np

from ..header_utils import *


class Mesh:
    def __init__(self, field_type, coords, deltas, staggered):
        self.field_type = field_type
        self.x, self.y, self.z = np.meshgrid(*coords, indexing='ij')
        self.deltas = np.array(deltas)
        self.offsets = np.zeros((3, 3), dtype=int)

        if staggered:
            if field_type in ['E', 'A', 'eps']:
                self.offsets = 0.5 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
            elif field_type in ['B', 'mu']:
                self.offsets = 0.5 *np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=int)
            elif field_type == 'geom':
                self.offsets = 0.5 * np.array([[1, 1, 1]], dtype=int)


@dataclass
class Domain:
    dims: tuple
    x_range: tuple
    y_range: tuple
    z_range: tuple
    deltas: tuple = None
    boundaries: tuple = (0, 0, 0, 0, 0, 0)

    def __post_init__(self):
        dims = np.array(self.dims)
        assert np.all(dims >= 2)

        xr = np.array(self.x_range)
        yr = np.array(self.y_range)
        zr = np.array(self.z_range)
        bcs = np.array(self.boundaries)

        dx = (xr[1] - xr[0]) / (dims[0] - 1)
        dy = (yr[1] - yr[0]) / (dims[1] - 1)
        dz = (zr[1] - zr[0]) / (dims[2] - 1)

        x_min = xr[0] - bcs[0] * dx
        x_max = xr[1] + bcs[1] * dx
        x_dims = dims[0] + bcs[0] + bcs[1]

        y_min = yr[0] - bcs[2] * dy
        y_max = yr[1] + bcs[3] * dy
        y_dims = dims[1] + bcs[2] + bcs[3]

        z_min = zr[0] - bcs[4] * dz
        z_max = zr[1] + bcs[5] * dz
        z_dims = dims[2] + bcs[4] + bcs[5]

        self.dims = (x_dims, y_dims, z_dims)
        self.deltas = (dx, dy, dz)
        self.x_range = (x_min, x_max)
        self.y_range = (y_min, y_max)
        self.z_range = (z_min, z_max)


    def coords(self):
        return (np.linspace(self.x_range[0], self.x_range[1], self.dims[0], endpoint=True),
                np.linspace(self.y_range[0], self.y_range[1], self.dims[1], endpoint=True),
                np.linspace(self.z_range[0], self.z_range[1], self.dims[2], endpoint=True))

    def get_mesh(self, field_type, staggered):
        return Mesh(field_type, self.coords(), self.deltas, staggered)

    def __repr__(self):
        return str(
            section_label('Domain Parameters') +
            constexpr_declaration('x_collapsed', str(self.dims[0] == 2).lower()) +
            constexpr_declaration('y_collapsed', str(self.dims[1] == 2).lower()) +
            constexpr_declaration('z_collapsed', str(self.dims[2] == 2).lower()) +
            '\n' +
            constexpr_declaration('Nx', f'{self.dims[0]}zu') +
            constexpr_declaration('Ny', f'{self.dims[1]}zu') +
            constexpr_declaration('Nz', f'{self.dims[2]}zu') +
            '\n' +
            constexpr_declaration('x_range', tuple_to_array(self.x_range), typestr='std::array') +
            constexpr_declaration('y_range', tuple_to_array(self.y_range), typestr='std::array') +
            constexpr_declaration('z_range', tuple_to_array(self.z_range), typestr='std::array') +
            '\n' +
            constexpr_declaration('dx', f'{self.deltas[0]}zu') +
            constexpr_declaration('dy', f'{self.deltas[1]}zu') +
            constexpr_declaration('dz', f'{self.deltas[2]}zu')
        )
