from dataclasses import dataclass
import numpy as np

from ..params.header_utils import *


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

    def __post_init__(self):
        assert self.dims[0] >= 2 and self.dims[1] >= 2 and self.dims[2] >= 2
        calc_delta = lambda a, na: (a[1] - a[0]) / na
        if self.deltas is None:
            self.deltas = (
                calc_delta(self.x_range, self.dims[0] - 1),
                calc_delta(self.y_range, self.dims[1] - 1),
                calc_delta(self.z_range, self.dims[2] - 1)
            )

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
