import numpy as np
from enum import Enum

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

class Boundary(Enum):
    PER=0       # adds nHalo nodes, defaults to 2
    PML=1
    REF=2

def update_domain(domain_shape, em_boundary, x_range, y_range, z_range, pml_depth = 10, nhalo = 2):
    em_boundary_values = np.array([x.value for x in em_boundary],int)

    pml_boundary = em_boundary_values==Boundary.PML.value
    per_boundary = em_boundary_values==Boundary.PER.value

    # compute the actual simulation mesh shape in nodes (hence the +1 at the end)
    shape = domain_shape + pml_boundary.reshape((3,2)).sum(axis=1) * pml_depth + per_boundary.reshape((3,2)).sum(axis=1) * nhalo + 1

    # todo: make this work for variable cell sizes
    dx = (x_range[1] - x_range[0]) / domain_shape[0]
    dy = (y_range[1] - y_range[0]) / domain_shape[1]
    dz = (z_range[1] - z_range[0]) / domain_shape[2]

    # compute the full mesh ranges
    x_range_full = (x_range[0] - dx * (pml_boundary[0] * pml_depth + per_boundary[0] * nhalo), x_range[1] + dx * (pml_boundary[1] * pml_depth + per_boundary[1] * nhalo))
    y_range_full = (y_range[0] - dy * (pml_boundary[2] * pml_depth + per_boundary[2] * nhalo), y_range[1] + dy * (pml_boundary[3] * pml_depth + per_boundary[3] * nhalo))
    z_range_full = (z_range[0] - dz * (pml_boundary[4] * pml_depth + per_boundary[4] * nhalo), z_range[1] + dz * (pml_boundary[5] * pml_depth + per_boundary[5] * nhalo))

    return shape, (dx, dy, dz), (x_range_full, y_range_full, z_range_full), em_boundary_values