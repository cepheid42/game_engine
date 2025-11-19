from collections.abc import Callable, Iterable

import numpy as np
from scipy import constants
from adios2 import Stream


# ===== Position Sampling Utilities =====
# =======================================
def generate_normalized_positions_1d(ppc):
    return np.linspace(1.0 / ppc, 1.0, ppc) - 0.5 / ppc

def generate_normalized_positions_3d(ppc):
    return [points.flatten() for points in np.meshgrid(generate_normalized_positions_1d(ppc[0]),
                                                       generate_normalized_positions_1d(ppc[1]),
                                                       generate_normalized_positions_1d(ppc[2]),
                                                       indexing='ij')]


def generate_density_mapping(density : float | Callable | Iterable, xcenter, ycenter, zcenter):
    assert(xcenter.shape == ycenter.shape == zcenter.shape) # check shapes for self-consistency

    # if centers aren't already meshed, mesh them
    if xcenter.ndim == 1:
        xmesh, ymesh, zmesh = np.meshgrid(xcenter, ycenter, zcenter, sparse=True)
    elif xcenter.ndim == 3:
        xmesh, ymesh, zmesh = xcenter, ycenter, zcenter
    else:
        raise ValueError

    if isinstance(density, Iterable):
        return density
    elif isinstance(density, Callable):
        return density(xmesh, ymesh, zmesh)
    else:
        return np.full(xmesh.shape, density)



