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


def gen_magnetic_mirror_fields(B0, L, z, shape):
    xs = np.broadcast_to(z[:-1], (shape[0], shape[1] - 1, shape[2] - 1))
    ys = np.broadcast_to(z[:-1], (shape[0] - 1, shape[1], shape[2] - 1))
    zs = np.broadcast_to(z, (shape[0] - 1, shape[1] - 1, shape[2]))
    x = -B0 * xs / L**2
    y = -B0 * ys / L**2
    z =  B0 * (1.0 + (zs / L)**2)
    return x, y, z


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