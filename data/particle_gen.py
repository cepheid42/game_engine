#!/usr/bin/env python3

import numpy as np
import scipy.constants as constants

ppc_x = 4
ppc_y = 4
ppc_z = 4
ppc = ppc_x * ppc_y * ppc_z
density = 1.0e17
temp = 4 # eV
velocity = 0.0

# Number of cells
nx = 50
ny = 50
nz = 50
nhalo = 0

xmin, xmax = 0.0, 0.005
ymin, ymax = 0.0, 0.005
zmin, zmax = 0.0, 0.005

dx = (xmax - xmin) / nx
dy = (ymax - ymin) / ny
dz = (zmax - zmin) / nz

px_min, px_max = 0.0015, 0.0035
py_min, py_max = 0.0015, 0.0035
pz_min, pz_max = 0.0015, 0.0035

pnx = int((px_max - px_min) / dx)
pny = int((py_max - py_min) / dy)
pnz = int((pz_max - pz_min) / dz)

print(f'Pnx/y/z {pnx}, {pny}, {pnz}')

nc = pnx * pny * pnz

print(f'Num Cells: {nc}')

pdx = dx / ppc_x
pdy = dy / ppc_y
pdz = dz / ppc_z

xs = np.arange(xmin, xmax, dx)
ys = np.arange(ymin, ymax, dy)
zs = np.arange(zmin, zmax, dz)

px0 = xs[np.searchsorted(xs, px_min)]
px1 = xs[np.searchsorted(xs, px_max)]

py0 = ys[np.searchsorted(ys, py_min)]
py1 = ys[np.searchsorted(ys, py_max)]

pz0 = zs[np.searchsorted(zs, pz_min)]
pz1 = zs[np.searchsorted(zs, pz_max)]

cell_volume = dx * dy * dz
num_particles = ppc * nc

print(f'Number of particles per group: {num_particles}')

# ===== Velocities =====
def make_velocities(mass, T_M):
    rng = np.random.default_rng()
    v_thermal = np.sqrt(constants.elementary_charge * T_M / mass)
    velocities = rng.normal(velocity, v_thermal, (num_particles, 3))
    return velocities

def maxwell_juttner(mass, T_M):
    rng = np.random.default_rng()
    T_norm = constants.elementary_charge * T_M / (mass * constants.c**2)
    a, b, R0 = 0.56, 0.35, 0.95
    root2 = np.sqrt(2)
    w3 = np.sqrt(np.pi)
    w4 = a * np.sqrt(2 * T_norm)
    w5 = 1.5 * np.sqrt(np.pi) * b * T_norm
    w6 = (2 * T_norm)**1.5
    s_sum = w3 + w4 + w5 + w6
    pi3 = w3 / s_sum
    pi4 = w4 / s_sum
    pi5 = w5 / s_sum
    def sample():
        def R(xx):
            return (1 + xx) * np.sqrt(xx + 2) / (root2 + (a * np.sqrt(xx)) + (b * root2 * xx) + xx ** (3 / 2))
        while True:
            X1 = rng.uniform()
            X2 = rng.uniform()
            i = 6
            if X1 < pi3:
                i = 3
            elif X1 < pi3 + pi4:
                i = 4
            elif X1 < pi3 + pi4 + pi5:
                i = 5
            x = rng.gamma(i / 2, T_norm)
            if X2 < R0 or X2 < R(x):
                break

        X3 = rng.uniform()
        X4 = rng.uniform()
        u_mag = constants.c * np.sqrt(x * (x + 2))
        u = u_mag * np.asarray([(2 * X3 - 1),
                                 2.0 * np.sqrt(X3 * (1.0 - X3)) * np.cos(2.0 * np.pi * X4),
                                 2.0 * np.sqrt(X3 * (1.0 - X3)) * np.sin(2.0 * np.pi * X4)])
        return u / np.sqrt(1 + (u @ u) / constants.c**2)

    velocities = np.zeros((num_particles, 3))
    for p in range(num_particles):
        velocities[p, :] = sample()
    return velocities

# ===== Weights =====
weights = np.full((num_particles, 1), density * cell_volume / ppc)

#===== Positions =====
pxs = np.arange(px0 + pdx / 2, px1, pdx)
pys = np.arange(py0 + pdy / 2, py1, pdy)
pzs = np.arange(pz0 + pdz / 2, pz1, pdz)

xc, yc, zc = np.meshgrid(pxs, pys, pzs)
positions = np.vstack((xc.flatten(), yc.flatten(), zc.flatten())).T

# v_e = maxwell_juttner(constants.m_e, temp)
# v_i = maxwell_juttner(100.0 * constants.m_e, temp)

v_e = make_velocities(constants.m_e, temp)
v_i = make_velocities(constants.m_e, temp)

# positions[:, 1] = ymax / 2.0

# import matplotlib.pyplot as plt
# plt.scatter(positions[:, 0], positions[:, 1], s=2)
# plt.xticks(xs)
# plt.yticks(zs)
# plt.grid()
# # plt.xlim([px_min, px_max])
# # plt.ylim([pz_min, px_max])
# plt.show()

# ===== Output =====
electrons = np.hstack((positions, v_e, weights))
ions = np.hstack((positions, v_i, weights))

np.savetxt('/home/cepheid/TriForce/game_engine/data/electrons.dat', electrons, delimiter=' ')
np.savetxt('/home/cepheid/TriForce/game_engine/data/ion_slab.dat', ions, delimiter=' ')

