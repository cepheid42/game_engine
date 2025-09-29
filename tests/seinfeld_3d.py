#!/usr/bin/env python3

import numpy as np
from scipy import constants

sim_name = 'seinfeld_3d'

nx = 45
ny = 2
nz = 101
nhalo = 2

xmin, xmax = -0.11, 0.11
ymin, ymax = 0.0, 0.005
zmin, zmax = -0.25, 0.25

dx = (xmax - xmin) / (nx - 1)
dy = (ymax - ymin) / (ny - 1)
dz = (zmax - zmin) / (nz - 1)

cfl = 0.95
dt = 2.5e-12
t_end = 4000 * dt
nt = int(t_end / dt) + 1

save_interval = 40
nthreads = 24
interp_order = 1

# ===== EM =====
# Periodic = 0, PML = 1, Reflecting = 2
BCS = [1, 1, 2, 2, 1, 1]
bc_str = f'{BCS[0]}zu, {BCS[1]}zu, {BCS[2]}zu, {BCS[3]}zu, {BCS[4]}zu, {BCS[5]}zu'
PMLDepth = 15
PMLGrade = 3.5
PMLAlphaMax = 0.2

# ===== Particles =====
init_temp = 10000 # eV
init_density = 8.5e27 # m^3
ppc = (10, 1, 10)
p_range = ((-5e-7, 5e-7), # pxmin, pxmax
           (ymin, ymax),  # pymin, pymax
           (-1e-5, 1e-5)) # pzmin, pzmax
particleBCs = 1

program_params = (
    '#ifndef PROGRAM_PARAM_HPP\n'
    '#define PROGRAM_PARAM_HPP\n'
    '\n'
    '#include <array>\n'
    '\n'
    f'inline constexpr auto nThreads = {nthreads};\n'
    '\n'
    f'inline constexpr auto Nx = {nx}zu;\n'
    f'inline constexpr auto Ny = {ny}zu;\n'
    f'inline constexpr auto Nz = {nz}zu;\n'
    '\n'
    f'inline constexpr std::array x_range = {{{xmin}, {xmax}}};\n'
    f'inline constexpr std::array y_range = {{{ymin}, {ymax}}};\n'
    f'inline constexpr std::array z_range = {{{zmin}, {zmax}}};\n'
    '\n'
    f'inline constexpr auto dx = {dx};\n'
    f'inline constexpr auto dy = {dy};\n'
    f'inline constexpr auto dz = {dz};\n'
    '\n'
    f'inline constexpr auto cfl   = {cfl};\n'
    f'inline constexpr auto dt    = {dt};\n'
    f'inline constexpr auto t_end = {t_end};\n'
    f'inline constexpr auto Nt    = {nt}zu;\n'
    '\n'
    f'inline constexpr auto save_interval = {save_interval}zu;\n'
    '\n'
    '/ *--------------------------------------------------------------- /\n'
    '/ -                        EM Parameters                         - /\n'
    '/ ---------------------------------------------------------------* /\n'
    'enum class EMFace { X, Y, Z };\n'
    'enum class EMSide { Lo, Hi };\n'
    '\n'
    f'inline constexpr auto PMLDepth    = {PMLDepth}zu;\n'
    f'inline constexpr auto PMLGrade    = {PMLGrade};\n'
    f'inline constexpr auto PMLAlphaMax = {PMLAlphaMax};\n'
    '//inline constexpr auto PMLKappaMax = 1.0;\n'
    '\n'
    f'inline constexpr auto nHalo = 2zu;\n'
    '\n'
    '// Periodic = 0, PML = 1, Reflecting = 2\n'
    f'inline constexpr std::array BCSelect = {{{bc_str}}};\n'
    '\n'
    '/ *--------------------------------------------------------------- /\n'
    '/ -                     Particle Parameters                      - /\n'
    '/ ---------------------------------------------------------------* /\n'
    'enum class ParticleBCType { Periodic, Outflow };\n'
    '\n'
    f'inline constexpr auto interpolation_order = {interp_order}zu;\n'
    '\n'
    f'inline constexpr auto PBCSelect = {'ParticleBCType::Periodic' if particleBCs == 0 else 'ParticleBCType::Outflow'};\n'
    '\n'
    '#endif //PROGRAM_PARAM_HPP\n'
)


param_path = '/home/cepheid/TriForce/game_engine/params/'
with open(param_path + 'program_params.hpp', 'w+') as f:
    cur_header = f.read()
    if cur_header != program_params:
        f.write(program_params)


def make_velocities(mass, T_M, num_particles, velocity=0.0):
    rng = np.random.default_rng()
    v_thermal = np.sqrt(constants.elementary_charge * T_M / mass)
    velocities = rng.normal(velocity, v_thermal, (num_particles, 3))
    return velocities

def maxwell_juttner(mass, T_M, num_particles):
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

def create_particles(name, ppc, prange, density, temp, mass, charge):
    ppc_x, ppc_y, ppc_z = ppc
    px_min, px_max = prange[0]
    py_min, py_max = prange[1]
    pz_min, pz_max = prange[2]

    pnx = int((px_max - px_min) / dx)
    pny = int((py_max - py_min) / dy)
    pnz = int((pz_max - pz_min) / dz)

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

    ppc_total = ppc[0] * ppc[1] * ppc[2]
    num_particles = ppc_total * pnx * pny * pnz

    pxs = np.arange(px0 + pdx / 2, px1, pdx)
    pys = np.arange(py0 + pdy / 2, py1, pdy)
    pzs = np.arange(pz0 + pdz / 2, pz1, pdz)
    xc, yc, zc = np.meshgrid(pxs, pys, pzs)

    weights = np.full((num_particles, 1), density * dx * dy * dz / ppc_total)
    positions = np.vstack((xc.flatten(), yc.flatten(), zc.flatten())).T
    velocities = maxwell_juttner(mass, temp, num_particles)

    output = np.hstack((positions, velocities, weights))
    header = f'{name} {mass} {charge}'
    np.savetxt(f'/home/cepheid/TriForce/game_engine/data/{name}.dat', output, delimiter=' ', header=header)

    # import matplotlib.pyplot as plt
    # plt.scatter(positions[:, 0], positions[:, 1], s=2)
    # plt.xticks(xs)
    # plt.yticks(zs)
    # plt.grid()
    # # plt.xlim([px_min, px_max])
    # # plt.ylim([pz_min, px_max])
    # plt.show()


init_density = 1.0e17 # m^3
p_range = ((-0.04, 0.04), # pxmin, pxmax meters
           (0.0, 0.005), # pymin, pymax
           (-0.15, 0.15)) # pzmin, pzmax

create_particles('electrons', (3, 1, 3), p_range, init_density, 4.0, constants.m_e, -1.0 * constants.e)
create_particles('ions', (2, 1, 2), p_range, init_density, 1.0, constants.m_p, constants.e)
