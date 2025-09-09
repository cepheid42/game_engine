#!/usr/bin/env python3

# import subprocess
import math
from scipy import constants

sim_name = 'lsi'

nx = 11
ny = 11
nz = 11
nhalo = 0

xmin, xmax = 0.0, 0.001
ymin, ymax = 0.0, 0.001
zmin, zmax = 0.0, 0.001

dx = (xmax - xmin) / (nx - 1)
dy = (ymax - ymin) / (ny - 1)
dz = (zmax - zmin) / (nz - 1)

cfl = 0.95
dt = cfl / (constants.c * math.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2))
nt = 400
t_end = nt * dt

# dt = 4.0e-17
# t_end = 3.0e-13
# nt = int(t_end / dt) + 1
# cfl = constants.c * dt * math.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

save_interval = 1
nthreads = 1
interp_order = 2

PMLDepth = 0
PMLGrade = 3.5
PMLAlphaMax = 0.2

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
    f'inline constexpr auto NHalo = {nhalo}zu;\n'
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
    f'inline constexpr auto interpolation_order = {interp_order}zu;\n'
    '\n'
    f'inline constexpr auto Ncx = {nx - 1}zu;\n'
    f'inline constexpr auto Ncy = {ny - 1}zu;\n'
    f'inline constexpr auto Ncz = {nz - 1}zu;\n'
    '\n'
    f'inline constexpr auto PMLDepth    = {PMLDepth}zu;\n'
    f'inline constexpr auto PMLGrade    = {PMLGrade};\n'
    f'inline constexpr auto PMLAlphaMax = {PMLAlphaMax};\n'
    f'//inline constexpr auto PMLKappaMax = 1.0;\n'
    '\n'
    '#endif //PROGRAM_PARAM_HPP\n'
)

param_path = '/home/cepheid/TriForce/game_engine/params/program_params.hpp'
with open(param_path, 'w+') as f:
    cur_header = f.read()
    if cur_header != program_params:
        f.write(program_params)

# build_path = '/home/cepheid/TriForce/game_engine/buildDir'
# meson = ['meson', 'compile', '-C', build_path]
# subprocess.run(meson)
