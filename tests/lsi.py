#!/usr/bin/env python3

# import subprocess
import math
from scipy import constants

sim_name = 'lsi'

nx = 1501 # 751
ny = 2
nz = 1501 # 401
nhalo = 0

xmin, xmax = -15.0e-6, 15.0e-6
ymin, ymax = 0.0, 0.01
zmin, zmax = -15.0e-6, 15.0e-6

dx = (xmax - xmin) / (nx - 1)
dy = (ymax - ymin) / (ny - 1)
dz = (zmax - zmin) / (nz - 1)

dt = 4.0e-17
t_end = 3.0e-13
nt = int(t_end / dt) + 1
cfl = constants.c * dt * math.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

save_interval = 100
nthreads = 32
interp_order = 2

PMLDepth = 25
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
