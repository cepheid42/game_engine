#!/usr/bin/env python3

# import subprocess
import math

nx = 16
ny = 2
nz = 16

xmin, xmax = -0.03, 0.03
ymin, ymax = 0.0, 0.004
zmin, zmax = -0.03, 0.03

dx = (xmax - xmin) / (nx - 1)
dy = (ymax - ymin) / (ny - 1)
dz = (zmax - zmin) / (nz - 1)

cfl = 0.95 / math.sqrt(3)
dt = 1.0e-10
t_end = 5.275e-7
nt = int(t_end / dt) + 1

save_interval = 1
nthreads = 1
interp_order = 1

PMLDepth = 1
PMLGrade = 3.5
PMLAlphaMax = 0.2

program_params = (
    '#ifndef PROGRAM_PARAM_HPP\n'
    '#define PROGRAM_PARAM_HPP\n'
    '\n'
    '#include "compute_type.hpp"\n'
    '\n'
    '#include <array>\n'
    '\n'
    f'inline constexpr auto nThreads = {nthreads};\n'
    '\n'
    f'inline constexpr auto Nx = {nx}zu;\n'
    f'inline constexpr auto Ny = {ny}zu;\n'
    f'inline constexpr auto Nz = {nz}zu;\n'
    '\n'
    f'inline constexpr std::array x_range = {{{xmin}_fp, {xmax}_fp}};\n'
    f'inline constexpr std::array y_range = {{{ymin}_fp, {ymax}_fp}};\n'
    f'inline constexpr std::array z_range = {{{zmin}_fp, {zmax}_fp}};\n'
    '\n'
    f'inline constexpr auto dx = {dx}_fp;\n'
    f'inline constexpr auto dy = {dy}_fp;\n'
    f'inline constexpr auto dz = {dz}_fp;\n'
    '\n'
    f'inline constexpr auto cfl   = {cfl}_fp;\n'
    f'inline constexpr auto dt    = {dt}_fp;\n'
    f'inline constexpr auto t_end = {t_end}_fp;\n'
    f'inline constexpr auto Nt    = {nt}zu;\n'
    '\n'
    f'inline constexpr auto save_interval = {save_interval}zu;\n'
    f'inline constexpr auto interpolation_order = {interp_order}zu;\n'
    '\n'
    f'inline constexpr auto Ncx = Nx - 1zu;\n'
    f'inline constexpr auto Ncy = Ny - 1zu;\n'
    f'inline constexpr auto Ncz = Nz - 1zu;\n'
    '\n'
    f'inline constexpr auto PMLDepth    = {PMLDepth}zu;\n'
    f'inline constexpr auto PMLGrade    = {PMLGrade}_fp;\n'
    f'inline constexpr auto PMLAlphaMax = {PMLAlphaMax}_fp;\n'
    f'//inline constexpr auto PMLKappaMax = 1.0_fp;\n'
    '\n'
    '#endif //PROGRAM_PARAM_HPP\n'
)

param_path = '/home/cepheid/TriForce/game_engine/params/program_params.hpp'
write_new = False
with open(param_path, 'r') as f:
    cur_header = f.read()
    if cur_header != program_params:
        write_new = True

if write_new:
    with open(param_path, 'w+') as f:
        f.write(program_params)

