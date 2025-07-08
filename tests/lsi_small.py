#!/usr/bin/env python3

# import subprocess
import math

sim_name = 'lsi'

nx = 101
ny = 2
nz = 101

xmin, xmax = -1.0e-6, 1.0e-6
ymin, ymax = 0.0, 2.0e-8
zmin, zmax = -1.0e-6, 1.0e-6

dx = (xmax - xmin) / (nx - 1)
dy = (ymax - ymin) / (ny - 1)
dz = (zmax - zmin) / (nz - 1)

cfl = 0.848 / math.sqrt(3)
dt = 4.0e-17 # 0.04 fs
t_end = 2.0e-14 # 300 fs
nt = int(t_end / dt) + 1 # ~7500

save_interval = 5
nthreads = 1
interp_order = 1

PMLDepth = 10
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
    f'inline constexpr auto Ncx = {nx - 1}zu;\n'
    f'inline constexpr auto Ncy = {ny - 1}zu;\n'
    f'inline constexpr auto Ncz = {nz - 1}zu;\n'
    '\n'
    f'inline constexpr auto PMLDepth    = {PMLDepth}zu;\n'
    f'inline constexpr auto PMLGrade    = {PMLGrade}_fp;\n'
    f'inline constexpr auto PMLAlphaMax = {PMLAlphaMax}_fp;\n'
    f'//inline constexpr auto PMLKappaMax = 1.0_fp;\n'
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
