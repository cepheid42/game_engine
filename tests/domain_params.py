# import numpy as np
# from scipy import constants
import os
from pathlib import Path
from dataclasses import dataclass

@dataclass
class Simulation:
    name: str
    shape: tuple
    save_interval: int
    nthreads: int
    particle_bcs: int
    interp_order: int
    em_bcs: tuple
    pml_depth: int
    dt: float
    t_end: float
    nt: int
    pml_grade: float = 3.5
    pml_alpha_max: float = 0.2
    nhalo: int = 0
    cfl: float = 1.0
    x_range: tuple = ()
    y_range: tuple = ()
    z_range: tuple = ()
    deltas: tuple = ()
    particle_data: tuple = ()


def update_header(params, project_path):
    assert params.dt != 1.0
    print('Updating header...', end=' ')
    nx, ny, nz = params.shape
    xmin, xmax = params.x_range
    ymin, ymax = params.y_range
    zmin, zmax = params.z_range
    dx, dy, dz = params.deltas
    em_bcs = params.em_bcs

    project_path = Path(project_path)

    # Get TriForce project root directory
    # project_path = Path(os.getenv('TRIFORCE_ROOT', ''))
    # if not project_path.exists():
    #     raise EnvironmentError('Could not find TRIFORCE_ROOT env variable.')

    header_path = project_path / "params/program_params.hpp"

    # Check if 2D in Y direction
    is_2d_xz = 'true' if ny == 2 else 'false'

    bc_str = f'{em_bcs[0]}lu, {em_bcs[1]}lu, {em_bcs[2]}lu, {em_bcs[3]}lu, {em_bcs[4]}lu, {em_bcs[5]}lu'
    particle_bc = 'ParticleBCType::Periodic' if params.particle_bcs == 0 else 'ParticleBCType::Outflow'
    particle_data = ', '.join(['"/data/' + p + '.bp"' for p in params.particle_data])

    program_params = (
        '#ifndef PROGRAM_PARAM_HPP\n'
        '#define PROGRAM_PARAM_HPP\n'
        '\n'
        '#include <array>\n'
        '\n'
        f'inline constexpr auto nThreads = {params.nthreads};\n'
        '\n'
        f'inline constexpr auto is_2D_XZ = {is_2d_xz};\n'
        '\n'
        f'inline constexpr auto Nx = {nx}lu;\n'
        f'inline constexpr auto Ny = {ny}lu;\n'
        f'inline constexpr auto Nz = {nz}lu;\n'
        '\n'
        f'inline constexpr std::array x_range = {{{xmin}, {xmax}}};\n'
        f'inline constexpr std::array y_range = {{{ymin}, {ymax}}};\n'
        f'inline constexpr std::array z_range = {{{zmin}, {zmax}}};\n'
        '\n'
        f'inline constexpr auto dx = {dx};\n'
        f'inline constexpr auto dy = {dy};\n'
        f'inline constexpr auto dz = {dz};\n'
        '\n'
        f'inline constexpr auto cfl   = {params.cfl};\n'
        f'inline constexpr auto dt    = {params.dt};\n'
        f'inline constexpr auto t_end = {params.t_end};\n'
        f'inline constexpr auto Nt    = {params.nt}lu;\n'
        '\n'
        f'inline constexpr auto save_interval = {params.save_interval}lu;\n'
        '\n'
        f'inline constexpr auto sim_name = "{params.name}";\n'
        f'inline constexpr auto sim_path = "{project_path}";\n'
        '\n'
        '/*---------------------------------------------------------------/\n'
        '/-                        EM Parameters                         -/\n'
        '/---------------------------------------------------------------*/\n'
        'enum class EMFace { X, Y, Z };\n'
        'enum class EMSide { Lo, Hi };\n'
        '\n'
        f'inline constexpr auto PMLDepth    = {params.pml_depth}lu;\n'
        f'inline constexpr auto PMLGrade    = {params.pml_grade};\n'
        f'inline constexpr auto PMLAlphaMax = {params.pml_alpha_max};\n'
        '//inline constexpr auto PMLKappaMax = 1.0;\n'
        '\n'
        f'inline constexpr auto nHalo = 2lu;\n'
        '\n'
        '// Periodic = 0, PML = 1, Reflecting = 2\n'
        f'inline constexpr std::array BCSelect = {{{bc_str}}};\n'
        '\n'
        '/*---------------------------------------------------------------/\n'
        '/-                     Particle Parameters                      -/\n'
        '/---------------------------------------------------------------*/\n'
        'enum class ParticleBCType { Periodic, Outflow };\n'
        '\n'
        f'inline constexpr auto interpolation_order = {params.interp_order}lu;\n'
        '\n'
        f'inline constexpr auto PBCSelect = {particle_bc};\n'
        '\n'
        f'inline constexpr std::array particle_data = {{{particle_data}}};\n'
        '#endif //PROGRAM_PARAM_HPP\n'
    )

    with open(header_path, 'w+') as f:
        cur_header = f.read()
        if cur_header != program_params:
            f.write(program_params)
    print('Done.')
