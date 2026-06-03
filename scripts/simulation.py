from pathlib import Path
from dataclasses import dataclass, field
from enum import StrEnum
from scipy import constants
import numpy as np
import shutil

from scripts.particles import ParticleParams


class MetricType(StrEnum):
    ParticleDump = 'MetricType::ParticleDump'
    ParticleDiagnostics = 'MetricType::ParticleDiag'
    ParticleEnergy = 'MetricType::ParticleEnergy'
    FieldDump = 'MetricType::FieldDump'
    FieldEnergy = 'MetricType::FieldEnergy'


@dataclass
class EMParams:
    save_interval: int = 10
    nhalo: int = 0
    pml_depth: int = 10
    pml_grade: float = 3.5
    pml_alpha_max: float = 0.2
    em_bcs: tuple = (2, 2, 2, 2, 2, 2)
    applied_fields: str = ''
    laser_enabled: bool = False # todo: these should be specs instead


@dataclass
class Metrics:
    data_path: str = ''
    metrics: tuple = ()

    def __repr__(self):
        return ',\n\t'.join(m for m in self.metrics)


@dataclass
class Simulation:
    name: str
    shape: tuple
    nthreads: int
    nt: int
    dt: float
    t_end: float
    em_params: EMParams = field(default_factory=EMParams)
    particle_params: ParticleParams = field(default_factory=ParticleParams)
    metric_params: Metrics = field(default_factory=Metrics)
    x_range: tuple = ()
    y_range: tuple = ()
    z_range: tuple = ()
    deltas: tuple = ()
    em_enabled: bool = True
    push_enabled: bool = True
    jdep_enabled: bool = True
    collisions_enabled: bool = True
    velocity_backstep_enabled: bool = True
    applied_fields_only: bool = False
    ionization_test_enabled: bool = False


def update_header(params, project_path, data_path):
    print('Updating header...', end=' ')
    nx, ny, nz = params.shape
    xmin, xmax = params.x_range
    ymin, ymax = params.y_range
    zmin, zmax = params.z_range
    dx, dy, dz = params.deltas
    em_params = params.em_params
    em_bcs = em_params.em_bcs
    particles = params.particle_params
    metrics = params.metric_params

    project_path = Path(project_path)
    data_path = Path(data_path)
    header_path = project_path / "params/program_params.hpp"

    # Check if various dimensions are collapsed
    x_collapsed = nx == 2
    y_collapsed = ny == 2
    z_collapsed = nz == 2

    bc_str = f'{em_bcs[0]}zu, {em_bcs[1]}zu, {em_bcs[2]}zu, {em_bcs[3]}zu, {em_bcs[4]}zu, {em_bcs[5]}zu'

    dt_em = 0.99 / (constants.c * np.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2))
    num_subcycles_em = max(1, np.ceil(params.dt / dt_em))
    if num_subcycles_em == 1:
        dt_em = params.dt

    for p in particles.particle_data:
        p.file_path = f'/data/{params.name}'

    particle_types = ',\n'.join([str(p) for p in particles.particle_data])
    collision_types = ',\n'.join([str(c) for c in particles.collisions])

    if params.applied_fields_only:
        assert em_params.applied_fields != ''

    program_params = (
        '#ifndef PROGRAM_PARAM_HPP\n'
        '#define PROGRAM_PARAM_HPP\n'
        '\n'
        '#include "particle_spec.hpp"\n'
        '\n'
        '#include <array>\n'
        '\n'
        f'inline constexpr auto nThreads = {params.nthreads};\n'
        '\n'
        f'inline constexpr auto x_collapsed = {str(x_collapsed).lower()};\n'
        f'inline constexpr auto y_collapsed = {str(y_collapsed).lower()};\n'
        f'inline constexpr auto z_collapsed = {str(z_collapsed).lower()};\n'
        '\n'
        f'inline constexpr auto Nx = {int(nx)}zu;\n'
        f'inline constexpr auto Ny = {int(ny)}zu;\n'
        f'inline constexpr auto Nz = {int(nz)}zu;\n'
        '\n'
        f'inline constexpr std::array x_range = {{{float(xmin)}, {float(xmax)}}};\n'
        f'inline constexpr std::array y_range = {{{float(ymin)}, {float(ymax)}}};\n'
        f'inline constexpr std::array z_range = {{{float(zmin)}, {float(zmax)}}};\n'
        '\n'
        f'inline constexpr auto dx = {float(dx)};\n'
        f'inline constexpr auto dy = {float(dy)};\n'
        f'inline constexpr auto dz = {float(dz)};\n'
        '\n'
        f'inline constexpr auto dt    = {float(params.dt)};\n'
        f'inline constexpr auto t_end = {float(params.t_end)};\n'
        f'inline constexpr auto Nt    = {int(params.nt)}zu;\n'
        '\n'
        f'inline constexpr auto sim_name = "{params.name}";\n'
        f'inline constexpr auto sim_path = "{project_path}";\n'
        '\n'
        f'inline constexpr auto   em_enabled = {str(params.em_enabled).lower()};\n'
        f'inline constexpr auto push_enabled = {str(params.push_enabled).lower()};\n'
        f'inline constexpr auto jdep_enabled = {str(params.jdep_enabled).lower()};\n'
        f'inline constexpr auto coll_enabled = {str(params.collisions_enabled).lower()};\n'
        f'inline constexpr auto applied_fields_only = {str(params.applied_fields_only).lower()};\n'
        f'inline constexpr auto velocity_backstep_enabled = {str(params.velocity_backstep_enabled).lower()};\n'
        f'inline constexpr auto ionization_test_enabled = {str(params.ionization_test_enabled).lower()};\n'
        '\n'
        '/*---------------------------------------------------------------/\n'
        '/-                        EM Parameters                         -/\n'
        '/---------------------------------------------------------------*/\n'
        'enum class EMFace { X, Y, Z };\n'
        'enum class EMSide { Lo, Hi };\n'
        '\n'
        f'inline constexpr auto em_save_interval = {int(em_params.save_interval)}zu;\n'
        f'inline constexpr auto em_subcycles = {int(num_subcycles_em)}zu;\n'
        f'inline constexpr auto dt_em = {float(dt_em)};\n'
        '\n'
        f'inline constexpr auto PMLDepth    = {em_params.pml_depth}zu;\n'
        f'inline constexpr auto PMLGrade    = {em_params.pml_grade};\n'
        f'inline constexpr auto PMLAlphaMax = {em_params.pml_alpha_max};\n'
        '//inline constexpr auto PMLKappaMax = 1.0;\n'
        '\n'
        f'inline constexpr auto nHalo = {em_params.nhalo}zu;\n'
        '\n'
        '// Periodic = 0, PML = 1, Reflecting = 2\n'
        f'inline constexpr std::array BCSelect = {{{bc_str}}};\n'
        '\n'
        f'inline constexpr auto laser_enabled = {str(em_params.laser_enabled).lower()};\n'
        f'inline constexpr auto applied_fields_path = "{em_params.applied_fields}";\n'
        '\n'
        '/*---------------------------------------------------------------/\n'
        '/-                     Particle Parameters                      -/\n'
        '/---------------------------------------------------------------*/\n'
        'enum class ParticleBCType { Reflecting, Periodic, Outflow };\n'
        'enum class ParticlePushType { Ballistic, Boris, HigueraCary };\n'
        '\n'
        f'inline constexpr auto particle_save_interval = {particles.save_interval}zu;\n'
        f'inline constexpr auto sort_frequency = {particles.sort_frequency}zu;\n'
        f'inline constexpr auto interpolation_order = {particles.interp_order}zu;\n'
        f'inline constexpr auto ParticlePushSelect = {str(particles.push_type)};\n'
        f'inline constexpr auto PBCSelect = {str(particles.particle_bcs)};\n'
        f'inline constexpr auto PBCDepth = {particles.bc_depth}zu;\n'
        '\n'
        f'inline constexpr std::array<ParticleGroupSpec, {len(particles.particle_data)}> particle_spec = {{\n'
        f'{particle_types}\n'
        '};\n'
        '\n'
        f'inline constexpr std::array<CollisionSpec, {len(particles.collisions)}> collision_spec = {{\n'
        f'{collision_types}\n'
        '};\n'
        '\n'
        '/*---------------------------------------------------------------/\n'
        '/-                      Metrics Parameters                      -/\n'
        '/---------------------------------------------------------------*/\n'
        'enum class MetricType { ParticleDump, ParticleDiag, ParticleEnergy, FieldDump, FieldEnergy };\n'
        '\n'
        f'inline constexpr auto metric_data_path = "{metrics.data_path}";\n'
        f'inline constexpr std::array<MetricType, {len(metrics.metrics)}> metric_spec = {{\n'
        f'\t{metrics}\n'
        '};\n'
        '\n'
        '#endif //PROGRAM_PARAM_HPP\n'
    )

    with open(header_path, 'w+') as f:
        cur_header = f.read()
        if cur_header != program_params:
            f.write(program_params)

    shutil.copy(header_path, data_path)
    print('Done.')
