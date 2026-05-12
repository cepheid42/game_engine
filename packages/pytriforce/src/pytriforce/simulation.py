import os
from dataclasses import dataclass, field
from pathlib import Path

from .em import EMParams
from .grid import Domain
from .particles import ParticleParams
from .metrics import Metrics
from .header_utils import *

@dataclass
class SimFlags:
    em_enabled      : bool = True
    push_enabled    : bool = True
    jdep_enabled    : bool = True
    metrics_enabled : bool = True
    velocity_backstep_enabled: bool = True

    apm_enabled        : bool = False
    collisions_enabled : bool = False
    applied_fields_only : bool = False
    ionization_test_enabled: bool = False

    def __repr__(self):
        return str(
            section_label('Compilation Flags') +
            constexpr_declaration('  em_enabled', bool2str(self.em_enabled)) +
            constexpr_declaration('push_enabled', bool2str(self.push_enabled)) +
            constexpr_declaration('jdep_enabled', bool2str(self.jdep_enabled)) +
            constexpr_declaration('coll_enabled', bool2str(self.collisions_enabled)) +
            constexpr_declaration(' apm_enabled', bool2str(self.apm_enabled)) +
            constexpr_declaration('metrics_enabled', bool2str(self.metrics_enabled)) +
            constexpr_declaration('velocity_backstep_enabled', bool2str(self.velocity_backstep_enabled)) +
            constexpr_declaration('ionization_test_enabled', bool2str(self.ionization_test_enabled))
        )

@dataclass
class Simulation:
    name : str
    project_path : str
    data_path : str
    nthreads : int
    nt: int
    dt: float
    t_end: float
    domain_params: Domain
    em_params: EMParams = field(default_factory=EMParams)
    particle_params: ParticleParams = field(default_factory=ParticleParams)
    metric_params: Metrics = field(default_factory=Metrics)
    compile_flags : SimFlags = field(default_factory=SimFlags)

    def __repr__(self):
        return str(
            constexpr_declaration('sim_name', self.name) +
            constexpr_declaration('sim_path', self.project_path) +
            '\n' +
            constexpr_declaration('dt', float(self.dt)) +
            constexpr_declaration('t_end', float(self.t_end)) +
            constexpr_declaration('Nt', f'{int(self.nt)}zu')
        )


def update_header(sim, project_path=None):
    print('Updating header...', end=' ')
    if project_path is None:
        project_path = os.environ['TRIFORCE_ROOT']

    project_path = Path(project_path)
    header_path = project_path / "params/program_params.hpp"
    #
    # # Check if various dimensions are collapsed
    # x_collapsed = nx == 2
    # y_collapsed = ny == 2
    # z_collapsed = nz == 2
    #
    # bc_str = f'{em_bcs[0]}zu, {em_bcs[1]}zu, {em_bcs[2]}zu, {em_bcs[3]}zu, {em_bcs[4]}zu, {em_bcs[5]}zu'
    #
    # dt_em = 0.99 / (constants.c * np.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2))
    # num_subcycles_em = max(1, np.ceil(params.dt / dt_em))
    #
    # for p in particles.particle_data:
    #     p.file_path = f'/data/{params.name}'
    #
    # particle_types = ',\n'.join([str(p) for p in particles.particle_data])
    # collision_types = ',\n'.join([str(c) for c in particles.collisions])
    #
    # if params.applied_fields_only:
    #     assert em_params.applied_fields != ''

    local_includes = [
        '#include "particle_spec.hpp"',
    ]

    global_includes = [
        '#include <array>',
        '#include <string_view>'
    ]

    program_params = (
        '#ifndef PROGRAM_PARAM_HPP\n'
        '#define PROGRAM_PARAM_HPP\n'
        '\n'
        f'{"\n".join(i for i in local_includes)}'
        '\n'
        f'{"\n".join(i for i in global_includes)}'
        '\n'
        f'{sim}'
        '\n'
        f'{sim.domain_params}'
        '\n'
        f'{sim.compile_flags}'
        '\n'
        f'{sim.em_params}'
        '\n'
        f'{sim.particle_params}'
        '\n'
        f'{sim.metric_params}'
        '\n'
        '#endif //PROGRAM_PARAM_HPP\n'
    )

    with open(header_path, 'w+') as f:
        f.write(program_params)

    # # Does comparing with the old version actually provide any benefit?
    # with open(header_path, 'w+') as f:
    #     cur_header = f.read()
    #     if cur_header != program_params:
    #         f.write(program_params)

    print('Done.')