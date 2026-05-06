from dataclasses import dataclass, field
from pathlib import Path

from .em import EMParams
from .particle import ParticleParams
from .metrics import Metric
from .params.header_utils import *

@dataclass
class SimFlags:
    em_enabled : bool = True
    push_enabled : bool = True
    jdep_enabled : bool = True
    metrics_enabled : bool = True
    collisions_enabled : bool = False
    apm_enabled : bool = False
    velocity_backstep_enabled: bool = False

    def __repr__(self):
        return str(
            section_label('Compilation Flags') +
            constexpr_declaration('  em_enabled', bool2str(self.em_enabled)) +
            constexpr_declaration('push_enabled', bool2str(self.push_enabled)) +
            constexpr_declaration('jdep_enabled', bool2str(self.jdep_enabled)) +
            constexpr_declaration('coll_enabled', bool2str(self.collisions_enabled)) +
            constexpr_declaration(' apm_enabled', bool2str(self.apm_enabled)) +
            constexpr_declaration('metrics_enabled', bool2str(self.metrics_enabled)) +
            constexpr_declaration('velocity_backstep_enabled', bool2str(self.velocity_backstep_enabled))
        )

@dataclass
class Simulation:
    name : str       # simulation name
    project_path : str # root path for tf code
    sim_path : str   # path to simulation working directory
    data_path : str  # path to input data directory
    nthreads : int   # number of threads for omp parallelization

    flags : SimFlags
    em_params: EMParams = field(default_factory=EMParams)
    particle_params: ParticleParams = field(default_factory=ParticleParams)
    metric_params: Metrics = field(default_factory=Metrics)

    def __post_init__(self):
        self.finalize_grid()

    def finalize_grid(self):
        for i, d in enumerate(('x', 'y', 'z')):
            dim = self.grid_params.__getattribute__(d)

            if dim is None:
                self.metrics.__setattr__(f'_{d}lo', 0)
                self.metrics.__setattr__(f'_{d}hi', 1)
                continue

            lo_bc = self.em_params.em_bcs[2 * i]
            new_nodes_lo = 0
            if lo_bc is EMBCType.Periodic:
                new_nodes_lo += self.em_params.nhalo
            if lo_bc is EMBCType.PML:
                new_nodes_lo += self.em_params.pml_depth
            newmin = dim.minval - dim.delta * new_nodes_lo

            hi_bc = self.em_params.em_bcs[2 * i + 1]
            new_nodes_hi = 0
            if hi_bc is EMBCType.Periodic:
                new_nodes_hi += self.em_params.nhalo
            if hi_bc is EMBCType.PML:
                new_nodes_hi += self.em_params.pml_depth
            newmax = dim.maxval + dim.delta * new_nodes_hi

            newsize = dim.size + new_nodes_lo + new_nodes_hi

            if newsize != dim.size:
                self.grid_params.__setattr__(d, GridDimension(newmin, newmax, delta=dim.delta))
                print(f'Finalized dimension {d} : min = {newmin}, max = {newmax}, delta = {self.grid_params.__getattribute__(d).delta}, size = {self.grid_params.__getattribute__(d).size}')
            
            self.metrics.__setattr__(f'_{d}lo', new_nodes_lo)
            self.metrics.__setattr__(f'_{d}hi', self.grid_params.__getattribute__(d).size - new_nodes_hi)

    def update_header(self):

        program_params = (
            start_header_guard()
            + '\n'
            + include_string('array')
            + include_string('string')
            + include_string('tuple')
            + '\n'
            + '#include "particle_spec.hpp"\n'
            + '\n'
            + str(self.general_params)
            + '\n'
            + str(self.flags)
            + '\n'
            + str(self.grid_params)
            + '\n'
            + str(self.time_params)
            + '\n'
            + str(self.intervals)
            + '\n'
            + str(self.em_params)
            + '\n'
            + str(self.particle_params)
            + '\n'
            + str(self.metrics)
            + '\n'
            + end_header_guard()
        )

        print(program_params)

        header_path = Path(self.general_params.project_path) / "params/program_params.hpp"
        with open(header_path, 'w+') as f:
            print('Updating header...', end=' ')

            cur_header = f.read()
            if cur_header != program_params:
                f.write(program_params)
        print('Done.')