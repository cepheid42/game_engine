from dataclasses import dataclass

from .header_utils import *

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
class GeneralParams:
    name: str
    nthreads: int
    nt: int
    dt: float
    t_end: float
    domain: field(default_factory=)

    def __repr__(self):
        return str(
             constexpr_declaration('sim_name', f'"{self.name}"') + 
             constexpr_declaration('sim_path', f'"{self.sim_path}"') + 
             constexpr_declaration('data_path', f'"{self.data_path}"') + 
             constexpr_declaration('nThreads', f'{self.nthreads}lu')
        )
