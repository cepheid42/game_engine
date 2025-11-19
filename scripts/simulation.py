from pathlib import Path
from dataclasses import dataclass, field

from particles import *
from collisions import *
from electromagnetics import *


@dataclass
class Simulation:
    name: str
    shape: tuple
    nthreads: int
    nt: int
    dt: float
    t_end: float
    em_params: EMParams = field(default_factory=EMParams)
    particle_params: Particles = field(default_factory=Particles)
    cfl: float = 1.0
    x_range: tuple = ()
    y_range: tuple = ()
    z_range: tuple = ()
    deltas: tuple = ()
    em_enabled: bool = True
    push_enabled: bool = True
    jdep_enabled: bool = True
    coll_enabled: bool = True
    
    def __repr__(self):
        return (f'inline constexpr auto nThreads = {self.nthreads};\n'
                '\n'
                f'inline constexpr auto x_collapsed = {str(self.shape[0] <= 2).lower()};\n'
                f'inline constexpr auto y_collapsed = {str(self.shape[1] <= 2).lower()};\n'
                f'inline constexpr auto z_collapsed = {str(self.shape[2] <= 2).lower()};\n'
                '\n'
                f'inline constexpr auto Nx = {self.shape[0]}zu;\n'
                f'inline constexpr auto Ny = {self.shape[1]}zu;\n'
                f'inline constexpr auto Nz = {self.shape[2]}zu;\n'
                '\n'
                f'inline constexpr std::array x_range = {{{self.x_range[0]}, {self.x_range[1]}}};\n'
                f'inline constexpr std::array y_range = {{{self.y_range[0]}, {self.y_range[1]}}};\n'
                f'inline constexpr std::array z_range = {{{self.z_range[0]}, {self.z_range[1]}}};\n'
                '\n'
                f'inline constexpr auto dx = {self.deltas[0]};\n'
                f'inline constexpr auto dy = {self.deltas[0]};\n'
                f'inline constexpr auto dz = {self.deltas[0]};\n'
                '\n'
                f'inline constexpr auto cfl   = {self.cfl};\n'
                f'inline constexpr auto dt    = {self.dt};\n'
                f'inline constexpr auto t_end = {self.t_end};\n'
                f'inline constexpr auto Nt    = {self.nt}zu;\n'
                '\n'
                f'inline constexpr auto sim_name = "{self.name}";\n'
                # f'inline constexpr auto sim_path = "{project_path}";\n'
                '\n'
                f'inline constexpr auto   em_enabled = {str(self.em_enabled).lower()};\n'
                f'inline constexpr auto push_enabled = {str(self.push_enabled).lower()};\n'
                f'inline constexpr auto jdep_enabled = {str(self.jdep_enabled).lower()};\n'
                f'inline constexpr auto coll_enabled = {str(self.coll_enabled).lower()};\n')

