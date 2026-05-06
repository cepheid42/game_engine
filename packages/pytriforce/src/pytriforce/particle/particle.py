from dataclasses import dataclass
from enum import StrEnum

from ..params.header_utils import *
from ..geometry.geometry import ExtrudedShape


class ParticleDistributionType(StrEnum):
    NoInit = 'NoInit'
    NonRelativistic = 'NonRelativistic'
    Relativistic = 'Relativistic'


class ParticlePusherType(StrEnum, NameGetter):
    Ballistic = 'ParticlePusherType::Ballistic'
    Boris = 'ParticlePusherType::Boris'
    HigueraCary = 'ParticlePusherType::HigueraCary'


class ParticleBCType(StrEnum, NameGetter):
    Reflecting = 'ParticleBCType::Reflecting'
    Periodic = 'ParticleBCType::Periodic'
    Outflow = 'ParticleBCType::Outflow'


@dataclass
class ParticleGroup:
    name : str
    mass : float
    charge : float
    temp : tuple
    density : tuple
    ppc : tuple
    atomic_number : int = 0
    tracer_fraction : float = 0.0
    file_path : str = None
    distribution : ParticleDistributionType = ParticleDistributionType.Relativistic
    geometry: ExtrudedShape = None

    def __repr__(self):
        filestr = f'{self.file_path}/{self.name.lower()}.bp' if self.distribution != ParticleDistributionType.NoInit else ''
        return (
            '   ParticleGroupSpec{\n'
            f'      .name = "{self.name}",\n'
            f'      .filepath = "{filestr}",\n'
            f'      .mass = {self.mass},\n'
            f'      .charge = {float(self.charge)},\n'
            f'      .atomic_number = {self.atomic_number},\n'
            f'      .tracer = {bool2str(self.tracer_fraction == 1.0)}\n'
            '   }'
        )


@dataclass
class ParticleParams:
    push_type : ParticlePusherType = ParticlePusherType.Boris
    particle_bcs : ParticleBCType = ParticleBCType.Outflow
    bc_depth : int = 1
    sort_frequency : int = 100
    interp_order : int = 1
    particle_groups : tuple = ()
    collisions : tuple = ()
    geometry_file : str = ""

    def __repr__(self) -> str:
        groupstr = ',\n  '.join([str(g) for g in self.particle_groups])
        simple_particle_boundary = len(self.geometry_file) == 0
        return str(
            section_label('Particles Params') +
            enum_declaration('ParticlePusherType', ParticlePusherType) +
            enum_declaration('ParticleBCType', ParticleBCType) +
            '\n' +
            constexpr_declaration('ParticlePusherSelect', self.push_type) +
            '\n' +
            constexpr_declaration('PBCSelect', self.particle_bcs) +
            constexpr_declaration('PBCDepth', f'{self.bc_depth}lu') +
            '\n' +
            constexpr_declaration('sort_frequency', f'{self.sort_frequency}lu') +
            constexpr_declaration('interpolation_order', f'{self.interp_order}lu') +
            '\n' +
            constexpr_declaration('std::array', 'particle_spec', f'{{\n  {groupstr}\n}}') +
            '\n' +
            constexpr_declaration('geometry_file', f'"{self.geometry_file}"sv') +
            '\n' +
            constexpr_declaration('simple_particle_boundary', f'{bool2str(simple_particle_boundary)}')
        )
