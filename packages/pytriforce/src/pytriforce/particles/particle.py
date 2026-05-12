from dataclasses import dataclass
from enum import StrEnum

from ..header_utils import *
# from ..geometry.geometry import ExtrudedShape


class ParticleDistributionType(StrEnum):
    NoInit = 'NoInit'
    NonRelativistic = 'NonRelativistic'
    Relativistic = 'Relativistic'
    SP_UniformE = 'SP_UniformE'
    SP_UniformB = 'SP_UniformB'
    SP_ForceFree = 'SP_ForceFree'
    SP_PerpFields = 'SP_PerpFields'


class ParticlePusherType(NameGetter, StrEnum):
    Ballistic = 'ParticlePusherType::Ballistic'
    Boris = 'ParticlePusherType::Boris'
    HigueraCary = 'ParticlePusherType::HigueraCary'


# todo: add velocity and position enums or whatever



class ParticleBCType(NameGetter, StrEnum):
    Reflecting = 'ParticleBCType::Reflecting'
    Periodic = 'ParticleBCType::Periodic'
    Outflow = 'ParticleBCType::Outflow'


@dataclass
class ParticleGroup:
    name : str
    file_path : str
    mass : float
    charge : float
    density : float
    temp : tuple
    ppc : tuple
    atomic_number : int = 0
    tracer_fraction : float = 0.0
    # geometry: ExtrudedShape = None
    distribution : ParticleDistributionType = ParticleDistributionType.Relativistic

    def __repr__(self):
        filestr = ''
        if self.distribution != ParticleDistributionType.NoInit:
            filestr = f'{self.file_path}/{self.name.lower()}.bp'

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
    save_interval : int = 100
    bc_depth : int = 1
    sort_frequency : int = 100
    interp_order : int = 1
    particle_groups : tuple = ()
    collisions : tuple = ()
    geometry_file : str = ''

    def __repr__(self):
        groupstr = ',\n  '.join([str(g) for g in self.particle_groups])
        return str(
            section_label('Particles Params') +
            enum_declaration('ParticlePusherType', ParticlePusherType) +
            enum_declaration('ParticleBCType', ParticleBCType) +
            '\n' +
            constexpr_declaration('ParticlePusherSelect', self.push_type) +
            constexpr_declaration('PBCSelect', self.particle_bcs) +
            constexpr_declaration('PBCDepth', f'{self.bc_depth}lu') +
            '\n' +
            constexpr_declaration('sort_frequency', f'{self.sort_frequency}lu') +
            constexpr_declaration('interpolation_order', f'{self.interp_order}lu') +
            '\n' +
            constexpr_declaration('geometry_file', f'"{self.geometry_file}"sv') +
            '\n' +
            constexpr_declaration('simple_particle_boundary', f'{bool2str(self.geometry_file == '')}') +
            '\n' +
            constexpr_declaration('std::array', 'particle_spec', f'{{\n  {groupstr}\n}}')
        )
