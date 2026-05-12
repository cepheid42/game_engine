from dataclasses import dataclass
from enum import StrEnum


class ParticlePushType(StrEnum):
    Ballistic = 'ParticlePushType::Ballistic'
    Boris = 'ParticlePushType::Boris'
    HC = 'ParticlePushType::HigueraCary'

    def get_name(self):
        return str(self).split(':')[-1]


class ParticleBCType(StrEnum):
    Reflecting = 'ParticleBCType::Reflecting'
    Periodic = 'ParticleBCType::Periodic'
    Outflow = 'ParticleBCType::Outflow'


@dataclass
class ParticleParams:
    save_interval: int = 1
    push_type: ParticlePushType = ParticlePushType.Boris
    particle_bcs: ParticleBCType = ParticleBCType.Outflow
    bc_depth : int = 3
    interp_order: int = 2
    sort_frequency: int = 100
    particle_data: tuple = ()
    collisions: tuple = ()


@dataclass
class Particles:
    name: str
    mass: float
    atomic_number: int
    charge: float
    temp: tuple
    density: float
    ppc: tuple
    file_path: str = ''
    distribution: str = 'relativistic'
    tracer: bool = False
    px_range: tuple = ()
    py_range: tuple = ()
    pz_range: tuple = ()

    def __repr__(self):
        filestr = f'{self.file_path}/{self.name.lower()}.bp' if self.distribution != 'none' else ''
        return (
            '   ParticleGroupSpec{\n'
            f'      .name = "{self.name}",\n'
            f'      .filepath = "{filestr}",\n'
            f'      .mass = {self.mass},\n'
            f'      .charge = {float(self.charge)},\n'
            f'      .atomic_number = {self.atomic_number},\n'
            f'      .tracer = {str(self.tracer).lower()}\n'
            '   }'
        )
