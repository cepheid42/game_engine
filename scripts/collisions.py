from dataclasses import dataclass, field

from scripts.particles import Particles


@dataclass
class CoulombParams:
    coulomb_log: float = 10.0
    rate_mult: float = 1.0


@dataclass
class IonizationParams:
    products: tuple = ()
    ionization_energy: float = 0.0
    rate_multiplier: float = 1.0
    production_multiplier: float = 1.0
    rejection_multiplier: float = 1.0
    constant_cross_section: float = 0.0
    cross_section_file: str = ''


@dataclass
class FusionParams:
    products: tuple = ()
    energy_gain: float = 0.0
    rate_multiplier: float = 1.0
    production_multiplier: float = 1.0
    constant_cross_section: float = 0.0
    cross_section_file: str = ''


@dataclass
class RadiationParams:
    products: str = ''
    reduce_electron_energy: bool = False
    production_multiplier: float = 1.0
    cross_section_file: str = ''


@dataclass
class Collision:
    groups: tuple = ()
    channels: tuple = ()
    self_scatter: bool = False
    step_interval: int = 1
    coulomb: CoulombParams = field(default_factory=CoulombParams)
    ionization: IonizationParams = field(default_factory=IonizationParams)
    fusion: FusionParams = field(default_factory=FusionParams)
    radiation: RadiationParams = field(default_factory=RadiationParams)

    def __repr__(self):
        channels = ', '.join([t for t in self.channels])
        coulomb = ''
        ionization = ''
        fusion = ''
        radiation = ''
        if 'coulomb' in self.channels:
            coulomb = f'.coulomb = {{.coulomb_log = {float(self.coulomb.coulomb_log)}, .rate_multiplier = {self.coulomb.rate_mult}}},'

        if 'ionization' in self.channels:
            ionization = (
                '\n      .ionization = {\n'
                f'         .product1 = "{self.ionization.products[0]}",\n'
                f'         .product2 = "{self.ionization.products[1]}",\n'
                f'         .cross_section_file = "{self.ionization.cross_section_file}",\n'
                f'         .ionization_energy = {self.ionization.ionization_energy},\n'
                f'         .rate_multiplier = {self.ionization.rate_multiplier},\n'
                f'         .production_multiplier = {self.ionization.production_multiplier},\n'
                f'         .rejection_multiplier = {self.ionization.rejection_multiplier},\n'
                f'         .constant_cross_section = {self.ionization.constant_cross_section},\n'
                '      },'
            )

        if 'fusion' in self.channels:
            fusion = (
                '\n      .fusion = {\n'
                f'         .product1 = "{self.fusion.products[0]}",\n'
                f'         .product2 = "{self.fusion.products[1]}",\n'
                f'         .cross_section_file = "{self.fusion.cross_section_file}",\n'
                f'         .energy_gain = {self.fusion.energy_gain},\n'
                f'         .rate_multiplier = {self.fusion.rate_multiplier},\n'
                f'         .production_multiplier = {self.fusion.production_multiplier},\n'
                f'         .constant_cross_section = {self.fusion.constant_cross_section},\n'
                '      },'
            )

        if 'radiation' in self.channels:
            reduce_energy = str(self.radiation.reduce_electron_energy).lower()
            radiation = (
                '\n      .radiation = {\n'
                f'         .product1 = "{self.radiation.products}",\n'
                f'         .cross_section_file = "{self.radiation.cross_section_file}",\n'
                f'         .production_multiplier = {self.radiation.production_multiplier},\n'
                f'         .reduce_electron_energy = {reduce_energy},\n'
                '      },'
            )

        channel_spec = '\t\t'.join([coulomb, ionization, fusion, radiation]).lstrip()

        return (
            '   CollisionSpec{\n'
            f'      .group1 = "{self.groups[0].name}",\n'
            f'      .group2 = "{self.groups[1].name}",\n'
            f'      .channels = {{"{channels}"}},\n'
            f'      .step_interval = {self.step_interval},\n'
            f'      .probability_search_area = {1.0},\n'
            f'      .self_scatter = {str(self.self_scatter).lower()},\n'
            f'      {channel_spec + '\n'}'
            '   }'
        )