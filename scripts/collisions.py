from dataclasses import dataclass, field

from scripts.particles import Particles


@dataclass
class CoulombParams:
    coulomb_log: float = 10.0
    rate_mult: float = 1.0

    def __repr__(self):
        return f'CoulombSpec{{.coulomb_log = {float(self.coulomb_log)}, .rate_multiplier = {self.rate_mult}}}'


@dataclass
class IonizationParams:
    products: tuple = ()
    ionization_energy: float = 0.0
    rate_multiplier: float = 1.0
    production_multiplier: float = 1.0
    rejection_multiplier: float = 1.0
    constant_cross_section: float = 0.0
    cross_section_file: str = ''
    
    def __repr__(self):
        return (
            'IonizationSpec{\n'
            f'         .product1 = "{self.products[0].name}",\n'
            f'         .product2 = "{self.products[1].name}",\n'
            f'         .cross_section_file = "{self.cross_section_file}",\n'
            f'         .ionization_energy = {self.ionization_energy},\n'
            f'         .rate_multiplier = {self.rate_multiplier},\n'
            f'         .production_multiplier = {self.production_multiplier},\n'
            f'         .rejection_multiplier = {self.rejection_multiplier},\n'
            f'         .constant_cross_section = {self.constant_cross_section}\n'
            '      },'
        )


@dataclass
class FusionParams:
    products: tuple = ()
    energy_gain: float = 0.0
    rate_multiplier: float = 1.0
    production_multiplier: float = 1.0
    constant_cross_section: float = 0.0
    cross_section_file: str = ''

    def __repr__(self):
        return (
            '         FusionSpec{\n'
            f'            .product1 = "{self.products[0].name}",\n'
            f'            .product2 = "{self.products[1].name}",\n'
            f'            .cross_section_file = "{self.cross_section_file}",\n'
            f'            .energy_gain = {self.energy_gain},\n'
            f'            .rate_multiplier = {self.rate_multiplier},\n'
            f'            .production_multiplier = {self.production_multiplier},\n'
            f'            .constant_cross_section = {self.constant_cross_section}\n'
            '         },\n'
        )


@dataclass
class RadiationParams:
    products: Particles = None
    reduce_electron_energy: bool = False
    production_multiplier: float = 1.0
    min_energy: float = 0.0
    max_energy: float = 0.0
    cross_section_file: str = ''
    use_TFD : bool = False

    def __repr__(self):
        if self.use_TFD:
            assert self.min_energy != 0.0 and self.max_energy != 0.0
        else:
            assert self.cross_section_file != ''

        reduce_energy = str(self.reduce_electron_energy).lower()
        return (
            'RadiationSpec{\n'
            f'         .product1 = "{self.products.name}",\n'
            f'         .cross_section_file = "{self.cross_section_file}",\n'
            f'         .production_multiplier = {self.production_multiplier},\n'
            f'         .min_energy = {self.min_energy},\n'
            f'         .max_energy = {self.max_energy},\n'
            f'         .reduce_electron_energy = {reduce_energy},\n'
            f'         .use_TFD = {str(self.use_TFD).lower()}\n'
            '      },'
        )


@dataclass
class Collision:
    groups: tuple = ()
    channels: tuple = ()
    self_scatter: bool = False
    step_interval: int = 1
    coulomb: CoulombParams = field(default_factory=CoulombParams)
    ionization: IonizationParams = field(default_factory=IonizationParams)
    fusion: tuple = ()
    radiation: RadiationParams = field(default_factory=RadiationParams)
    # inverse_radiation: InverseRadiationParams = field(default_factory=InverseRadiationParams)

    def __repr__(self):
        channels = ''
        for c in self.channels:
            channels += f'"{c}"'
            if c != self.channels[-1]:
                channels += ', '
        coulomb = ''
        ionization = ''
        fusion = ''
        radiation = ''
        # inv_radiation = ''
        if 'coulomb' in self.channels:
            coulomb = f'.coulomb = {self.coulomb},\n'

        if 'ionization' in self.channels:
            ionization = f'.ionization = {self.ionization}\n'

        if 'fusion' in self.channels:
            fusion = '.fusion = {\n'
            for f in self.fusion:
                fusion += str(f)
            fusion += '      },\n'

        if 'radiation' in self.channels:
            radiation = f'.radiation = {self.radiation}\n'

        # if 'inverse_radiation' in self.channels:
        #     inv_radiation = f'.inverse_radiation = {self.inverse_radiation}\n'

        channel_spec = '\t'.join([coulomb, ionization, fusion, radiation]).lstrip()

        return (
            '   CollisionSpec{\n'
            f'      .group1 = "{self.groups[0].name}",\n'
            f'      .group2 = "{self.groups[1].name}",\n'
            f'      .channels = {{{channels}}},\n'
            f'      .step_interval = {self.step_interval},\n'
            f'      .probability_search_area = {1.0},\n'
            f'      .self_scatter = {str(self.self_scatter).lower()},\n'
            f'      {channel_spec}\n'
            '   }'
        )