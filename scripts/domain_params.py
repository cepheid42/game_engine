from pathlib import Path
from dataclasses import dataclass, field
from enum import StrEnum

@dataclass
class EMParams:
    save_interval: int = 1
    nhalo: int = 0
    pml_depth: int = 10
    pml_grade: float = 3.5
    pml_alpha_max: float = 0.2
    em_bcs: tuple = (2, 2, 2, 2, 2, 2)
    applied_fields: str = ''

@dataclass
class ParticleParams:
    save_interval: int = 1
    particle_bcs: str = 'outflow'
    interp_order: int = 2
    particle_data: tuple = ()
    collisions: tuple = ()

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
            f'      .group1 = "{self.groups[0]}",\n'
            f'      .group2 = "{self.groups[1]}",\n'
            f'      .channels = {{"{channels}"}},\n'
            f'      .step_interval = {self.step_interval},\n'
            f'      .probability_search_area = {1.0},\n'
            f'      .self_scatter = {str(self.self_scatter).lower()},\n'
            f'      {channel_spec + '\n'}'
            '   }'
        )

@dataclass
class Particles:
    name: str
    mass: float
    atomic_number: int
    charge: float
    temp: tuple
    density: float
    ppc: tuple
    distribution: str = 'relativistic'
    px_range: tuple = ()
    py_range: tuple = ()
    pz_range: tuple = ()

    def __repr__(self):
        filestr = f'/data/{self.name.lower()}.bp' if self.distribution != 'none' else ''
        return (
            '   ParticleGroupSpec{\n'
            f'      .name = "{self.name}",\n'
            f'      .filepath = "{filestr}",\n'
            f'      .mass = {self.mass},\n'
            f'      .charge = {float(self.charge)},\n'
            f'      .atomic_number = {self.atomic_number}\n'
            '   }'
        )


class MetricType(StrEnum):
    ParticleDump = 'MetricType::ParticleDump'
    ParticleDiagnostics = 'MetricType::ParticleDiag'
    ParticleEnergy = 'MetricType::ParticleEnergy'
    FieldDump = 'MetricType::FieldDump'
    FieldEnergy = 'MetricType::FieldEnergy'

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
    cfl: float = 1.0
    x_range: tuple = ()
    y_range: tuple = ()
    z_range: tuple = ()
    deltas: tuple = ()
    em_enabled: bool = True
    push_enabled: bool = True
    jdep_enabled: bool = True
    coll_enabled: bool = True


def update_header(params: Simulation, project_path: str, ionization_test_override: bool=False):
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
    header_path = project_path / "params/program_params.hpp"

    # Check if various dimensions are collapsed
    x_collapsed = nx == 2
    y_collapsed = ny == 2
    z_collapsed = nz == 2

    particle_bcs = {
        'reflecting': 'ParticleBCType::Reflecting',
        'periodic': 'ParticleBCType::Periodic',
        'outflow': 'ParticleBCType::Outflow'
    }

    bc_str = f'{em_bcs[0]}zu, {em_bcs[1]}zu, {em_bcs[2]}zu, {em_bcs[3]}zu, {em_bcs[4]}zu, {em_bcs[5]}zu'

    particle_types = ',\n'.join([str(p) for p in particles.particle_data])
    collision_types = ',\n'.join([str(c) for c in particles.collisions])
    ionization_test = "#define IONIZATION_TEST_OVERRIDE\n\n" if ionization_test_override else ""

    program_params = (
        '#ifndef PROGRAM_PARAM_HPP\n'
        '#define PROGRAM_PARAM_HPP\n'
        '\n'
        '#include "particle_spec.hpp"\n'
        '\n'
        '#include <array>\n'
        '\n'
        f'{ionization_test}'
        f'inline constexpr auto nThreads = {params.nthreads};\n'
        '\n'
        f'inline constexpr auto x_collapsed = {str(x_collapsed).lower()};\n'
        f'inline constexpr auto y_collapsed = {str(y_collapsed).lower()};\n'
        f'inline constexpr auto z_collapsed = {str(z_collapsed).lower()};\n'
        '\n'
        f'inline constexpr auto Nx = {nx}zu;\n'
        f'inline constexpr auto Ny = {ny}zu;\n'
        f'inline constexpr auto Nz = {nz}zu;\n'
        '\n'
        f'inline constexpr std::array x_range = {{{xmin}, {xmax}}};\n'
        f'inline constexpr std::array y_range = {{{ymin}, {ymax}}};\n'
        f'inline constexpr std::array z_range = {{{zmin}, {zmax}}};\n'
        '\n'
        f'inline constexpr auto dx = {dx};\n'
        f'inline constexpr auto dy = {dy};\n'
        f'inline constexpr auto dz = {dz};\n'
        '\n'
        f'inline constexpr auto cfl   = {params.cfl};\n'
        f'inline constexpr auto dt    = {params.dt};\n'
        f'inline constexpr auto t_end = {params.t_end};\n'
        f'inline constexpr auto Nt    = {params.nt}zu;\n'
        '\n'
        f'inline constexpr auto sim_name = "{params.name}";\n'
        f'inline constexpr auto sim_path = "{project_path}";\n'
        '\n'
        f'inline constexpr auto   em_enabled = {str(params.em_enabled).lower()};\n'
        f'inline constexpr auto push_enabled = {str(params.push_enabled).lower()};\n'
        f'inline constexpr auto jdep_enabled = {str(params.jdep_enabled).lower()};\n'
        f'inline constexpr auto coll_enabled = {str(params.coll_enabled).lower()};\n'
        '\n'
        '/*---------------------------------------------------------------/\n'
        '/-                        EM Parameters                         -/\n'
        '/---------------------------------------------------------------*/\n'
        'enum class EMFace { X, Y, Z };\n'
        'enum class EMSide { Lo, Hi };\n'
        '\n'
        f'inline constexpr auto em_save_interval = {em_params.save_interval}zu;\n'
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
        f'inline constexpr auto applied_fields_path = "{em_params.applied_fields}";\n'
        '\n'
        '/*---------------------------------------------------------------/\n'
        '/-                     Particle Parameters                      -/\n'
        '/---------------------------------------------------------------*/\n'
        'enum class ParticleBCType { Reflecting, Periodic, Outflow };\n'
        '\n'
        f'inline constexpr auto particle_save_interval = {particles.save_interval}zu;\n'
        f'inline constexpr auto interpolation_order = {particles.interp_order}zu;\n'
        '\n'
        f'inline constexpr auto PBCSelect = {particle_bcs[particles.particle_bcs]};\n'
        '\n'
        'inline constexpr std::array particle_spec = {\n'
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
    print('Done.')
