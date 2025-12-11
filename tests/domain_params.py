from pathlib import Path
from dataclasses import dataclass, field

@dataclass
class EMParams:
    save_interval: int = 1
    nhalo: int = 0
    pml_depth: int = 10
    pml_grade: float = 3.5
    pml_alpha_max: float = 0.2
    em_bcs: tuple = (2, 2, 2, 2, 2, 2)

@dataclass
class ParticleParams:
    save_interval: int = 1
    particle_bcs: str = 'outflow'
    interp_order: int = 2
    beam_source: str = ''
    particle_data: tuple = ()
    collisions: tuple = ()

@dataclass
class Collision:
    groups: tuple = ()
    products: tuple = ()
    types: tuple = ()
    coulomb_log: float = -1.0
    rate_mult: float = 1.0
    self_scatter: bool = False
    step_interval: int = 1
    ionization_energy: float = 0.0
    cross_section_file: str = ''

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
    init: str = 'file'

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
    cfl: float = 1.0
    x_range: tuple = ()
    y_range: tuple = ()
    z_range: tuple = ()
    deltas: tuple = ()
    em_enabled: bool = True
    push_enabled: bool = True
    jdep_enabled: bool = True
    coll_enabled: bool = True


def update_header(params: Simulation, project_path: str):
    print('Updating header...', end=' ')
    nx, ny, nz = params.shape
    xmin, xmax = params.x_range
    ymin, ymax = params.y_range
    zmin, zmax = params.z_range
    dx, dy, dz = params.deltas
    em_params = params.em_params
    em_bcs = em_params.em_bcs
    particles = params.particle_params

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
    particle_data = ', '.join(['"/data/' + p + ('.bp"' if i == 'file' else '.empty"') for p, i in particles.particle_data])

    collision_params = params.particle_params.collisions
    collision_str = []
    for c in collision_params:
        result = f'   std::tuple("{c.groups[0]}", "{c.groups[1]}", "{c.products[0]}", "{c.products[1]}", {c.coulomb_log}, {c.rate_mult}, {c.step_interval}, {str(c.self_scatter).lower()}, {c.ionization_energy}, "{c.cross_section_file}"),'
        collision_str.append(result)

    collisions = ''
    if collision_str:
        collisions = '\n' + '\n'.join(c for c in collision_str) + '\n'

    program_params = (
        '#ifndef PROGRAM_PARAM_HPP\n'
        '#define PROGRAM_PARAM_HPP\n'
        '\n'
        '#include <array>\n'
        '#include <string>\n'
        '#include <tuple>\n'
        '\n'
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
        '/*---------------------------------------------------------------/\n'
        '/-                     Particle Parameters                      -/\n'
        '/---------------------------------------------------------------*/\n'
        'using collision_spec = std::tuple<std::string, std::string, std::string, std::string, double, double, int, bool, double, std::string>;\n'
        'enum class ParticleBCType { Static, Reflecting, Periodic, Outflow };\n'
        '\n'
        f'inline constexpr auto particle_save_interval = {particles.save_interval}zu;\n'
        f'inline constexpr auto interpolation_order = {particles.interp_order}zu;\n'
        '\n'
        f'inline constexpr auto PBCSelect = {particle_bcs[particles.particle_bcs]};\n'
        '\n'
        f'inline constexpr std::array particle_data = {{{particle_data}}};\n'
        f'inline constexpr auto particle_beam_file = "{particles.beam_source}";'
        '\n'
        f'inline constexpr std::array<collision_spec, {len(collision_str)}> collision_params = {{'
        f'{collisions}'
        '};\n'
        '\n'
        '#endif //PROGRAM_PARAM_HPP\n'
    )

    with open(header_path, 'w+') as f:
        cur_header = f.read()
        if cur_header != program_params:
            f.write(program_params)
    print('Done.')
