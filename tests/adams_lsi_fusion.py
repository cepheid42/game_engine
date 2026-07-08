#!/usr/bin/env python3

# import sys
# sys.path.append('/home/cepheid/TriForce/game_engine')

import numpy as np
from scipy import constants

from scripts.pyforce import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'rlsi'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'

shape = (551, 2, 151)

xmin, xmax = -15.5e-6, -4.5e-6
ymin, ymax = 0.0, 0.01
zmin, zmax = -15.5e-6, 15.5e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 2.0e-17
t_end = 0.5e-13
nt = int(t_end / dt) + 1

save_interval = 50

# =====================
# ===== Particles =====
# =====================
px_range = (-11.0e-6, -9.0e-6) # meters
py_range = (ymin, ymax)
pz_range = (-10.0e-6, 10.0e-6)

ppc = (4, 1, 4)
density = 1.0e29 #m^-3
temp_eV = 10

deuterium = Particles(
    name='deuterium',
    mass=2.0141017778 * constants.atomic_mass,
    charge=1,
    atomic_number=1,
    temp=tuple(3 * [temp_eV]),
    density=density,
    ppc=ppc,
    distribution='relativistic',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-1,
    atomic_number=0,
    temp=tuple(3 * [temp_eV]),
    density=density,
    ppc=ppc,
    distribution='relativistic',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

neutrons = Particles(
    name='neutrons',
    mass=constants.m_n,
    charge=0,
    atomic_number=0,
    temp=tuple(3 * [0.0]), # eV
    density=0.0, # m^-3,
    ppc=tuple(3 * [0.0]),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

helium3 = Particles(
    name='helium3',
    mass=3.016029322 * constants.atomic_mass,
    charge=2,
    atomic_number=2,
    temp=tuple(3 * [0.0]), # eV
    density=0.0, # m^-3,
    ppc=tuple(3 * [0.0]),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

protons = Particles(
    name='protons',
    mass=constants.m_p,
    charge=1,
    atomic_number=1,
    temp=tuple(3 * [0.0]),
    density=0.0,
    ppc=tuple(3 * [0.0]),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

tritium = Particles(
    name='tritium',
    mass=3.01605 * constants.atomic_mass,
    charge=1,
    atomic_number=1,
    temp=tuple(3 * [0.0]), # eV
    density=0.0, # m^-3,
    ppc=tuple(3 * [0.0]),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

photons = Particles(
    name='photons',
    mass=0.0,
    charge=0.0,
    atomic_number=0,
    temp=tuple(3 * [0.0]), # eV
    density=0.0, # m^-3,
    ppc=tuple(3 * [0.0]),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# ==========================================
# ===== Collisions and Particle Params =====
# ==========================================
coll_interval = 1
production_mult = 1.0e8
particle_params = ParticleParams(
    save_interval=save_interval,
    particle_bcs=ParticleBCType.Outflow,
    push_type=ParticlePushType.Boris,
    interp_order=1,
    particle_data=(deuterium, neutrons, helium3, tritium, protons, electrons, photons),
    collisions=(
        Collision(
            groups=(electrons, electrons),
            channels=('coulomb',),
            coulomb=CoulombParams(0.0, 1.0), # set LnLambda = 0 to calculate it on the fly
            self_scatter=True,
            step_interval=coll_interval
        ),
        Collision(
            groups=(electrons, deuterium),
            channels=('coulomb', 'radiation'),
            coulomb=CoulombParams(0.0, 1.0), # set LnLambda = 0 to calculate it on the fly
            radiation=RadiationParams(
                products=photons,
                reduce_electron_energy=True,
                production_multiplier=production_mult,
                # min_energy=1e3,
                # max_energy=1e10,
                # use_TFD=True
                cross_section_file=project_path + '/tests/cross_section_data/SB_G4_Z1_kdsdk_MeV_barns.csv',
                use_TFD=False
            ),
            self_scatter=False,
            step_interval=coll_interval
        ),
        Collision(
            groups=(deuterium, deuterium),
            channels=('coulomb', 'fusion'),
            coulomb=CoulombParams(0.0, 1.0), # set LnLambda = 0 to calculate it on the fly
            fusion=(
                FusionParams(
                    products=(neutrons, helium3),
                    energy_gain=3.269e6,
                    production_multiplier=production_mult,
                    cross_section_file=project_path + '/tests/cross_section_data/DD_nHe3_BH_eV_m2.txt'
                ),
                FusionParams(
                    products=(tritium, protons),
                    energy_gain=4.03e6,
                    production_multiplier=production_mult,
                    cross_section_file=project_path + '/tests/cross_section_data/DD_pT_BH_eV_m2.txt'
                ),
            ),
            self_scatter=True,
            step_interval=coll_interval
        ),
    )
)


# ==================================
# ===== Electromagnetic Params =====
# ==================================
em_params = EMParams(
    save_interval=save_interval,
    pml_depth=15,
    em_bcs=(1, 1, 2, 2, 1, 1),
    laser_spec=Laser(8.0e-7, -2.75e13, 2.5479e-6, 5.0e-6, 0.60454),
)

# ==========================
# ===== Metrics Params =====
# ==========================
metric_params = Metrics(
    data_path,
    (
        # MetricType.ParticleEnergy,
        # MetricType.FieldEnergy,
        # MetricType.FieldDump,
        # MetricType.ParticleDump,
        MetricType.ParticleDiagnostics,
    )
)

# ============================
# ===== Simulation Class =====
# ============================
sim_params = Simulation(
    name=sim_name,
    shape=shape,
    nthreads=48,
    dt=dt,
    t_end=t_end,
    nt=nt,
    x_range=(xmin, xmax),
    y_range=(ymin, ymax),
    z_range=(zmin, zmax),
    deltas=(dx, dy, dz),
    em_params=em_params,
    particle_params=particle_params,
    metric_params=metric_params,
    collisions_enabled=True,
    push_enabled=True,
    jdep_enabled=True,
    em_enabled=True,
    velocity_backstep_enabled=True,
)

# ===========================
# ===== Compile and Run =====
# ===========================
run = True
# run = False

if run:
    print(f'Setting up "{sim_name}"')
    create_data_dir(data_path)
    create_particles(sim_params, (deuterium, tritium, protons, neutrons, helium3, electrons, photons), data_path)
    update_header(sim_params, project_path=project_path, data_path=data_path)

    compile_project(build_path, output=True)
    run_project(build_path + '/game_engine', output=True)

# ===========================
# ===== Post Processing =====
# ===========================
xs = np.linspace(xmin, xmax, shape[0])
zs = np.linspace(zmin, zmax, shape[2])

block = True
# block = False

data_path = project_path + f'/data/{sim_name}_longrun'
plot_step = 15000

# plot_density(['neutrons', 'protons', 'helium3', 'tritium', 'electrons', 'deuterium', 'photons'], plot_step, data_path, xs, zs, block=block)
# plot_temperature(['neutrons', 'protons', 'helium3', 'tritium', 'electrons', 'deuterium'], plot_step, data_path, xs, zs, block=block)
# plot_total_particle_yield(data_path, ['protons', 'neutrons', 'helium3', 'tritium', 'electrons', 'photons'], [0, nt - 2 * save_interval, save_interval])

# plot_density(['photons'], plot_step, data_path, xs, zs, block=block)
# plot_temperature(['electrons', 'deuterium'], plot_step, data_path, xs, zs, block=block)

