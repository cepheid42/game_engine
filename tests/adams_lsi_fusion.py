#!/usr/bin/env python3

# import sys
# sys.path.append('/home/cepheid/TriForce/game_engine')

import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from adios2 import FileReader

from scripts.pyforce import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'rlsi_fusion'
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
t_end = 1.5e-13
nt = int(t_end / dt) + 1

save_interval = 10

# =====================
# ===== Particles =====
# =====================
px_range = (-11.0e-6, -9.0e-6) # meters
py_range = (ymin, ymax)
pz_range = (-10.0e-6, 10.0e-6)

ppc = (5, 1, 5)
density = 1.0e29 #m^-3
Tcold = 1
Thot = 10

deuterium = Particles(
    name='deuterium',
    mass=2.0141017778 * constants.atomic_mass,
    charge=1,
    atomic_number=1,
    temp=tuple(3 * [Thot]),
    density=density,
    ppc=ppc,
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# neutrons = Particles(
#     name='neutrons',
#     mass=constants.m_n,
#     charge=0,
#     atomic_number=0,
#     temp=tuple(3 * [0.0]), # eV
#     density=0.0, # m^-3,
#     ppc=tuple(3 * [0.0]),
#     distribution='none',
#     px_range=px_range,
#     py_range=py_range,
#     pz_range=pz_range
# )
#
# protons = Particles(
#     name='protons',
#     mass=constants.m_p,
#     charge=1,
#     atomic_number=1,
#     temp=tuple(3 * [0.0]),
#     density=0.0,
#     ppc=tuple(3 * [0.0]),
#     distribution='none',
#     px_range=px_range,
#     py_range=py_range,
#     pz_range=pz_range
# )
# helium3 = Particles(
#     name='helium3',
#     mass=3.016029 * constants.atomic_mass,
#     charge=2,
#     atomic_number=2,
#     temp=tuple(3 * [0.0]), # eV
#     density=0.0, # m^-3,
#     ppc=tuple(3 * [0.0]),
#     distribution='none',
#     px_range=px_range,
#     py_range=py_range,
#     pz_range=pz_range
# )
#
# tritium = Particles(
#     name='tritium',
#     mass=3.01605 * constants.atomic_mass,
#     charge=1,
#     atomic_number=1,
#     temp=tuple(3 * [0.0]), # eV
#     density=0.0, # m^-3,
#     ppc=tuple(3 * [0.0]),
#     distribution='none',
#     px_range=px_range,
#     py_range=py_range,
#     pz_range=pz_range
# )

# ==========================================
# ===== Collisions and Particle Params =====
# ==========================================
coll_interval = 1

particle_params = ParticleParams(
    save_interval=save_interval,
    particle_bcs=ParticleBCType.Outflow,
    push_type=ParticlePushType.Boris,
    interp_order=1,
    particle_data=(deuterium, ),#tritium, protons, neutrons, helium3),
    collisions=(
        Collision(
            groups=(deuterium, deuterium),
            channels=('coulomb', ),#'fusion'),
            coulomb=CoulombParams(0.0, 1.0), # set LnLambda = 0 to calculate it on the fly
            # fusion=(
            #     FusionParams(
            #         products=(neutrons, helium3),
            #         energy_gain=3.269e6,
            #         cross_section_file=project_path + '/tests/collision_tests/cross_section_data/DD_nHe3_BH_eV_m2.txt'
            #     ),
            #     FusionParams(
            #         products=(tritium, protons),
            #         energy_gain=4.03e6,
            #         cross_section_file=project_path + '/tests/collision_tests/cross_section_data/DD_pT_BH_eV_m2.txt'
            #     ),
            # ),
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
        MetricType.ParticleEnergy,
        MetricType.FieldEnergy,
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
    # create_particles(sim_params, (deuterium, tritium, protons, neutrons, helium3), data_path)
    create_particles(sim_params, (deuterium,), data_path)
    update_header(sim_params, project_path=project_path, data_path=data_path)

    compile_project(build_path, output=True)
    run_project(build_path + '/game_engine', output=True)

# # ===========================
# # ===== Post Processing =====
# # ===========================
xs = np.linspace(xmin, xmax, shape[0])
zs = np.linspace(zmin, zmax, shape[2])

# block = True
block = False

plot_density('Deuterium', 7500, data_path, xs, zs, block=block)
plot_temperature('Deuterium', 7500, data_path, xs, zs, block=block)

# plot_density('electrons', 7500, data_path, xs, zs, block=block)
# plot_temperature('electrons', 7500, data_path, xs, zs, block=block)
# plot_density('ions', 7500, data_path, xs, zs, block=block)
# plot_temperature('ions', 7500, data_path, xs, zs, block=block)
plot_field_energy(data_path, block=block)
plot_particle_energy(data_path)

# plot_total_neutron_yield(data_path, [0, nt, save_interval])
