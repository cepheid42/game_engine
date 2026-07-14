#!/usr/bin/env python3

# import sys
# sys.path.append('/home/cepheid/TriForce/game_engine')

import numpy as np
from scipy import constants

from scripts.pyforce import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'lsi_smith_reduced'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'


shape = (1551, 2, 401)

xmin, xmax = -15.5e-6, 15.5e-6
ymin, ymax = 0.0, 0.01
zmin, zmax = -15.5e-6, 15.5e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 2.0e-17
t_end = 3.0e-13
nt = int(t_end / dt) + 1

save_interval = 150

# =====================
# ===== Particles =====
# =====================
px_range = (-5e-7, 5e-7) # meters
py_range = (ymin, ymax)
pz_range = (-1e-5, 1e-5)

ppc = (10, 1, 10)
density = 8.5e28 #m^-3
temp_eV = 10000

hydrogen = Particles(
    name='hydrogen',
    mass=1.008 * constants.atomic_mass,
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
    particle_data=(hydrogen, electrons),
    collisions=(
        Collision(
            groups=(electrons, electrons),
            channels=('coulomb',),
            coulomb=CoulombParams(0.0, 1.0), # set LnLambda = 0 to calculate it on the fly
            self_scatter=True,
            step_interval=coll_interval
        ),
        Collision(
            groups=(electrons, hydrogen),
            channels=('coulomb',),
            coulomb=CoulombParams(0.0, 1.0), # set LnLambda = 0 to calculate it on the fly
            self_scatter=False,
            step_interval=coll_interval
        ),
        Collision(
            groups=(hydrogen, hydrogen),
            channels=('coulomb',),
            coulomb=CoulombParams(0.0, 1.0), # set LnLambda = 0 to calculate it on the fly
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
    laser_spec=Laser(8.0e-7, -2.75e13, 2.5479e-6, 15.0e-6, 1.28855495),
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
    nthreads=64,
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
    create_particles(sim_params, (hydrogen, electrons), data_path)
    update_header(sim_params, project_path=project_path, data_path=data_path)

    compile_project(build_path, output=True)
    run_project(build_path + '/game_engine', output=True)

# ===========================
# ===== Post Processing =====
# ===========================
xs = np.linspace(xmin, xmax, shape[0])
zs = np.linspace(zmin, zmax, shape[2])

# data_path = project_path + f'/data/lsi_smith'

block = True
# block = False

plot_step = 7500

# plot_density(['electrons', 'hydrogen'], plot_step, data_path, xs, zs, block=block)
# plot_temperature(['electrons', 'hydrogen'], plot_step, data_path, xs, zs, block=block)

# plot_field_energy(data_path, block=block)
# plot_particle_energy(data_path, block=block)
