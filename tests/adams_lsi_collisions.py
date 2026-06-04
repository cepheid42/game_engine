#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from adios2 import FileReader

from scripts.plotting import plot_temperature
from scripts.pyforce import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'rlsi_allocator'
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

electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-1,
    atomic_number=0,
    temp=tuple(3 * [Thot]),
    density=density,
    ppc=ppc,
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

ions = Particles(
    name='ions',
    mass=constants.m_p,
    charge=1,
    atomic_number=1,
    temp=tuple(3 * [Tcold]),
    density=density,
    ppc=ppc,
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# ==========================================
# ===== Collisions and Particle Params =====
# ==========================================
coll_interval = 1
coulombParam = CoulombParams(0.0, 1) # set LnLambda = 0 to calculate it on the fly

particle_params = ParticleParams(
    save_interval=save_interval,
    particle_bcs=ParticleBCType.Outflow,
    push_type=ParticlePushType.Boris,
    interp_order=1,
    particle_data=(electrons, ions),
    collisions=(
        Collision(groups=(electrons, ions), channels=('coulomb',), step_interval=coll_interval),
        Collision(groups=(ions, ions), channels=('coulomb',), self_scatter=True, step_interval=coll_interval),
        Collision(groups=(electrons, electrons), channels=('coulomb',), self_scatter=True, step_interval=coll_interval),
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

# EM: 00:00:11.422479796
# Push: 00:00:15.103554543
# Jdep: 00:00:09.096954454
# IO: 00:00:04.699395502
# Collisions: 00:08:00.456872175
# Total: 00:08:40.903755126

if run:
    print(f'Setting up "{sim_name}"')
    create_data_dir(data_path)
    create_particles(sim_params, [electrons, ions], data_path)
    update_header(sim_params, project_path=project_path, data_path=data_path)

    compile_project(build_path, output=True)
    run_project(build_path + '/game_engine', output=True)

# # ===========================
# # ===== Post Processing =====
# # ===========================
xs = np.linspace(xmin, xmax, shape[0])
zs = np.linspace(zmin, zmax, shape[2])

block = True
# plot_density('electrons', 7500, data_path, xs, zs, block=block)
plot_temperature('electrons', 7500, data_path, xs, zs, block=block)
# plot_density('ions', 7500, data_path, xs, zs, block=block)
# plot_temperature('ions', 7500, data_path, xs, zs, block=block)
# plot_field_energy(data_path, block=block)
# plot_particle_energy(data_path)
