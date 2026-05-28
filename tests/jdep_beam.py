# #!/usr/bin/env python3
#

# from scripts.simulation import Simulation

import matplotlib.pyplot as plt
import numpy as np
from adios2 import FileReader, Stream
from scipy import constants
#
from scripts.pyforce import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'jdep_beam'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'


shape = (51, 51, 51)

xmin, xmax = -0.001, 0.001 # meters
ymin, ymax = -0.001, 0.001
zmin, zmax = -0.001, 0.001

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 1.0e-13 # seconds
t_end = 100 * dt #1.0e-10 # seconds
nt = int(t_end / dt) + 1

save_interval = 1

# =====================
# ===== Particles =====
# =====================
px_range = (0.0, dx) # meters
py_range = (0.0, dy)
pz_range = (0.0, dz)

single_particle = Particles(
    name='sp',
    mass=constants.m_e,
    charge=-1.0,
    atomic_number=0,
    tracer=True,
    temp=(1.0e7, 0.0, 0.0), # velocity
    density=1.0, # m^-3,
    ppc=(1, 1, 1),
    distribution='sp_jdep',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# ==========================================
# ===== Collisions and Particle Params =====
# ==========================================
particle_params = ParticleParams(
    save_interval=save_interval,
    bc_depth=2,
    push_type=ParticlePushType.Ballistic,
    particle_bcs=ParticleBCType.Periodic,
    interp_order=2,
    particle_data=(single_particle,)
)

# ==================================
# ===== Electromagnetic Params =====
# ==================================
em_params = EMParams(
    save_interval=save_interval
)

# ==========================
# ===== Metrics Params =====
# ==========================
metric_params = Metrics(
    data_path,
    (
        # MetricType.ParticleEnergy,
        # MetricType.FieldEnergy,
        MetricType.FieldDump,
        MetricType.ParticleDump,
        # MetricType.ParticleDiagnostics,
    )
)

# ============================
# ===== Simulation Class =====
# ============================
sim_params = Simulation(
    name=sim_name,
    shape=shape,
    nthreads=1,
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
    collisions_enabled=False,
    # push_enabled=False,
    jdep_enabled=True,
    # em_enabled=False
)

# # ===========================
# # ===== Compile and Run =====
# # ===========================
# print(f'Setting up "{sim_name}"')
# create_data_dir(data_path)
# create_particles(sim_params, single_particle, data_path)
# update_header(sim_params, project_path=project_path)
#
# compile_project(build_path, output=True)
# run_project(build_path + '/game_engine', output=True)

# ===========================
# ===== Post Processing =====
# ===========================
# with FileReader(data_path + f'/sp_tracer_{nt:010d}.bp') as f:
#     position = f.read("Position").reshape(-1, 3)

xs = np.linspace(xmin, xmax, shape[0], endpoint=True)
zs = np.linspace(zmin, zmax, shape[2], endpoint=True)

# Jxs = []
# Jys = []
# Jzs = []
# time = []
for n in range(0, nt, 1):
    with FileReader(data_path + f'/fields_{n:010d}.bp') as f:
        Jxs = f.read("Jx")[20:-20, 20:-20, shape[2] // 2].T

    fig, ax = plt.subplots(1, 1, figsize=(6, 12), layout='constrained')

    a0 = ax.pcolormesh(xs[20:-21], zs[20:-20], Jxs, cmap='jet', shading='nearest')
    ax.set_xticks(xs[20:-21])
    ax.set_yticks(zs[20:-20])

    ax.grid(which='both')
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_title('Jx')


    plt.colorbar(a0, ax=ax)

    plt.savefig(data_path + f'/jdep_{n}.png')
    plt.close(fig)
    # plt.show()