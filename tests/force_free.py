#!/usr/bin/env python3

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from adios2 import FileReader, Stream
from scipy import constants

from scripts.particle_generation import create_particles
from scripts.domain_params import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'force_free'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'

'''
Tests from https://iopscience.iop.org/article/10.3847/1538-4365/aab114
'''

shape = (16, 16, 16)

xmin, xmax = -1.0, 3.1e17 # meters
ymin, ymax = -1.0, 3.1e17
zmin, zmax = -1.0, 3.1e17

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 0.01 # seconds
t_end = 1e5 # seconds
nt = int(t_end / dt) + 1
cfl = constants.c * dt * np.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

save_interval = 1000

# =====================
# ===== Particles =====
# =====================
mass = 1.0
charge = 1.0 / constants.e

gamma_init = 1.0e6
vy = constants.c * np.sqrt(1 - 1 / gamma_init**2)
Bz_amp = 1.0
Ex_amp = -vy * Bz_amp

px_range = (0.0, dx) # meters
py_range = (0.0, dy)
pz_range = (0.0, dz)

singleton = Particles(
    name='singleton',
    mass=mass,
    charge=charge,
    atomic_number=0,
    temp=(0.0, vy, 0.0), # eV
    density=1.0, # m^-3,
    ppc=(1, 1, 1),
    distribution='sp_forcefree',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# ==========================================
# ===== Collisions and Particle Params =====
# ==========================================
particle_params = ParticleParams(
    save_interval=save_interval,
    bc_depth=0,
    interp_order=1,
    particle_data=(singleton,)
)

# ============================
# ===== Simulation Class =====
# ============================
pushers = [ParticlePushType.Boris, ParticlePushType.HC]
sim_names = [sim_name + '_' + pusher.split(':')[-1] for pusher in pushers]

for pusher, name in zip(pushers, sim_names):
    data_path = project_path + f'/data/{name}'
    fields_path = data_path + f'/{name}_applied_fields.bp'

    Ex_applied = np.full((shape[0] - 1, shape[1], shape[2]), Ex_amp)
    Bz_applied = np.full((shape[0] - 1, shape[1] - 1, shape[2]), Bz_amp)
    with Stream(fields_path, 'w') as f:
        f.write('Ex', Ex_applied, Ex_applied.shape, (0, 0, 0), Ex_applied.shape)
        f.write('Bz', Bz_applied, Bz_applied.shape, (0, 0, 0), Bz_applied.shape)

    em_params = EMParams(
        save_interval=save_interval,
        applied_fields=fields_path
    )

    metric_params = Metrics(
        data_path,
        (MetricType.ParticleDump,)
    )

    particle_params.push_type = pusher
    sim_params = Simulation(
        name=name,
        shape=shape,
        nthreads=1,
        dt=dt,
        t_end=t_end,
        nt=nt,
        cfl=cfl,
        x_range=(xmin, xmax),
        y_range=(ymin, ymax),
        z_range=(zmin, zmax),
        deltas=(dx, dy, dz),
        em_params=em_params,
        particle_params=particle_params,
        metric_params=metric_params,
        em_enabled=False,
        jdep_enabled=False,
        collisions_enabled=False
    )

    # ===========================
    # ===== Compile and Run =====
    # ===========================
    print(f'Setting up "{name}"')

    if not os.path.exists(data_path):
        print(f'Creating simulation data directory "{data_path}"...')
        os.makedirs(data_path)

    create_particles(sim_params, singleton, data_path)
    update_header(sim_params, project_path=project_path)

    subprocess.run(
        ['meson', 'compile', '-C', build_path, '-j4'],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

    subprocess.run(build_path + '/game_engine').check_returncode()

# ===========================
# ===== Post Processing =====
# ===========================
boris_path = project_path + f'/data/{sim_names[0]}'
boris_vel = []
boris_pos = []
time = []
for n in range(0, nt, save_interval):
    file = f'/singleton_dump_{n:010d}.bp'
    with FileReader(boris_path + file) as f:
        boris_vel.append(f.read('Velocity'))
        boris_pos.append(f.read("Position"))
        time.append(f.read("Time"))

boris_vel = np.array(boris_vel).reshape(-1, 3)
boris_pos = np.array(boris_pos).reshape(-1, 3)
time = np.array(time).flatten()

hc_path = project_path + f'/data/{sim_names[1]}'
hc_vel = []
hc_pos = []
for n in range(0, nt, save_interval):
    file = f'/singleton_dump_{n:010d}.bp'
    with FileReader(hc_path + file) as f:
        hc_vel.append(f.read('Velocity'))
        hc_pos.append(f.read("Position"))

hc_vel = np.array(hc_vel).reshape(-1, 3)
hc_pos = np.array(hc_pos).reshape(-1, 3)

time = np.linspace(0, t_end, boris_vel.shape[0])

fig, ax = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')

mark_every = 500
ax[0].set_xlabel('t')
ax[0].set_ylabel('x')
ax[0].set_xlim([0, 1e5])
# ax[0].set_ylim([-1.4e-3, 0])
ax[0].plot(time, -boris_pos[:, 0], c='r', marker='P', ms=8, markevery=mark_every, label='Boris')
ax[0].plot(time, -hc_pos[:, 0], c='b', marker='D', ms=8, markevery=mark_every, fillstyle='none', label='HC')
ax[0].legend()

ax[1].set_xlabel('t')
ax[1].set_ylabel(r'$v_x$')
ax[1].set_xlim([0, 1e5])
# ax[1].set_ylim([-2.5e-8, 0])
ax[1].plot(time, -boris_vel[:, 0], ls='--', c='r', marker='P', ms=8, markevery=mark_every, label='Boris')
ax[1].plot(time, -hc_vel[:, 0], ls='--', c='b', marker='D', ms=8, markevery=mark_every, fillstyle='none', label='HC')
ax[1].legend()

plt.show()
