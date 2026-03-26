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
sim_name = 'uniform_E_field'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'

'''
Tests from https://iopscience.iop.org/article/10.3847/1538-4365/aab114
'''

shape = (16, 16, 16)

xmin, xmax = -1.0, 4e17 # meters
ymin, ymax = -1.0, 4e17
zmin, zmax = -1.0, 4e17

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 1e3 # seconds
t_end = 1e9 # seconds
nt = int(t_end / dt) + 1
cfl = constants.c * dt * np.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

save_interval = 200

# =====================
# ===== Particles =====
# =====================
mass = 1.0
charge = 1.0 / constants.e
Ex_amp = constants.c

gamma_half = np.sqrt(1.0 + (charge * constants.e * Ex_amp * 0.5 * dt / (mass * constants.c))**2)
x_half = (mass * constants.c**2) / (charge * constants.e * Ex_amp) * (gamma_half - 1.0)

px_range = (x_half, dx) # meters
py_range = (ymax / 2, dy)
pz_range = (zmax / 2, dz)

single_particle = Particles(
    name='sp',
    mass=mass,
    charge=charge,
    atomic_number=0,
    tracer=True,
    temp=(0.0, 0.0, 0.0), # eV
    density=1.0, # m^-3,
    ppc=(1, 1, 1),
    distribution='sp_uniformE',
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
    particle_data=(single_particle,)
)

# ============================
# ===== Simulation Class =====
# ============================
pushers = [
    ParticlePushType.Boris,
    ParticlePushType.HC
]
sim_names = [sim_name + '_' + pusher.split(':')[-1] for pusher in pushers]
Ex_applied = np.full((shape[0] - 1, shape[1], shape[2]), Ex_amp)

# for pusher, name in zip(pushers, sim_names):
#     data_path = project_path + f'/data/{name}'
#     fields_path = data_path + f'/{name}_applied_fields.bp'
#
#     with Stream(fields_path, 'w') as f:
#         f.write('Ex', Ex_applied, Ex_applied.shape, (0, 0, 0), Ex_applied.shape)
#
#     em_params = EMParams(save_interval=save_interval, applied_fields=fields_path)
#     metric_params = Metrics(data_path, (MetricType.ParticleDump,))
#     particle_params.push_type = pusher
#
#     sim_params = Simulation(
#         name=name,
#         shape=shape,
#         nthreads=1,
#         dt=dt,
#         t_end=t_end,
#         nt=nt,
#         x_range=(xmin, xmax),
#         y_range=(ymin, ymax),
#         z_range=(zmin, zmax),
#         deltas=(dx, dy, dz),
#         em_params=em_params,
#         particle_params=particle_params,
#         metric_params=metric_params,
#         em_enabled=False,
#         jdep_enabled=False,
#         collisions_enabled=False
#     )
#
#     # ===========================
#     # ===== Compile and Run =====
#     # ===========================
#     print(f'Setting up "{name}"')
#
#     if not os.path.exists(data_path):
#         print(f'Creating simulation data directory "{data_path}"...')
#         os.makedirs(data_path)
#
#     create_particles(sim_params, single_particle, data_path)
#     update_header(sim_params, project_path=project_path)
#
#     subprocess.run(
#         ['meson', 'compile', '-C', build_path, '-j4'],
#         check=True,
#         # stdout=subprocess.PIPE,
#         # stderr=subprocess.STDOUT
#     )
#
#     subprocess.run(
#         build_path + '/game_engine',
#         check=True,
#         # stdout=subprocess.PIPE,
#         # stderr=subprocess.STDOUT
#     )

# ===========================
# ===== Post Processing =====
# ===========================
boris_path = project_path + f'/data/{sim_names[0]}'
boris_gammas = []
boris_pos = []
time = []
for n in range(save_interval, nt, save_interval):
    file = f'/{single_particle.name}_tracer_{n:010d}.bp'
    with FileReader(boris_path + file) as f:
        boris_gammas.append(f.read('Gamma'))
        boris_pos.append(f.read("Position"))
        time.append(f.read("Time"))

boris_gammas = np.array(boris_gammas).flatten()
boris_pos = np.array(boris_pos).reshape(-1, 3)
time = np.array(time).flatten()

hc_path = project_path + f'/data/{sim_names[1]}'
hc_gammas = []
hc_pos = []
for n in range(save_interval, nt, save_interval):
    file = f'/{single_particle.name}_tracer_{n:010d}.bp'
    with FileReader(hc_path + file) as f:
        hc_gammas.append(f.read('Gamma'))
        hc_pos.append(f.read("Position"))

hc_gammas = np.array(hc_gammas).flatten()
hc_pos = np.array(hc_pos).reshape(-1, 3)

gamma_an = np.sqrt(1.0 + (charge * constants.e * Ex_amp * time / (mass * constants.c))**2)
x_an = (mass * constants.c**2) / (charge * constants.e * Ex_amp) * (gamma_an - 1.0)
v_an = (charge * constants.e * Ex_amp * time) / (mass * gamma_an)

fig, ax = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')

boris_pos = np.abs(boris_pos[:, 0] - x_an) / np.abs(x_an)
hc_pos = np.abs(hc_pos[:, 0] - x_an) / np.abs(x_an)

mark_every = time.shape[0] // 20

ax[0].set_xlabel('time')
ax[0].set_ylabel(r'|x - x_{an}| / |x_{an}|')
ax[0].set_yscale('log')
ax[0].set_xlim([0, 10e8])
ax[0].set_ylim([1e-15, 1])
ax[0].plot(time, boris_pos, ls='--', c='r', marker='P', ms=8, markevery=mark_every, label='Boris')
ax[0].plot(time, hc_pos, ls='--', c='b', marker='D', ms=8, markevery=mark_every, fillstyle='none', label='HC')
ax[0].legend()

boris_gammas = np.abs(boris_gammas - gamma_an) / gamma_an
hc_gammas = np.abs(hc_gammas - gamma_an) / gamma_an

ax[1].set_xlabel('time')
ax[1].set_ylabel(r'$|\gamma - \gamma_{an}|$ / $\gamma_{an}$')
# ax[1].set_yscale('log')
ax[1].set_xlim([0, 10e8])
# ax[1].set_ylim([1e-16, 1e-10])
ax[1].plot(time, boris_gammas, ls='--', c='r', marker='P', ms=8, markevery=mark_every, label='Boris')
# ax[1].plot(time, hc_gammas, ls='--', c='b', marker='D', ms=8, markevery=mark_every, fillstyle='none', label='HC')
ax[1].legend()

plt.show()
