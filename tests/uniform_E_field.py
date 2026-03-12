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
data_path = project_path + f'/data/{sim_name}'

'''
Tests from https://iopscience.iop.org/article/10.3847/1538-4365/aab114
'''

shape = (16, 16, 16)

xmin, xmax = -1.0, 4e17 # meters
ymin, ymax = 0.0, 4e17
zmin, zmax = 0.0, 4e17

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 1e3 # seconds
t_end = 1e9 # seconds
nt = int(t_end / dt) + 1
cfl = constants.c * dt * np.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

save_interval = nt // 5000

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

singleton = Particles(
    name='singleton',
    mass=mass,
    charge=charge,
    atomic_number=0,
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
    particle_bcs='outflow',
    bc_depth=0,
    interp_order=1,
    particle_data=(singleton,)
)

# ==================================
# ===== Electromagnetic Params =====
# ==================================

Ex_applied = np.full((shape[0] - 1, shape[1] - 1, shape[2]), Ex_amp)
with Stream(data_path + f'/{sim_name}_applied_fields.bp', 'w') as f:
    f.write('Ex', Ex_applied, Ex_applied.shape, (0, 0, 0), Ex_applied.shape)

em_params = EMParams(
    save_interval=save_interval,
    applied_fields=data_path + f'/{sim_name}_applied_fields.bp'
)

# ==========================
# ===== Metrics Params =====
# ==========================
metric_params = Metrics(
    data_path,
    (MetricType.ParticleDump,)
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
print(f'Setting up "{sim_name}"')

# if not os.path.exists(data_path):
#     print(f'Creating simulation data directory "{data_path}"...')
#     os.makedirs(data_path)
#
# create_particles(sim_params, singleton, data_path)
# update_header(sim_params, project_path=project_path)
#
# subprocess.run(
#     ['meson', 'compile', '-C', build_path, '-j4'],
#     stdout=subprocess.DEVNULL,
#     stderr=subprocess.DEVNULL
# ).check_returncode()
#
# subprocess.run(build_path + '/game_engine').check_returncode()

# ===========================
# ===== Post Processing =====
# ===========================
gammas = []
vel = []
pos = []
for n in range(0, nt, save_interval):
    file = f'/singleton_dump_{n:010d}.bp'
    with FileReader(data_path + file) as f:
        gammas.append(f.read('Gamma'))
        vel.append(f.read("Velocity"))
        pos.append(f.read("Position"))

gammas = np.array(gammas).flatten()
pos = np.array(pos).reshape(-1, 3)
vel = np.array(vel).reshape(-1, 3)

time = np.linspace(0, t_end, nt // save_interval + 1)

gamma_an = np.sqrt(1.0 + (charge * constants.e * Ex_amp * time / (mass * constants.c))**2)
x_an = (mass * constants.c**2) / (charge * constants.e * Ex_amp) * (gamma_an - 1.0)
v_an = (charge * constants.e * Ex_amp * time) / (mass * gamma_an)

# x_error = np.mean((pos[:, 0] - x_an)**2)
# v_error = np.mean((vel[:, 0] - v_an)**2)
# print(f'MSE position = {x_error}')
# print(f'MSE velocity = {v_error}')

fig, ax = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')

pos_data = np.abs(pos[:, 0] - x_an) / np.abs(x_an)
ax[0].set_xlabel('time (ns)')
ax[0].set_ylabel('x (m)')
ax[0].set_yscale('log')
ax[0].set_xlim([0, 10e8])
# ax[0].set_ylim([1e-15, 1])
# ax[0].plot(time, x_an, label='Analytic')
ax[0].plot(time, pos_data, label='TF')
ax[0].legend()

gamma_data = np.abs(gammas - gamma_an) / gamma_an
ax[1].set_xlabel('time (ns)')
ax[1].set_ylabel(r'$\gamma$')
ax[1].set_yscale('log')
ax[1].set_xlim([0, 10e8])
# ax[1].set_ylim([1e-16, 1e-10])
# ax[1].plot(time, gamma_an, label='Analytic')
# ax[1].plot(time, gammas, label='TF')
ax[1].plot(time, gamma_data, label='TF')
ax[1].legend()

# ax[2].set_xlabel('time (ns)')
# ax[2].set_ylabel(r'$v_x$ (m/s)')
# ax[2].plot(time, v_an, label='Analytic')
# ax[2].plot(time, vel[:, 0], label='TF')
# ax[2].legend()

plt.show()
