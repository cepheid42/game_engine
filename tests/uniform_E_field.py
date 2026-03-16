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
# data_path = project_path + f'/data/{sim_name}'

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

    Ex_applied = np.full((shape[0] - 1, shape[1] - 1, shape[2]), Ex_amp)
    with Stream(fields_path, 'w') as f:
        f.write('Ex', Ex_applied, Ex_applied.shape, (0, 0, 0), Ex_applied.shape)

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
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    ).check_returncode()

    subprocess.run(build_path + '/game_engine').check_returncode()

    break

# ===========================
# ===== Post Processing =====
# ===========================
boris_path = project_path + f'/data/{sim_names[0]}'
boris_gammas = []
boris_pos = []
for n in range(0, nt, save_interval):
    file = f'/singleton_dump_{n:010d}.bp'
    with FileReader(boris_path + file) as f:
        boris_gammas.append(f.read('Gamma'))
        boris_pos.append(f.read("Position"))

boris_gammas = np.array(boris_gammas).flatten()
boris_pos = np.array(boris_pos).reshape(-1, 3)

hc_path = project_path + f'/data/{sim_names[1]}'
hc_gammas = []
hc_pos = []
for n in range(0, nt, save_interval):
    file = f'/singleton_dump_{n:010d}.bp'
    with FileReader(hc_path + file) as f:
        hc_gammas.append(f.read('Gamma'))
        hc_pos.append(f.read("Position"))

hc_gammas = np.array(hc_gammas).flatten()
hc_pos = np.array(hc_pos).reshape(-1, 3)

time = np.linspace(0, t_end, boris_gammas.shape[0])

gamma_an = np.sqrt(1.0 + (charge * constants.e * Ex_amp * time / (mass * constants.c))**2)
x_an = (mass * constants.c**2) / (charge * constants.e * Ex_amp) * (gamma_an - 1.0)
v_an = (charge * constants.e * Ex_amp * time) / (mass * gamma_an)

fig, ax = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')

# for b, a in zip(boris_gammas, gamma_an):
#     print(f'{b:25.25f} | {a:25.25f}')
boris_pos = np.abs(boris_pos[:, 0] - x_an) / np.abs(x_an)
hc_pos = np.abs(hc_pos[:, 0] - x_an) / np.abs(x_an)

ax[0].set_xlabel('time (ns)')
ax[1].set_ylabel(r'|x - x_{an}| / |x_{an}|')
ax[0].set_yscale('log')
ax[0].set_xlim([0, 10e8])
# ax[0].set_ylim([1e-15, 1])
ax[0].plot(time, boris_pos, ls='--', c='r', marker='P', ms=8, markevery=100, label='Boris')
ax[0].plot(time, hc_pos, ls='--', c='b', marker='D', ms=8, markevery=100, fillstyle='none', label='HC')
ax[0].legend()

boris_gammas = np.abs(boris_gammas - gamma_an) / gamma_an
hc_gammas = np.abs(hc_gammas - gamma_an) / gamma_an

ax[1].set_xlabel('time (ns)')
ax[1].set_ylabel(r'$|\gamma - \gamma_{an}|$ / $\gamma_{an}$')
ax[1].set_yscale('log')
ax[1].set_xlim([0, 10e8])
# ax[1].set_ylim([1e-16, 1e-10])
ax[1].plot(time, boris_gammas, ls='--', c='r', marker='P', ms=8, markevery=100, label='Boris')
ax[1].plot(time, hc_gammas, ls='--', c='b', marker='D', ms=8, markevery=100, fillstyle='none', label='HC')
ax[1].legend()

plt.show()
