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
sim_name = 'uniform_B_field'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
# data_path = project_path + f'/data/{sim_name}'

'''
Tests from https://iopscience.iop.org/article/10.3847/1538-4365/aab114
'''

shape = (16, 16, 16)

mass = 1.0
charge = 1.0 / constants.e # charge * e == 1
gamma_an = 1e6
v_perp = constants.c * np.sqrt(1 - 1 / gamma_an**2)
Bz_amp = gamma_an * v_perp
omega_c = Bz_amp / gamma_an
Tc = 2 * np.pi / omega_c

xmin, xmax = -1.1, 1.1 # meters
ymin, ymax = -1.1, 1.1
zmin, zmax = -1.1, 1.1

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = Tc / 100 # seconds
t_end = 10000 * dt # seconds
nt = int(t_end / dt) + 1
cfl = constants.c * dt * np.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

save_interval = 5

# =====================
# ===== Particles =====
# =====================
px_range = (-1.0, 0.0) # meters
py_range = (0.0, 0.0)
pz_range = (0.0, 0.0)

singleton = Particles(
    name='singleton',
    mass=mass,
    charge=charge,
    atomic_number=0,
    # actually velocity for sp_uniformB distribution
    temp=(0.0, -v_perp, 0.0),
    density=1.0, # m^-3,
    ppc=(1, 1, 1),
    distribution='sp_uniformB',
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

    Bz_applied = np.full((shape[0] - 1, shape[1] - 1, shape[2]), -Bz_amp)
    with Stream(fields_path, 'w') as f:
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
        # stdout=subprocess.PIPE,
        # stderr=subprocess.STDOUT
        # stdout=subprocess.DEVNULL,
        # stderr=subprocess.DEVNULL
    )

#     subprocess.run(build_path + '/game_engine').check_returncode()
#
# # ===========================
# # ===== Post Processing =====
# # ===========================
# boris_path = project_path + f'/data/{sim_names[0]}'
# boris_gammas = []
# boris_pos = []
# for n in range(0, nt, save_interval):
#     file = f'/singleton_dump_{n:010d}.bp'
#     with FileReader(boris_path + file) as f:
#         boris_gammas.append(f.read('Gamma'))
#         boris_pos.append(f.read("Position"))
#
# boris_gammas = np.array(boris_gammas).flatten()
# boris_pos = np.array(boris_pos).reshape(-1, 3)
#
# hc_path = project_path + f'/data/{sim_names[1]}'
# hc_gammas = []
# hc_pos = []
# for n in range(0, nt, save_interval):
#     file = f'/singleton_dump_{n:010d}.bp'
#     with FileReader(hc_path + file) as f:
#         hc_gammas.append(f.read('Gamma'))
#         hc_pos.append(f.read("Position"))
#
# hc_gammas = np.array(hc_gammas).flatten()
# hc_pos = np.array(hc_pos).reshape(-1, 3)
#
# time = np.linspace(0, t_end, boris_gammas.shape[0])
#
# theta_an = -omega_c * time
# Rc_an = gamma_an * v_perp / Bz_amp
# p_Rc = plt.Circle((0, 0), Rc_an, fill=False)
#
# boris_theta_c = -np.unwrap(np.arctan2(boris_pos[:, 1], boris_pos[:, 0])) - np.pi
# hc_theta_c = -np.unwrap(np.arctan2(hc_pos[:, 1], hc_pos[:, 0])) - np.pi
#
# boris_thetas = np.abs(boris_theta_c - theta_an)
# hc_thetas = np.abs(hc_theta_c - theta_an)
#
# boris_gammas = np.abs(boris_gammas - gamma_an) / gamma_an
# hc_gammas = np.abs(hc_gammas - gamma_an) / gamma_an
#
# fig, ax = plt.subplots(1, 3, figsize=(18, 6), layout='constrained')
#
# ax[0].set_xlabel('x')
# ax[0].set_ylabel('y')
# ax[0].set_xlim([xmin, xmax])
# ax[0].set_ylim([ymin, ymax])
# ax[0].set_xlim([-5, 5])
# ax[0].set_ylim([-5, 5])
# ax[0].add_patch(p_Rc)
# ax[0].plot(boris_pos[:, 0], boris_pos[:, 1], c='r', label='Boris')
# ax[0].plot(hc_pos[:, 0], hc_pos[:, 1], c='b', label='HC')
# ax[0].legend()
#
# ax[1].set_xlabel(r't / $T_c$')
# ax[1].set_ylabel(r'$|\gamma - \gamma_{an}|$ / $\gamma_{an}$')
# ax[1].set_xlim([0, 100])
# # ax[1].set_ylim([1e-16, 1e-13])
# ax[1].set_yscale('log')
# ax[1].plot(time[::20] / Tc, boris_gammas[::20], c='r', marker='P', ms=8, markevery=40, label='Boris')
# ax[1].plot(time[::20] / Tc, hc_gammas[::20], c='b', marker='D', ms=8, markevery=40, fillstyle='none', label='HC')
# ax[1].legend()
#
# ax[2].set_xlabel(r't / $T_c$')
# ax[2].set_ylabel(r'$|\theta_c - \theta_{an}|$')
# ax[2].set_xlim([0, 100])
# ax[2].set_ylim([0, 0.25])
# ax[2].plot(time[::20] / Tc, boris_thetas[::20], c='r', marker='P', ms=8, markevery=40, label='Boris')
# ax[2].plot(time[::20] / Tc, hc_thetas[::20], c='b', marker='D', ms=8, markevery=40, fillstyle='none', label='HC')
# ax[2].legend()
#
# plt.show()