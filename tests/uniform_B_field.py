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
data_path = project_path + f'/data/{sim_name}'

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

save_interval = 20

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
    particle_bcs='outflow',
    bc_depth=0,
    interp_order=1,
    particle_data=(singleton,)
)

# ==================================
# ===== Electromagnetic Params =====
# ==================================
Bz_applied = np.full((shape[0] - 1, shape[1] - 1, shape[2]), -Bz_amp)
with Stream(data_path + f'/{sim_name}_applied_fields.bp', 'w') as f:
    f.write('Bz', Bz_applied, Bz_applied.shape, (0, 0, 0), Bz_applied.shape)

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

# ===========================
# ===== Post Processing =====
# ===========================
pos = []
gammas = []
time = []
for n in range(0, nt, save_interval):
    file = f'/singleton_dump_{n:010d}.bp'
    with FileReader(data_path + file) as f:
        pos.append(f.read("Position"))
        gammas.append(f.read("Gamma"))
        time.append(f.read("Time")[0])

time = np.array(time)
pos = np.array(pos).reshape(-1, 3)
gammas = np.array(gammas).reshape(-1, 1)

Rc_an = gamma_an * v_perp / Bz_amp
p_Rc = plt.Circle((0, 0), Rc_an, fill=False)

theta_an = -omega_c * time
theta_c = -np.unwrap(np.arctan2(pos[:, 1], pos[:, 0])) + np.pi
gamma_data = np.abs(gammas - gamma_an) / gamma_an

fig, ax = plt.subplots(1, 3, figsize=(18, 6), layout='constrained')

ax[0].set_xlabel('x')
ax[0].set_ylabel('y')
ax[0].set_xlim([xmin, xmax])
ax[0].set_ylim([ymin, ymax])
ax[0].add_patch(p_Rc)
ax[0].plot(pos[:, 0], pos[:, 1])

ax[1].set_xlabel(r't / $T_c$')
ax[1].set_ylabel(r'$|\gamma - \gamma_{an}|$ / $\gamma_{an}$')
ax[1].set_xlim([0, 100])
ax[1].set_ylim([1e-16, 1e-13])
ax[1].set_yscale('log')
ax[1].plot(time / Tc, gamma_data)

ax[2].set_xlabel(r't / $T_c$')
ax[2].set_ylabel(r'$|\theta_c - \theta_{an}|$')
ax[2].set_xlim([0, 100])
ax[2].set_ylim([0, 0.25])
ax[2].plot(time / Tc, np.abs(theta_c - theta_an))
# ax[2].plot(time / Tc, theta_an)
# ax[2].plot(time / Tc, theta_c)

plt.show()