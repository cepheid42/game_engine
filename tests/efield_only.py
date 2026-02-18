#!/usr/bin/env python3

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from adios2 import FileReader
from scipy import constants

from scripts.particle_generation import create_particles
from scripts.domain_params import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'efield_only'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'

'''
Test 2 from 
https://pubs.aip.org/aip/pop/article/32/1/013902/3330730/Step-by-step-verification-of-particle-in-cell
'''

shape = (207, 207, 2)

dx = dy = 5.0e-3 # meters
dz = 1.0

xmin, xmax = -3*dx, 1.0 + 3*dx # meters
ymin, ymax = -3*dy, 1.0 + 3*dy
zmin, zmax = 0.0, 1.0

dt = 2.5e-9 # seconds
t_end = 5.0e-6 # seconds
nt = int(t_end / dt) + 1
cfl = constants.c * dt * np.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

save_interval = 50

# =====================
# ===== Particles =====
# =====================
e_temp = 1.0 # eV, ~11600K
e_den = 5.0e11 # m^-3
i_temp = 0.02585 # eV, kT -> T = 300K
i_den = 5.0e11 # m^-3
ppc = (5, 5, 1) # (10, 10, 1)
px_range = (0.0, 1.0) # meters
py_range = (0.0, 1.0)
pz_range = (zmin, zmax)

electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-1,
    atomic_number=0,
    temp=(e_temp, e_temp, e_temp), # eV
    density=e_den, # m^-3,
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
    temp=(i_temp, i_temp, i_temp), # eV
    density=i_den, # m^-3,
    ppc=ppc,
    distribution='thermal',
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
    interp_order=2,
    particle_data=(electrons, ions)
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
    nthreads=48,
    dt=dt,
    t_end=t_end,
    nt=nt,
    cfl=cfl,
    x_range=(xmin, xmax),
    y_range=(ymin, ymax),
    z_range=(zmin, zmax),
    deltas=(dx, dy, dz),
    particle_params=particle_params,
    metric_params=metric_params,
    # em_enabled=False,
    # jdep_enabled=False,
    collisions_enabled=False
)

# ===========================
# ===== Compile and Run =====
# ===========================
print(f'Setting up "{sim_name}"')

create_particles(sim_params, electrons, data_path)
create_particles(sim_params, ions, data_path)
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
eden_tf = []
# eKE_tf = []
for n in range(0, nt, save_interval):
    file = f'/electrons_dump_{n:010d}.bp'
    with FileReader(data_path + file) as f:
        weight = f.read('Weight')
        # gammas = f.read('Gamma')
        # energy = weight * (gammas - 1.0) * constants.m_e * constants.c**2
        eden_tf.append(weight.sum())
        # eKE_tf.append(energy.sum())

eden_tf = np.array(eden_tf) / e_den
# eKE_tf = np.array(eKE_tf) * 1e7

iden_tf = []
# iKE_tf = []
for n in range(0, nt, save_interval):
    file = f'/ions_dump_{n:010d}.bp'
    with FileReader(data_path + file) as f:
        weight = f.read('Weight')
        # gammas = f.read('Gamma')
        # energy = weight * (gammas - 1.0) * constants.m_p * constants.c**2
        iden_tf.append(weight.sum())
        # iKE_tf.append(energy.sum())

iden_tf = np.array(iden_tf) / i_den
# iKE_tf = np.array(iKE_tf) * 1e7

time = np.linspace(0, t_end, nt // save_interval + 1) * 1.0e6

# Load Parodi 2025 plot data
eden_pp = np.genfromtxt('./data/efield_only_eden.txt')
iden_pp = np.genfromtxt('./data/efield_only_iden.txt')

fig, ax = plt.subplots(2, 2, figsize=(8, 8), layout='constrained', sharex=True)

ax[0, 0].plot(time, iden_pp, ls='--', c='b', marker='v', ms=6, markevery=3, label='Parodi')
ax[0, 0].plot(time, iden_tf, ls='--', c='r', marker='X', ms=6, markevery=5, label='Ions (TF)')
ax[0, 0].legend()
ax[0, 0].set_ylabel('Norm. density [-]')
ax[0, 0].set_ylim([0.85, 1.0])

plt.show()
