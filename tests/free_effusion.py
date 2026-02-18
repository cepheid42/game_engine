#!/usr/bin/env python3

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
sim_name = 'free_effusion'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'

'''
Test 1 from 
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
# i_temp = 0.03878 # eV, 3/2 kT -> T = 300K
i_den = 5.0e11 # m^-3
ppc = (10, 10, 1)
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
    interp_order=1,
    particle_data=(electrons, ions)
)

# ==========================
# ===== Metrics Params =====
# ==========================
metric_params = Metrics(
    data_path,
    (MetricType.ParticleDump, MetricType.FieldDump)
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
    em_enabled=False,
    jdep_enabled=False,
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
e_weights = []
e_energies = []
for n in range(0, nt, save_interval):
    file = f'/electrons_dump_{n:010d}.bp'
    with FileReader(data_path + file) as f:
        weight = f.read('Weight')
        gammas = f.read('Gamma')
        energy = weight * (gammas - 1.0) * constants.m_e * constants.c**2
        # velocity = f.read('Velocity')
        # v = weight * velocity**2
        # energy = 0.5 * constants.m_e * v
        e_weights.append(weight.sum())
        e_energies.append(energy.sum())

e_weights = np.array(e_weights) / e_den
e_energies = np.array(e_energies) * 1e7

i_weights = []
i_energies = []
for n in range(0, nt, save_interval):
    file = f'/ions_dump_{n:010d}.bp'
    with FileReader(data_path + file) as f:
        weight = f.read('Weight')
        gammas = f.read('Gamma')
        energy = weight * (gammas - 1.0) * constants.m_p * constants.c**2
        # velocity = f.read('Velocity')
        # v = weight * velocity**2
        # energy = 0.5 * constants.m_p * v
        i_weights.append(weight.sum())
        i_energies.append(energy.sum())

i_weights = np.array(i_weights) / i_den
i_energies = np.array(i_energies) * 1e7

T_0e = 11600.0 # Kelvin ~= 1 eV
T_0i = 300.0

wall_area = 4.0 # 1 m^2 for each wall
e_theory_slope = 0.25 * e_den * wall_area * np.sqrt(8.0 * constants.k * T_0e / (np.pi * constants.m_e))
i_theory_slope = 0.25 * i_den * wall_area * np.sqrt(8.0 * constants.k * T_0i / (np.pi * constants.m_p))

time = np.linspace(0, t_end, e_weights.shape[0]) * 1.0e6

e_theory = 1.0 - (e_theory_slope * time) / (e_den * 1e6)
i_theory = 1.0 - (i_theory_slope * time) / (i_den * 1e6)

e_energy_theory = np.genfromtxt('./data/effusion_electron_energy.csv', delimiter=',')
i_energy_theory = np.genfromtxt('./data/effusion_ion_energy.csv', delimiter=',')

fig, ax = plt.subplots(2, 2, figsize=(8, 8), layout='constrained', sharex=True)

ax[0, 0].plot(time, i_theory, 'k--', label='Theory')
ax[0, 0].plot(time, i_weights, 'r--', label='TriForce')
# ax[0, 1].set_ylim([0, 1])
ax[0, 0].legend()
ax[0, 0].set_ylabel('Norm. density ions [-]')

ax[0, 1].plot(time, e_theory, 'k--', label='Theory')
ax[0, 1].plot(time, e_weights, 'r--', label='TriForce')
ax[0, 1].set_ylim([0, 1])
ax[0, 1].legend()
ax[0, 1].set_ylabel('Norm. density electrons [-]')

ax[1, 0].plot(i_energy_theory[:, 0], i_energy_theory[:, 1], 'k--', label='Parodi')
ax[1, 0].plot(time, i_energies, 'r--', label='TriForce')
ax[1, 0].legend()
ax[1, 0].set_ylabel(r'Energy ions [$10^{-7}$ J]')
ax[1, 0].set_xlabel(r'Time [$\mu$s]')
ax[1, 0].set_xlim([0, time[-1]])

ax[1, 1].plot(e_energy_theory[:, 0], e_energy_theory[:, 1], 'k--', label='Parodi')
ax[1, 1].plot(time, e_energies, 'r--', label='TriForce')
ax[1, 1].set_ylim([0, 1.2])
ax[1, 1].legend()
ax[1, 1].set_ylabel(r'Energy electrons [$10^{-7}$ J]')
ax[1, 1].set_xlabel(r'Time [$\mu$s]')
ax[1, 1].set_xlim([0, time[-1]])

plt.show()
