#!/usr/bin/env python3

import numpy as np
import subprocess
import matplotlib.pyplot as plt
from scipy import constants
from adios2 import FileReader

from scripts.domain_params import *
from scripts.particle_generation import create_particles

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'ionization'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'

shape = (2, 2, 2)

xmin, xmax = 0.0, 1.0e-6
ymin, ymax = 0.0, 1.0e-6
zmin, zmax = 0.0, 1.0e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 5.0e-18
t_end = 3.18e-15
nt = int(t_end / dt) + 1
cfl = constants.c * dt * np.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

save_interval = 10

# =====================
# ===== Particles =====
# =====================
px_range=(xmin, xmax)
py_range=(ymin, ymax)
pz_range=(zmin, zmax)

electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-1,
    atomic_number=0,
    temp=(500, 0, 0), # eV
    density=1.1e27, # m^-3,
    ppc=(100, 10, 10),
    distribution='constant',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

electron_products = Particles(
    name='electron_products',
    mass=constants.m_e,
    charge=-1,
    atomic_number=0,
    temp=(0, 0, 0), # eV
    density=0.0, # m^-3,
    ppc=(0, 0, 0),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

neutral_aluminum = Particles(
    name='Al',
    mass=26.9815384 * constants.atomic_mass + 13.0 * constants.m_e,
    charge=0,
    atomic_number=13,
    temp=(0, 0, 0), # eV
    density=6.6e28, # m^-3,
    ppc=(20, 20, 20),
    distribution='constant',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

ionized_aluminum = Particles(
    name='Al_products',
    mass=26.9815384 * constants.atomic_mass + 12.0 * constants.m_e,
    charge=1,
    atomic_number=13,
    temp=(0, 0, 0), # eV
    density=0, # m^-3,
    ppc=(0, 0, 0),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range,
)

# ==========================================
# ===== Collisions and Particle Params =====
# ==========================================
ionization_params = IonizationParams(
    products=('electron_products', 'Al_products'),
    ionization_energy=5.9858,
    cross_section_file='/tests/cross_section_data/eAl0_ionization_eV_m2.txt'
)

particle_params = ParticleParams(
    save_interval=save_interval,
    particle_bcs='periodic',
    particle_data=(electrons, electron_products, neutral_aluminum, ionized_aluminum),
    collisions=(
        Collision(
            groups=('electrons', 'Al'),
            channels=('ionization',),
            ionization=ionization_params
        ),
    )
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
    nthreads=4,
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
    push_enabled=False,
    jdep_enabled=False
)

# ===========================
# ===== Compile and Run =====
# ===========================
print(f'Setting up "{sim_name}"')
create_particles(sim_params, electrons, data_path)
create_particles(sim_params, neutral_aluminum, data_path)
update_header(sim_params, project_path=project_path, ionization_test_override=True)

subprocess.run(
    ['meson', 'compile', '-C', build_path, '-j4'],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL
).check_returncode()

subprocess.run(build_path + '/game_engine').check_returncode()

# ===========================
# ===== Post Processing =====
# ===========================
step = save_interval
start = step
stop = nt

ionized = []
for n in range(start, stop, step):
    file_name = f'/Al_products_dump_{n:010d}.bp'
    with FileReader(data_path + file_name) as f:
        ionized.append(f.read('Weight').sum())
ionized = np.array(ionized)

v_beam = 1.32523e7
sigma = 1.428e-20
e_den = 1.1e27
omega_p_inv = 0.53e-15

time_thry = np.linspace(0, 6, 100) * omega_p_inv
charge_thry = 1.0 - np.exp(-v_beam * e_den * sigma * time_thry)

time = dt * np.arange(start, stop, step)
mean_ion_charge = ionized / 6.6e10

plt.style.use('/home/cepheid/TriForce/game_engine/data/triforce.mplstyle')
good_colors = ['#db6d00', '#006ddb', '#920000', '#52a736', '#9B30FF']

fig, ax = plt.subplots(figsize=(8, 8), layout='constrained')
ax.set_title('Mean ion charge', fontsize=16, loc='left')
# ax[1].set_title('(b) Particle count', fontsize=16, loc='left')
ax.set_ylim([0.0, 0.5])
ax.set_xlabel('Time (s)')
ax.tick_params(axis='x', pad=8)
ax.grid()

ax.plot(time_thry, charge_thry, '-k', label="Theory")
ax.plot(time, mean_ion_charge, c=good_colors[0], label=r'$M_{\mathrm{p}}=10^0$',
        linestyle='none', marker='o', mec=good_colors[0], mew=1.5, ms=5)

ax.legend(loc=4)
# plt.savefig(data_path + f'/ionization_test_mean_ion_charge.png')
# plt.close(fig)
plt.show()
