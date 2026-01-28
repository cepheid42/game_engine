#!/usr/bin/env python3

import numpy as np
import subprocess
import matplotlib.pyplot as plt
from scipy import constants
from adios2 import FileReader

from scripts.domain_params import *
from scripts.particle_generation import create_particles

sim_name = 'ionization'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
particle_data = project_path + '/data/'

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

# ===== Particles =====
px_range=(xmin, xmax)
py_range=(ymin, ymax)
pz_range=(zmin, zmax)
save_interval = 10

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

# ===== Collisions and Particle Params =====
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
    em_enabled=False,
    push_enabled=False,
    jdep_enabled=False
)

print(f'Setting up "{sim_name}"')
create_particles(sim_params, electrons, particle_data)
create_particles(sim_params, neutral_aluminum, particle_data)
update_header(sim_params, project_path=project_path, ionization_test_override=True)

subprocess.run(
    ['meson', 'compile', '-C', build_path, '-j4'],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL
).check_returncode()

subprocess.run(build_path + '/game_engine').check_returncode()

step = save_interval
start = step
stop = nt

ionized = []
for n in range(start, stop, step):
    file_name = f'{sim_name}/Al_products_dump_{n:010d}.bp'
    with FileReader(particle_data + file_name) as f:
        ionized.append(f.read('Weight').sum())
ionized = np.asarray(ionized)

v_beam = 1.32523e7
sigma = 1.428e-20
e_den = 1.1e27
time_thry = np.linspace(0, 6, 100) * 0.53e-15
charge_thry = 1.0 - np.exp(-v_beam * e_den * sigma * time_thry)

time = dt * np.arange(start, stop, step)
mean_ion_charge = ionized / 6.6e10

fig, ax = plt.subplots(figsize=(8, 8), layout='constrained')
ax.plot(time_thry, charge_thry)
ax.plot(time, mean_ion_charge)

plt.show()
