#!/usr/bin/env python3

# import numpy as np
from scipy import constants
import math

from scripts.particle_generation import create_particles
from scripts.domain_params import *

project_path = '/home/cepheid/TriForce/game_engine'
particle_data = project_path + '/data'

shape = (751, 2, 751)

xmin, xmax = -15.0e-6, 15.0e-6
zmin, zmax = -15.0e-6, 15.0e-6

dx = (xmax - xmin) / (shape[0] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

ymin, ymax = 0.0, 0.01
dy = 0.01

dt = 8.0e-17
t_end = 4000 * dt #3.0e-13
nt = int(t_end / dt) + 1
cfl = constants.c * dt * math.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

# ===== Particles =====
ppc = (5, 1, 5)
px_range = (-5e-7, 5e-7)
py_range = (ymin, ymax)
pz_range = (-1e-5, 1e-5)
temp = 10000 / math.sqrt(3)

electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-1,
    atomic_number=0,
    temp=(temp, temp, temp), # eV
    density=8.5e27, # m^-3,
    ppc=ppc,
    distribution='relativistic',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

ions = Particles(
    name='ions',
    mass=constants.m_p,
    charge=1,
    atomic_number=1,
    temp=(temp, temp, temp), # eV
    density=8.5e27, # m^-3,
    ppc=ppc,
    distribution='relativistic',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

em_params = EMParams(
    save_interval=20,
    pml_depth=15,
    em_bcs=(1, 1, 2, 2, 1, 1),
)

coulomb_params = CoulombParams(coulomb_log=0.0)

ionization_params = IonizationParams(
    products=('electron', 'ions'),
    ionization_energy=13.6, # eV
    constant_cross_section=0.0
)

particle_params = ParticleParams(
    save_interval=20,
    particle_bcs='outflow',
    interp_order=2,
    particle_data=(electrons, ions),
    collisions=(
        Collision(
            groups=('electrons', 'ions'),
            channels=('coulomb', 'ionization'),
            coulomb=coulomb_params,
            ionization=ionization_params
        ),
    )
)

sim_params = Simulation(
    name='lsi_test',
    shape=shape,
    nthreads=32,
    dt=dt,
    t_end=t_end,
    nt=nt,
    cfl=cfl,
    x_range=(xmin, xmax),
    y_range=(ymin, ymax),
    z_range=(zmin, zmax),
    deltas=(dx, dy, dz),
    em_params=em_params,
    particle_params=particle_params
)

create_particles(sim_params, electrons, particle_data)
create_particles(sim_params, ions, particle_data)
update_header(sim_params, project_path=project_path)
