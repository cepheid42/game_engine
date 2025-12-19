#!/usr/bin/env python3

# import numpy as np
from scipy import constants
import math

from scripts.domain_params import *
from scripts.particle_generation import create_particles

project_path = '/home/cepheid/TriForce/game_engine'
particle_data = project_path + '/data'

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
cfl = constants.c * dt * math.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

# ===== Particles =====
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
    ppc=(20, 20, 20),
    distribution='constant',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

neutral_aluminum = Particles(
    name='Al',
    mass=26.9815 * constants.atomic_mass + 13.0 * constants.m_e,
    charge=0,
    atomic_number=13,
    temp=(0, 0, 0), # eV
    density=6.6e28, # m^-3,
    ppc=(20, 20, 20),
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

ionized_aluminum = Particles(
    name='Al+',
    mass=26.9815 * constants.atomic_mass + 12.0 * constants.m_e,
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

em_params = EMParams(
    pml_depth=0,
    em_bcs=(2, 2, 2, 2, 2, 2),
)

particle_params = ParticleParams(
    save_interval=10,
    particle_bcs='periodic',
    interp_order=1,
    particle_data=(electrons, neutral_aluminum, ionized_aluminum),
    collisions=(
        Collision(
            groups=('electrons', 'Al'),
            products=('electrons', 'Al+'),
            channels=('ionization',),
            coulomb_log=10,
            self_scatter=False,
            ionization_energy=5.9858,
            cross_section_file='/data/al0cs.txt'
        ),
    )
)

sim_params = Simulation(
    name='ionization',
    shape=shape,
    nthreads=8,
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
    em_enabled=False,
    jdep_enabled=False
)

create_particles(sim_params, electrons, particle_data)
create_particles(sim_params, neutral_aluminum, particle_data)

update_header(sim_params, project_path=project_path)