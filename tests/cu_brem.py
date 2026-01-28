#!/usr/bin/env python3

# import numpy as np
from scipy import constants
import math

from scripts.domain_params import *
from scripts.particle_generation import create_particles

project_path = '/home/cepheid/TriForce/game_engine'
particle_data = project_path + '/data'

shape = (2, 2, 2)

xmin, xmax = 0.0, 1.6e-6
ymin, ymax = 0.0, 1.6e-6
zmin, zmax = 0.0, 1.6e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 1.0e-15
t_end = 4.2e-14
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
    temp=(0.0, 0.0, 40.0e6), # eV
    density=1.0e27, # m^-3,
    ppc=(100, 10, 10),
    distribution='relativistic',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

copper = Particles(
    name='copper',
    mass=63.55 * constants.atomic_mass,
    charge=0,
    atomic_number=29,
    temp=(1.0, 1.0, 1.0), # eV
    density=8.0e28, # m^-3,
    ppc=(100, 10, 10),
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

photons = Particles(
    name='photons',
    mass=0.0,
    charge=0.0,
    atomic_number=0,
    temp=(0.0, 0.0, 0.0), # eV
    density=0.0, # m^-3,
    ppc=(0, 0, 0),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# ===== Collisions and Particle Params =====
radiation_params = RadiationParams(
    products='photons',
    reduce_electron_energy=False,
    production_multiplier=1.0e5,
    cross_section_file='/tests/cross_section_data/SB_G4_Z29_kdsdk_eV_m2.txt'
)

particle_params = ParticleParams(
    save_interval=1,
    particle_bcs='periodic',
    interp_order=1,
    particle_data=(electrons, copper, photons),
    collisions=(
        Collision(
            groups=('electrons', 'copper'),
            channels=('radiation',),
            radiation=radiation_params,
        ),
    )
)

sim_params = Simulation(
    name='bremmstrahlung',
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

create_particles(sim_params, electrons, particle_data)
create_particles(sim_params, copper, particle_data)
update_header(sim_params, project_path=project_path)
