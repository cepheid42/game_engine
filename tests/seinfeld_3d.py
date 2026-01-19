#!/usr/bin/env python3

# import numpy as np
from scipy import constants
import math

from scripts.particle_generation import create_particles
from scripts.domain_params import *

project_path = '/home/cepheid/TriForce/game_engine'
particle_data = project_path + '/data'

shape = (45, 2, 101)

xmin, xmax = -0.11, 0.11
zmin, zmax = -0.25, 0.25

dx = (xmax - xmin) / (shape[0] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

ymin, ymax = 0.0, 0.005
dy = (ymax - ymin) / (shape[1] - 1)

dt = 2.5e-12
t_end = 4000 * dt
nt = int(t_end / dt) + 1
cfl = constants.c * dt * math.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

# ===== Particles =====
px_range = (-0.04, 0.04)
py_range = (ymin, ymax)
pz_range = (-0.15, 0.15)

electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-1,
    atomic_number=0,
    temp=(4, 4, 4), # eV
    density=1.0e17, # m^-3,
    ppc=(3, 1, 3),
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
    temp=(1, 1, 1), # eV
    density=1.0e17, # m^-3,
    ppc=(2, 1, 2),
    distribution='relativistic',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# ===== Collisions and Particle Params =====
particle_params = ParticleParams(
    save_interval=10,
    particle_bcs='periodic',
    interp_order=2,
    particle_data=(electrons, ions)
)

# ===== Electromagnetic Params =====
em_params = EMParams(
    save_interval=10,
    pml_depth=8,
    em_bcs=(1, 1, 2, 2, 1, 1),
)

sim_params = Simulation(
    name='seinfeld3D',
    shape=shape,
    nthreads=24,
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
