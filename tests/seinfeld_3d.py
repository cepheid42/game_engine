#!/usr/bin/env python3

# import numpy as np
from scipy import constants
import math

from particle_generation import create_particles
from domain_params import *

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

em_params = EMParams(
    pml_depth=8,
    em_bcs=(1, 1, 2, 2, 1, 1),
)

particle_params = ParticleParams(
    particle_bcs=1,
    interp_order=2,
    particle_data=('electrons', 'ions')
)

sim_params = Simulation(
    name='seinfeld3D',
    shape=shape,
    save_interval=10,
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

# ===== Particles =====
px_range = (-0.04, 0.04)
pz_range = (-0.15, 0.15)


electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-constants.e,
    temp=4, # eV
    density=1.0e17, # m^-3,
    ppc=(3, 1, 3),
    distribution='relativistic',
    px_range=px_range,
    py_range=sim_params.y_range,
    pz_range=pz_range
)

ions = Particles(
    name='ions',
    mass=constants.m_p,
    charge=+constants.e,
    temp=1, # eV
    density=1.0e17, # m^-3,
    ppc=(2, 1, 2),
    distribution='relativistic',
    px_range=px_range,
    py_range=sim_params.y_range,
    pz_range=pz_range
)

create_particles(sim_params, electrons, particle_data)
create_particles(sim_params, ions, particle_data)
update_header(sim_params, project_path=project_path)
