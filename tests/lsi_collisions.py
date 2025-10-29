#!/usr/bin/env python3

# import numpy as np
from scipy import constants
import math

from particle_generation import create_particles
from domain_params import *

project_path = '/home/cepheid/TriForce/game_engine'
particle_data = project_path + '/data'


shape = (1501, 2, 1501)

xmin, xmax = -15.0e-6, 15.0e-6
zmin, zmax = -15.0e-6, 15.0e-6

dx = (xmax - xmin) / (shape[0] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dy = dx
ymin, ymax = 0.0, dx

t_end = 3e-13
dt = 4e-17
nt = int(t_end / dt) + 1
cfl = constants.c * dt * math.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)



particle_params = ParticleParams(
    particle_data=('electrons', 'ions'),
    collisions=(
        Collision(('electrons', 'ions'), coulomb_log=0.0, types=('coulomb',), step_interval=50, self_scatter=False),
        Collision(('electrons', 'electrons'), coulomb_log=0.0, types=('coulomb',), step_interval=50, self_scatter=True),
        Collision(('ions', 'ions'), coulomb_log=0.0, types=('coulomb',), step_interval=50, self_scatter=True),
    )
)

sim_params = Simulation(
    name='lsi_collisions',
    shape=shape,
    save_interval=100,
    nthreads=32,
    dt=dt,
    t_end=t_end,
    nt=nt,
    cfl=cfl,
    x_range=(xmin, xmax),
    y_range=(ymin, ymax),
    z_range=(zmin, zmax),
    deltas=(dx, dy, dz),
    particle_params=particle_params,
    push_enabled=True,
    jdep_enabled=True,
    em_enabled=True
)

# ===== Particles =====
px_range = (-5e-7, 5e-7)
py_range = (ymin, ymax)
pz_range = (-1e-5, 1e-5)

electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-1 * constants.e,
    temp=10000, # eV
    density=8.5e27, # m^-3,
    ppc=(10, 1, 10),
    distribution='relativistic',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

ions = Particles(
    name='ions',
    mass=constants.m_p,
    charge=1 * constants.e,
    temp=10000, # eV
    density=8.5e27, # m^-3,
    ppc=(10, 1, 10),
    distribution='relativistic',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

create_particles(sim_params, electrons, particle_data)
create_particles(sim_params, ions, particle_data)
update_header(sim_params, project_path=project_path)
