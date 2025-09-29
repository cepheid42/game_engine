#!/usr/bin/env python3

# import numpy as np
from scipy import constants
from particle_generation import Particles, create_particles
from domain_params import Simulation, update_header

shape = (1501, 2, 1501)

xmin, xmax = -15.0e-6, 15.0e-6
zmin, zmax = -15.0e-6, 15.0e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = dx
dz = (zmax - zmin) / (shape[2] - 1)

ymin, ymax = 0.0, dx

dt = 4.0e-17
t_end = 3.0e-13
nt = int(t_end / dt) + 1

sim_params = Simulation(
    name='lsi',
    shape=shape,
    save_interval=100,
    nthreads=32,
    particle_bcs=1,
    interp_order=1,
    em_bcs=(1, 1, 2, 2, 1, 1),
    pml_depth=15,
    dt=dt,
    t_end=t_end,
    nt=nt,
    x_range=(xmin, xmax),
    y_range=(ymin, ymax),
    z_range=(zmin, zmax),
    deltas=(dx, dy, dz)
)

# ===== Particles =====
ppc = (10, 1, 10)
px_range = (-5e-7, 5e-7)
pz_range = (-1e-5, 1e-5)

electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-1.0 * constants.e,
    temp=10000, # eV
    density=8.5e27, # m^-3,
    ppc=ppc,
    distribution='relativistic',
    px_range=px_range,
    py_range=sim_params.y_range,
    pz_range=pz_range
)

ions = Particles(
    name='ions',
    mass=constants.m_p,
    charge=constants.e,
    temp=10000, # eV
    density=8.5e27, # m^-3,
    ppc=ppc,
    distribution='relativistic',
    px_range=px_range,
    py_range=sim_params.y_range,
    pz_range=pz_range
)

data_path = '/home/cepheid/TriForce/game_engine/data'
create_particles(sim_params, electrons, data_path)
create_particles(sim_params, ions, data_path)
update_header(sim_params)
