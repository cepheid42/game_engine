#!/usr/bin/env python3

# import numpy as np
from scipy import constants
import math

from scripts.particle_generation import create_particles
from scripts.domain_params import *

project_path = '/home/cepheid/TriForce/game_engine'
particle_data = project_path + '/data'


shape = (2, 2, 2)

xmin, xmax = 0.0, 1e-6
ymin, ymax = 0.0, 1e-6
zmin, zmax = 0.0, 1e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

t_end = 5e-12
dt = 5e-16
nt = int(t_end / dt) + 1
cfl = constants.c * dt * math.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

# ===== Particles =====
px_range = (xmin, xmax)
py_range = (ymin, ymax)
pz_range = (zmin, zmax)
ppc = (20, 20, 20)
Z_carbon = 6
m_carbon = 1.9945e-26 # kg
n_carbon = 1e26 # m^3
T_hot = 1250 / math.sqrt(3)
T_cold = 250 / math.sqrt(3)
lnL = 10
carbon1 = Particles(
    name='carbon1',
    mass=m_carbon,
    atomic_number=Z_carbon,
    charge=Z_carbon,
    temp=(T_hot, T_hot, T_hot), # eV
    density=n_carbon, # m^-3,
    ppc=ppc,
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

carbon2 = Particles(
    name='carbon2',
    mass=m_carbon,
    atomic_number=Z_carbon,
    charge=Z_carbon,
    temp=(T_cold, T_cold, T_cold), # eV
    density=n_carbon, # m^-3,
    ppc=ppc,
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

particle_params = ParticleParams(
    save_interval=nt//100,
    particle_bcs='periodic',
    interp_order=1,
    particle_data=(carbon1, carbon2),
    collisions=(
        Collision(groups=('carbon1', 'carbon2'), products=('carbon1', 'carbon2'), channels=('coulomb',), coulomb_log=lnL, self_scatter=False),
        Collision(groups=('carbon1', 'carbon1'), products=('carbon1', 'carbon1'), channels=('coulomb',), coulomb_log=lnL, self_scatter=True),
        Collision(groups=('carbon2', 'carbon2'), products=('carbon2', 'carbon2'), channels=('coulomb',), coulomb_log=lnL, self_scatter=True),
    )
)

sim_params = Simulation(
    name='carbon_thermal_eq',
    shape=shape,
    nthreads=1,
    dt=dt,
    t_end=t_end,
    nt=nt,
    cfl=cfl,
    x_range=(xmin, xmax),
    y_range=(ymin, ymax),
    z_range=(zmin, zmax),
    deltas=(dx, dy, dz),
    particle_params=particle_params,
    push_enabled=False,
    jdep_enabled=False,
    em_enabled=False
)

create_particles(sim_params, carbon1, particle_data)
create_particles(sim_params, carbon2, particle_data)

update_header(sim_params, project_path=project_path)
