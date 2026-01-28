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

dt = 5.0e-14
t_end = 1.0e-12
nt = int(t_end / dt) + 1
cfl = constants.c * dt * math.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

# ===== Particles =====
px_range=(xmin, xmax)
py_range=(ymin, ymax)
pz_range=(zmin, zmax)
T_eV = 5.0e3

deuterium = Particles(
    name='deuterium',
    mass=2.0014 * constants.atomic_mass,
    charge=1,
    atomic_number=1,
    temp=(T_eV, T_eV, T_eV), # eV
    density=1.0e26, # m^-3,
    ppc=(25, 20, 20),
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

tritium = Particles(
    name='tritium',
    mass=3.01605 * constants.atomic_mass,
    charge=1,
    atomic_number=1,
    temp=(T_eV, T_eV, T_eV), # eV
    density=1.0e26, # m^-3,
    ppc=(25, 20, 20),
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

neutrons = Particles(
    name='neutrons',
    mass=constants.m_n,
    charge=0,
    atomic_number=0,
    temp=(0, 0, 0), # eV
    density=0.0, # m^-3,
    ppc=(0, 0, 0),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

helium3 = Particles(
    name='helium3',
    mass=3.016029 * constants.atomic_mass,
    charge=2,
    atomic_number=2,
    temp=(0, 0, 0), # eV
    density=0.0, # m^-3,
    ppc=(0, 0, 0),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

helium4 = Particles(
    name='helium4',
    mass=4.002603 * constants.atomic_mass,
    charge=2,
    atomic_number=2,
    temp=(0, 0, 0), # eV
    density=0.0, # m^-3,
    ppc=(0, 0, 0),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# ===== Collisions and Particle Params =====
DD_params = FusionParams(
    products=('neutrons', 'helium3'),
    energy_gain=3.269e6,
    production_multiplier=1.0e10,
    cross_section_file='/tests/cross_section_data/DD_nHe3_BH_eV_m2.txt'
)

DT_params = FusionParams(
    products=('neutrons', 'helium4'),
    energy_gain=17.589e6,
    production_multiplier=1.0e10,
    cross_section_file='/tests/cross_section_data/DT_nHe4_BH_eV_m2.txt'
)

particle_params = ParticleParams(
    save_interval=1,
    particle_bcs='periodic',
    interp_order=1,
    particle_data=(deuterium, neutrons, helium3),
    # particle_data=(tritium, deuterium, neutrons, helium4),
    collisions=(
        Collision(
            groups=('deuterium', 'deuterium'),
            channels=('fusion',),
            step_interval=1,
            self_scatter=True,
            fusion=DD_params
        ),
        # Collision(
        #     groups=('deuterium', 'tritium'),
        #     channels=('fusion',),
        #     step_interval=1,
        #     fusion=DT_params
        # ),
    )
)

sim_params = Simulation(
    name=f'DD_fusion_{int(T_eV * 1e-3)}keV',
    # name=f'DT_fusion_{int(T_eV * 1e-3)}keV',
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
    em_enabled=False,
    push_enabled=False,
    jdep_enabled=False
)

create_particles(sim_params, deuterium, particle_data)
create_particles(sim_params, tritium, particle_data)
update_header(sim_params, project_path=project_path)