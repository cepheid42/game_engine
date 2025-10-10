#!/usr/bin/env python3

# import numpy as np
from scipy import constants
import math

from particle_generation import Particles, create_particles
from domain_params import Simulation, update_header

shape = (11, 11, 11)

xmin, xmax = 0.0, 0.001
ymin, ymax = 0.0, 0.001
zmin, zmax = 0.0, 0.001

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 1.829541541469147e-13
t_end = 7.318166165876587e-11
nt = int(t_end / dt) + 1
cfl = constants.c * dt * math.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

sim_params = Simulation(
    name='lsi',
    shape=shape,
    save_interval=10,
    nthreads=1,
    particle_bcs=1,
    interp_order=1,
    em_bcs=(2, 2, 2, 2, 2, 2),
    pml_depth=15,
    dt=dt,
    t_end=t_end,
    nt=nt,
    x_range=(xmin, xmax),
    y_range=(ymin, ymax),
    z_range=(zmin, zmax),
    deltas=(dx, dy, dz),
    cfl=cfl
)

# constexpr auto q_e = constants::q_e<double>;
# constexpr auto m_e = constants::m_e<double>;
# ParticleGroup g1("electrons", m_e, -q_e);
# constexpr vec3 loc0{5.25, 5.25, 5.25};
# constexpr vec3 vel{1.8e6, 0.0, 0.0};
# constexpr auto weight = 3.5;
# const Particle p0 = {loc0, loc0, vel, weight, calculateGamma(vel), false};
# g1.particles.push_back(p0);

update_header(sim_params)
